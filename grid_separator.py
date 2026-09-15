import numpy as np
from scipy.spatial import cKDTree
import sys


def parse_pdb_coords(pdb_path):
    """
    Parses a PDB file and returns:
    1. A numpy array of (N, 3) coordinates.
    2. A list of the actual string lines from the file (to save them later).
    """
    coords = []
    lines = []

    print(f"Reading PDB: {pdb_path}...")
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        # PDB Fixed Width Parsing
                        # x: 30-38, y: 38-46, z: 46-54
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                        lines.append(line)
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"Error: File {pdb_path} not found.")
        sys.exit(1)

    return np.array(coords), lines


def load_msms_verts(vert_path):
    """
    Parses MSMS .vert file.
    Assumes format: x y z nx ny nz (ignoring extra cols).
    """
    print(f"Reading MSMS vertices: {vert_path}...")
    try:
        # Load data, ignoring lines starting with '#' (headers)
        data = np.loadtxt(vert_path, skiprows=3, usecols=(0, 1, 2, 3, 4, 5))
    except Exception as e:
        print(f"Error loading MSMS file: {e}")
        sys.exit(1)

    vertices = data[:, 0:3]
    normals = data[:, 3:6]
    return vertices, normals


def inject_bfactor(line, value):
    """
    Replaces the B-factor column (chars 60-66) in a PDB line
    with the signed distance value.
    """
    # Format as float with 2 decimal places, width 6
    b_factor_str = f"{value:6.2f}"

    # Handle overflow if distance is > 999 or < -99 to maintain column width
    if len(b_factor_str) > 6:
        # If it overflows, just cap it visually so the PDB doesn't break
        b_factor_str = "999.99" if value > 0 else "-99.99"

    # Reconstruct the line:
    # line[:60] is everything before B-factor
    # line[66:] is everything after (occupancy, etc.)
    return line[:60] + b_factor_str + line[66:]


def main(
    pdb_input,
    msms_vert_input,
    out_inside="inside_colored.pdb",
    out_outside="outside_colored.pdb",
):
    # --- 1. Load Data ---
    grid_points, grid_lines = parse_pdb_coords(pdb_input)
    vertices, normals = load_msms_verts(msms_vert_input)

    if len(grid_points) == 0:
        print("Error: No coordinates found in PDB file.")
        return

    print(f"Loaded {len(grid_points)} grid points.")
    print(f"Loaded {len(vertices)} surface vertices.")

    # --- 2. Build KD-Tree ---
    print("Building KD-Tree...")
    tree = cKDTree(vertices)

    # --- 3. Find Nearest Neighbors ---
    print("Querying nearest vertices...")
    # k=1 returns (distances, indices)
    distances, indices = tree.query(grid_points, k=1, workers=-1)

    # --- 4. Calculate Signed Distance ---
    print("Calculating signs...")

    nearest_verts = vertices[indices]
    nearest_norms = normals[indices]

    # Vector from surface vertex to grid point
    diff_vectors = grid_points - nearest_verts

    # Dot product:
    # dot > 0 : Outside (positive distance)
    # dot < 0 : Inside (negative distance)
    dot_products = np.sum(diff_vectors * nearest_norms, axis=1)

    # Determine sign (-1 or 1)
    signs = np.sign(dot_products)
    signs[signs == 0] = 1  # Treat exactly on surface as outside/positive

    # Signed Distance = Distance * Sign
    signed_distances = distances * signs

    # --- 5. Categorize & Write ---
    print("Writing output files with Signed Distance in B-factor column...")

    count_in = 0
    count_out = 0

    with open(out_inside, "w") as f_in, open(out_outside, "w") as f_out:
        for i, line in enumerate(grid_lines):
            sd = signed_distances[i]

            # Inject the signed distance into the B-factor column
            new_line = inject_bfactor(line, sd)

            if sd < 0:
                f_in.write(new_line)
                count_in += 1
            else:
                f_out.write(new_line)
                count_out += 1

    print("Done.")
    print(f"  -> Inside:  {count_in} points (saved to {out_inside})")
    print(f"  -> Outside: {count_out} points (saved to {out_outside})")


if __name__ == "__main__":
    # --- CONFIGURATION ---
    # Replace these with your actual filenames
    grid_pdb = sys.argv[1]
    msms_vert = sys.argv[2]

    main(grid_pdb, msms_vert)
