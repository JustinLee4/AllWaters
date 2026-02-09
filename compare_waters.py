import sys
import math
import argparse
import os

# --- GLOBAL CONFIGURATION ---
OUTPUT_DIR = "comparison"
# ----------------------------

def parse_pdb_points(filename):
    """
    Parses a PDB file using STRICT column slicing.
    Handles 'HETATM10000' merged columns correctly.
    """
    points = []
    # Target atom names (stripped)
    target_atoms = {'O', 'OW', 'XP'}
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                # 1. Check Record Type (Cols 1-6)
                #    We strictly check the first 6 chars. 
                #    This handles HETATM10000 because we don't look past index 6 yet.
                record_type = line[0:6].strip()
                
                if record_type == "ATOM" or record_type == "HETATM":
                    try:
                        # 2. Extract Atom Name (Cols 13-16 -> Index 12-16)
                        atom_name = line[12:16].strip()
                        
                        if atom_name in target_atoms:
                            # 3. Extract Serial Number (Cols 7-11 -> Index 6-11)
                            #    This is the key fix. Even if the line is "HETATM10000",
                            #    line[6:11] will correctly grab "10000".
                            serial_str = line[6:11].strip()
                            
                            # 4. Extract Coordinates (Standard Fixed Width)
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            
                            # 5. Extract Metadata
                            res_name = line[17:20].strip()
                            res_seq_str = line[22:26].strip()
                            
                            # Handle potentially empty/merged serials
                            serial = int(serial_str) if serial_str else 0
                            res_seq = int(res_seq_str) if res_seq_str else 0

                            points.append({
                                'coords': (x, y, z),
                                'atom_name': atom_name,
                                'res_name': res_name,
                                'res_seq': res_seq,
                                'serial': serial,
                                'matched': False
                            })
                            
                    except ValueError:
                        # Skip lines that look like atoms but have bad numbers
                        continue
                        
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        sys.exit(1)
        
    print(f"-> Parsed {len(points)} atoms from {filename}")
    return points

def get_grid_key(coords, cell_size=1.0):
    return (
        int(math.floor(coords[0] / cell_size)),
        int(math.floor(coords[1] / cell_size)),
        int(math.floor(coords[2] / cell_size))
    )

def build_spatial_grid(points, cell_size=1.0):
    grid = {}
    for p in points:
        key = get_grid_key(p['coords'], cell_size)
        if key not in grid:
            grid[key] = []
        grid[key].append(p)
    return grid

def get_neighbor_keys(center_key):
    cx, cy, cz = center_key
    neighbors = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                neighbors.append((cx + dx, cy + dy, cz + dz))
    return neighbors

def distance_sq(c1, c2):
    return (c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2

def format_pdb_line(atom, serial_override=None):
    """
    Formats atom data into the STRICT PDB standard.
    Matches the 'vectortopdb' C++ format exactly.
    """
    x, y, z = atom['coords']
    
    # Use override serial if provided, otherwise use atom's original serial
    serial = serial_override if serial_override is not None else atom['serial']
    res_seq = serial # Usually for waters, resSeq matches serial
    
    # Handle Atom Name Alignment (O must be " O  ")
    raw_name = atom['atom_name']
    if len(raw_name) == 1:
        fmt_name = f" {raw_name}  " # Column 14
    elif len(raw_name) == 2:
        fmt_name = f" {raw_name} "  # Column 14
    elif len(raw_name) == 3:
        fmt_name = f" {raw_name}"   # Column 14
    else:
        fmt_name = raw_name         # 4 chars fit perfectly
        
    # Python f-string formatting for fixed widths:
    # <  : Left align
    # >  : Right align
    # 8.3f : Float with 3 decimal places, width 8
    
    line = (
        f"HETATM"
        f"{serial:>5}"        # Serial (6-11)
        f" "                  # (12)
        f"{fmt_name:<4}"      # Name (13-16)
        f"{atom['res_name']:>3}" # ResName (17-20)
        f"  "                 # (21)
        f"A"                  # ChainID (22) -> Forced to 'A'
        f"{res_seq:>4}"       # ResSeq (23-26)
        f"    "               # (27-30)
        f"{x:>8.3f}"          # X (31-38)
        f"{y:>8.3f}"          # Y (38-46)
        f"{z:>8.3f}"          # Z (46-54)
        f"{1.00:>6.2f}"       # Occ (55-60)
        f"{0.00:>6.2f}"       # Temp (61-66)
        f"           "        # (67-77)
        f"{raw_name[0]:>1}"   # Element (78)
        f"\n"
    )
    return line

def write_pdb(filename, points):
    """Writes list of atoms using strict formatting."""
    try:
        with open(filename, 'w') as f:
            # f.write("REMARK Generated by Python Script\n")
            for i, p in enumerate(points):
                # Clean serial numbers (1, 2, 3...)
                f.write(format_pdb_line(p, serial_override=i+1))
            f.write("END\n")
        print(f"-> Wrote {len(points)} atoms to {filename}")
    except IOError as e:
        print(f"Error writing to {filename}: {e}")

def generate_output_name(input_path, output_dir, prefix="unique_to_"):
    filename = os.path.basename(input_path)
    name, ext = os.path.splitext(filename)
    new_filename = f"{prefix}{name}{ext}"
    return os.path.join(output_dir, new_filename)

def main():
    parser = argparse.ArgumentParser(description="Compare coordinates in two PDB files.")
    parser.add_argument("file1", help="Input PDB File 1")
    parser.add_argument("file2", help="Input PDB File 2")
    parser.add_argument("--tolerance", type=float, default=0.01, help="Tolerance (A)")
    
    args = parser.parse_args()

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    print(f"Reading {args.file1}...")
    points1 = parse_pdb_points(args.file1)
    
    print(f"Reading {args.file2}...")
    points2 = parse_pdb_points(args.file2)

    print(f"Comparing {len(points1)} vs {len(points2)} atoms...")

    # Optimization: Grid for File 2
    grid_size = 3.5 # Optimized for water search (matches C++ logic)
    grid2 = build_spatial_grid(points2, grid_size)
    
    common_points = []
    unique_file1 = []
    
    tol_sq = args.tolerance ** 2
    
    for p1 in points1:
        key = get_grid_key(p1['coords'], grid_size)
        possible_keys = get_neighbor_keys(key)
        
        match_found = False
        
        for n_key in possible_keys:
            if n_key in grid2:
                for p2 in grid2[n_key]:
                    if distance_sq(p1['coords'], p2['coords']) <= tol_sq:
                        p2['matched'] = True
                        match_found = True
                        common_points.append(p1)
                        break 
            if match_found:
                break
        
        if not match_found:
            unique_file1.append(p1)

    unique_file2 = [p for p in points2 if not p['matched']]

    out_name_1 = generate_output_name(args.file1, OUTPUT_DIR, prefix="unique_to_")
    out_name_2 = generate_output_name(args.file2, OUTPUT_DIR, prefix="unique_to_")
    out_name_common = os.path.join(OUTPUT_DIR, "common_points.pdb")

    print("\n--- Results ---")
    write_pdb(out_name_1, unique_file1)
    write_pdb(out_name_2, unique_file2)
    write_pdb(out_name_common, common_points)

if __name__ == "__main__":
    main()