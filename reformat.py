import sys
import os
import argparse

# --- CONFIGURATION ---
# 1. ONLY Water Oxygens
TARGET_ATOMS = {'O', 'OW', 'OH2'} 

# 2. ONLY Water Residues (Common names in PDBs)
TARGET_RESIDUES = {'HOH', 'WAT', 'SOL', 'TIP3', 'TIP4', 'SPC'}
# ---------------------

def parse_pdb_line_smart(line):
    """
    Parses a PDB line (even with broken/merged columns) and extracts
    data ONLY if it is a Water Oxygen.
    """
    parts = line.split()
    
    # Safety check: Line too short to be an atom?
    if len(parts) < 8:
        return None

    try:
        # --- STRATEGY 1: Parse Coordinates from the END (Right-to-Left) ---
        # This avoids issues where "A30010" shifts columns on the left.
        # Standard PDB lines end with: ... X Y Z Occ Temp Element
        # We grab X, Y, Z relative to the end of the line.
        
        # indices: ... [X] [Y] [Z] [Occ] [Temp] [Element?]
        # If element column exists (is a string), coords are at -6, -5, -4
        # If element column is missing (is a float), coords are at -5, -4, -3
        
        # Check if last item is a number (Temp) or String (Element)
        try:
            float(parts[-1]) 
            has_element_col = False
        except ValueError:
            has_element_col = True
            
        if has_element_col:
            z_idx, y_idx, x_idx = -4, -5, -6
        else:
            z_idx, y_idx, x_idx = -3, -4, -5
            
        z = float(parts[z_idx])
        y = float(parts[y_idx])
        x = float(parts[x_idx])

        # --- STRATEGY 2: Identify Atom and Residue Names ---
        # We need to handle merged columns like "HETATM30010" vs "HETATM 30010"
        
        # Case A: Merged Record+Serial (e.g., ['HETATM30010', 'O', 'HOH', ...])
        if len(parts[0]) > 6: 
            atom_name = parts[1]
            res_name = parts[2]
            
        # Case B: Split Record+Serial (e.g., ['HETATM', '30010', 'O', 'HOH', ...])
        else:
            atom_name = parts[2]
            res_name = parts[3]

        # --- FILTERING ---
        # Strictly check if this is a Water Oxygen
        if atom_name not in TARGET_ATOMS:
            return None
        if res_name not in TARGET_RESIDUES:
            return None

        return {
            'atom_name': atom_name,
            'res_name': res_name, 
            'coords': (x, y, z)
        }
        
    except (ValueError, IndexError):
        # Determine if line was just bad data or a header/footer
        return None

def format_pdb_line(atom_data, serial, res_seq):
    """
    Writes a CLEAN PDB line adhering to strict column standards.
    - Col 17 (Alt Loc) is left blank.
    - Col 18-20 is Residue Name.
    - Col 22 is Chain ID.
    """
    x, y, z = atom_data['coords']
    name = atom_data['atom_name']
    res_name = atom_data['res_name']
    
    # 1. Handle Residue Sequence Wrapping (1-9999)
    # We use 9999 (standard max) but ensure it doesn't merge with Chain ID
    safe_res_seq = ((res_seq - 1) % 9999) + 1
    
    # 2. Format Atom Name (Center/Left alignment rules)
    # Standard PDB: " O  " (Space, O, Space, Space) for Oxygen
    if len(name) == 1: fmt_name = f" {name}  "
    elif len(name) == 2: fmt_name = f" {name} "
    elif len(name) == 3: fmt_name = f" {name}"
    else: fmt_name = name[:4] # Truncate if >4 to prevent alignment break

    # 3. Format Residue Name
    # Ensure it is exactly 3 chars. If "TIP3" (4 chars), PDBs often use "TP3" or "TIP".
    # We take the first 3 chars to be safe.
    fmt_res = res_name[:3]

    line = (
        f"HETATM"             # Cols 1-6
        f"{serial:>5}"        # Cols 7-11 (Serial)
        f" "                  # Col 12    (Blank)
        f"{fmt_name:<4}"      # Cols 13-16 (Atom Name)
        f" "                  # Col 17    (Alt Loc - FIX: Added this space)
        f"{fmt_res:>3}"       # Cols 18-20 (Res Name)
        f" "                  # Col 21    (Blank)
        f"A"                  # Col 22    (Chain ID)
        f"{safe_res_seq:>4}"  # Cols 23-26 (Res Seq)
        f"    "               # Cols 27-30 (Code)
        f"{x:>8.3f}"          # Cols 31-38 (X)
        f"{y:>8.3f}"          # Cols 39-46 (Y)
        f"{z:>8.3f}"          # Cols 47-54 (Z)
        f"{1.00:>6.2f}"       # Cols 55-60 (Occ)
        f"{0.00:>6.2f}"       # Cols 61-66 (Temp)
        f"           "        # Spacing
        f"{name[0]:>1}"       # Element
        f"\n"
    )
    return line

def process_files(input_files, output_path):
    all_atoms = []
    print(f"Scanning {len(input_files)} files for Water Oxygens...")

    for fpath in input_files:
        if not os.path.exists(fpath):
            continue
            
        print(f"  -> Reading {os.path.basename(fpath)}...")
        count = 0
        try:
            with open(fpath, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        atom = parse_pdb_line_smart(line)
                        if atom:
                            all_atoms.append(atom)
                            count += 1
            print(f"     Found {count} water oxygens.")
        except Exception as e:
            print(f"Error reading {fpath}: {e}")

    # Write Output
    print(f"\nWriting {len(all_atoms)} combined water oxygens to {output_path}...")
    
    with open(output_path, 'w') as out:
        for i, atom in enumerate(all_atoms):
            # Serial: 1 to N
            # ResSeq: 1 to N (wrapped by formatter)
            serial = i + 1
            out.write(format_pdb_line(atom, serial, serial))
        out.write("END\n")
    print("Done.")

def main():
    parser = argparse.ArgumentParser(description="Combine PDBs (Water Oxygens Only).")
    parser.add_argument('inputs', nargs='+', help='List of input PDB files')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument('-n', '--name', default='reformat_waters.pdb', help='Output filename')
    
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    output_file = os.path.join(args.outdir, args.name)
    process_files(args.inputs, output_file)

if __name__ == "__main__":
    main()