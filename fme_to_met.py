import sys

def convert_fme_to_met(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            # Check if this line belongs to FME
            if "FME" in line:
                # 1. Extract the Atom Name (columns 12-16)
                atom_name = line[12:16].strip()
                
                # 2. DELETE the Formyl cap atoms (CN and O1)
                if atom_name == "CN" or atom_name == "O1":
                    continue
                
                # 3. REPLACE 'HETATM' with 'ATOM  ' (keep spacing)
                # PDB 'ATOM' is usually followed by 2 spaces. HETATM is 6 chars.
                line = line.replace("HETATM", "ATOM  ")
                
                # 4. REPLACE 'FME' with 'MET'
                line = line.replace("FME", "MET")
                
                f_out.write(line)
            else:
                # Write all other lines (standard residues) unchanged
                f_out.write(line)

    print(f"Done! Saved converted file to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fme_to_met.py input.pdb output.pdb")
    else:
        convert_fme_to_met(sys.argv[1], sys.argv[2])