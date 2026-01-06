import sys
import radii as source_data

def generate_cpp():
    # Write directly to the file
    with open("AtomicRadii.cpp", "w") as f:
        
        # Redirect standard output to the file
        original_stdout = sys.stdout
        sys.stdout = f
        
        try:
            print('#include "AtomicRadii.h"')
            print("")
            print("const std::map<AtomKey, AtomParams> ATOMIC_RADII = {")
            
            for (res, atom), values in source_data.radii.items():
                # 1. Handle Residue Name (Strip whitespace, convert None to "")
                if res is None:
                    res_str = '""'
                else:
                    res_str = f'"{res.strip()}"' # .strip() removes spaces
                
                # 2. Handle Atom Name (Strip whitespace)
                atom_str = f'"{atom.strip()}"'   # .strip() removes spaces
                
                # 3. Format the values
                params = f"{{ {values[0]}, {values[1]}, {values[2]}, {values[3]}, {values[4]}, {values[5]} }}"
                
                # 4. Print the C++ map entry
                print(f'    {{ {{{res_str}, {atom_str}}}, {params} }},')

            print("};")
            
        finally:
            sys.stdout = original_stdout

    print("Successfully generated AtomicRadii.cpp with stripped strings.")

if __name__ == "__main__":
    generate_cpp()