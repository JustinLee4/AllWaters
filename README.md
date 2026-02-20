# AllWaters

The workflow is as such:

1. from /pdbfiles, run fme_to_met and other necessary helpers to build a valid pdb file
2. Run convert.py on an Atomic Radii file to convert to a C++ Map
3. In the command line, run 
    $ make
    to build the allwaters program. The program should be located in bin/allwaters.exe
4. To call allwaters, run from parent folder allwaters
    $ bin/allwaters.exe -p <path to input pdb file> -v <path to .vert file> -o <path to output file (no .pdb extension required)>
    optional flags include:
        -r <value> (alternatively --radius, default 3.5 A)
            this determines which gridpoints are categorized as surface/internal - all gridpoints within <-r> Angstroms of the surface are categorized as surface
        -s <value> (alternatively --spacing, default 0.25 A)
            this determines the grid spacing of generated water molecules.
        -cluster
            including this flag skips all generation and assumes that you provide a pdb file of only waters for clustering purposes - this will use the <-s> argument to calculate neighbors
5. Results should be found in results/

compare_waters.py is included as a script to compare with other pdb files.