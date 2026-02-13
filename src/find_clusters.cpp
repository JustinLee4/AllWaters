#include "AtomicRadii.h"
#include "pdbtovector.h"
#include "map.h"
#include "internals.h"
#include "cluster.h"
#include "atom.h"

#include <iostream>
#include <iomanip>

double map_spacing = 0.25;

int main(int argc, char* argv[]){

    // 1. Check if the user provided exactly 3 arguments 
    // (argc is 4 because the program name itself is the 1st argument)
    if (argc != 3) {
        std::cerr << "Usage Error: Expected 2 arguments." << std::endl;
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <output_name>" << std::endl;
        return 1; // Exit with error code
    }

    // 2. Assign arguments to strings
    std::string input_file  = argv[1];
    std::string output_file = argv[2];

    // 3. (Optional) Basic verification - just for demonstration
    std::cout << "--- Job Started ---" << std::endl;
    std::cout << "Input PDB:  " << input_file << std::endl;
    std::cout << "Output:     " << output_file << std::endl;


    auto start_time = std::chrono::high_resolution_clock::now();

    std::tuple<std::vector<Atom>, double, double, double, double, double, double> newtuple = pdbtovector(input_file);

    // std::cout << "--- Tuple 1 (File 1) ---" << std::endl;
    // std::cout << "Atoms Found: " << std::get<0>(newtuple).size() << std::endl;
    // std::cout << "X Range: " << std::get<1>(newtuple) << " to " << std::get<2>(newtuple) << std::endl;
    // std::cout << "Y Range: " << std::get<3>(newtuple) << " to " << std::get<4>(newtuple) << std::endl;
    // std::cout << "Z Range: " << std::get<5>(newtuple) << " to " << std::get<6>(newtuple) << std::endl;


    std::vector<Atom> allAtoms = std::get<0>(newtuple);

    // 1. Cluster the atoms efficiently using the spatial grid
    std::cout << "-> Clustering " << allAtoms.size() << " points" << std::flush;
    std::vector<std::vector<Atom>> clusters = clusterAtoms(allAtoms, map_spacing);
    std::cout << "\n** Found " << clusters.size() << " clusters **" << std::endl;
    
    size_t totalClusteredAtoms = 0;
    for (const auto& cluster : clusters) {
        totalClusteredAtoms += cluster.size();
    }

    std::cout << "-> Total Clustered Atoms: " << totalClusteredAtoms << "  "
              << (totalClusteredAtoms == allAtoms.size() ? "PASS" : "FAIL") << std::endl;

    // 2. Write the clustered PDB
    std::cout << "-> Writing to " << output_file << ".pdb" << std::flush;
    writeClusteredPDB(clusters, output_file);
    std::cout << "\n-> Success" << std::endl;



    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    long long total_seconds = static_cast<long long>(elapsed.count());
    long long hours = total_seconds / 3600;
    long long minutes = (total_seconds % 3600) / 60;
    long long seconds = total_seconds % 60;
    std::cout << "Time: " 
          << std::setfill('0') << std::setw(2) << hours << ":"
          << std::setfill('0') << std::setw(2) << minutes << ":"
          << std::setfill('0') << std::setw(2) << seconds 
          << std::endl;

    return 0;
}