#include "AtomicRadii.h"
#include "cluster.h"
#include "common.h"
#include "internals.h"
#include "pdbtovector.h"
#include "map.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iomanip>


double hash_spacing = 3;
double grid_spacing = .25;
double water_diameter = 2.5;
double cutoff_distance = 5;
float shellradius = 2.5;


int main(int argc, char* argv[]){

    auto start_time = std::chrono::high_resolution_clock::now();

// ---------- input handling ----------
    std::string input_file = "";
    std::string vert_file = "";
    std::string output_file = "";
    std::string r_value = "3.5";
    bool only_cluster = false;


    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "-p" || arg == "--pdb") && i + 1 < argc) {
            input_file = argv[++i]; 
        } 
        else if ((arg == "-v" || arg == "--vert") && i + 1 < argc) {
            vert_file = argv[++i];
        } 
        else if ((arg == "-o" || arg == "--out") && i + 1 < argc) {
            output_file = std::string("results/") + argv[++i];
        }         
        else if ((arg == "-r" || arg == "--radius") && i + 1 < argc) {
            r_value = argv[++i];             
            try {
                double test_value = std::stod(r_value);
            } 
            catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid radius. '" << r_value << "' is not a valid number." << std::endl;
                return 1;
            } 
            catch (const std::out_of_range& e) {
                std::cerr << "Error: Radius '" << r_value << "' is too large/small for a double." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-s" || arg == "--spacing") && i + 1 < argc) {
            std::string test_grid_spacing = argv[++i];             
            try {
                grid_spacing = std::stod(test_grid_spacing);
            } 
            catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid grid spacing. '" << test_grid_spacing << "' is not a valid number." << std::endl;
                return 1;
            } 
            catch (const std::out_of_range& e) {
                std::cerr << "Error: Spacing '" << test_grid_spacing << "' is too large/small for a double." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-cluster")) {
            only_cluster = true;
        }
    

        else {
            std::cerr << "Error: Unknown or incomplete argument '" << arg << "'" << std::endl;
            std::cerr << "Usage: " << argv[0] << " -p <pdb> -v <vert> -o <out> [-r <value>] [-cluster]" << std::endl;
            return 1;
        }
    }

    if ((input_file.empty() || vert_file.empty() || output_file.empty())) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        std::cerr << "Usage: " << argv[0] << " -p <pdb> -v <vert> -o <out> [-r <value>] [-cluster]" << std::endl;
        return 1;
    }


    std::cout << "--- Files ---" << std::endl;
    std::cout << "Input PDB:  " << input_file << std::endl;
    std::cout << "Input Vert: " << vert_file << std::endl;
    std::cout << "Output:     " << output_file << std::endl;
    std::cout << "--- Params ---" << std::endl;
    std::cout << "Grid Spacing: " << grid_spacing << std::endl;
    if(!only_cluster) {
        std::cout << "Water Diameter: " << water_diameter << std::endl;
        std::cout << "Internal/External Shell Radius (Initial/Flood Fill): " << shellradius << std::endl;
        std::cout << "Internal/External Shell Radius (Secondary/Categorize): " << r_value << std::endl;
    }


    std::string response;
    while (true) {
        std::cout << "\nProceed? (y/n): ";
        std::cin >> response;
        if (response == "y" || response == "Y") {
            break;
        } else if (response == "n" || response == "N") {
            std::cout << "Job cancelled by user." << std::endl;
            return 0;
        } else {
            std::cout << "Invalid input. Please enter 'y' or 'n'." << std::endl;
        }
    }

    if(!only_cluster) { 
    // ---------- reading PDB into vector ----------

        std::cout << "-> Entering pdbtovector" << std::endl;

        std::tuple<std::vector<Atom>, double, double, double, double, double, double> newtuple = pdbtovector(input_file);

        std::vector<Atom> atomvector = std::get<0>(newtuple);

        double minx = std::get<1>(newtuple);
        double maxx = std::get<2>(newtuple);
        double miny = std::get<3>(newtuple);
        double maxy = std::get<4>(newtuple);
        double minz = std::get<5>(newtuple);
        double maxz = std::get<6>(newtuple);

    //---------- generating gridpoints ---------


    // ---------- build hashmap with vector ----------

        std::cout << "-> Building hashmap" << std::endl;

        std::unordered_map<GridKey, std::vector<int>> map = buildSpatialGrid(atomvector, hash_spacing);

        double start_x = std::floor(minx);
        double end_x   = std::ceil(maxx);

        double start_y = std::floor(miny);
        double end_y   = std::ceil(maxy);

        double start_z = std::floor(minz);
        double end_z   = std::ceil(maxz);

        int nx = (int)((end_x - start_x) / grid_spacing + 1e-5) + 1;
        int ny = (int)((end_y - start_y) / grid_spacing + 1e-5) + 1;
        int nz = (int)((end_z - start_z) / grid_spacing + 1e-5) + 1;


        long total_reps = (long)((nx+10) * (ny+10) * (nz+10));

    // ---------- separate surface and internal ----------

        std::cout << "-> Processing vertices" << std::endl;


        std::vector<Vertex> mySurface = vert_to_vector(vert_file);

        std::vector<Vec3> insidePoints;
        std::vector<Vec3> outsidePoints;
        
        Vec3 minB = {(float)start_x - 5, (float)start_y - 5, (float)start_z - 5};
        Vec3 maxB = {(float)end_x + 5, (float)end_y + 5, (float)end_z + 5};
        
        SeparateGridPoints(mySurface, minB, maxB, (float)grid_spacing, shellradius, insidePoints, outsidePoints);

        
        PRINT_LOG(WriteWaterPDB(insidePoints, output_file + "_in.pdb"));
        PRINT_LOG(WriteWaterPDB(outsidePoints, output_file + "_out.pdb"));


        std::cout << "-> Flood fill" << std::endl;

        std::vector<Vec3> allpoints;
        std::vector<std::vector<Vec3>> layers;

        FillInternalVoid(outsidePoints, insidePoints, minB, maxB, .25, allpoints, layers);

        std::cout  << "-> total gridpoints = " << std::scientific << std::setprecision(3) << (double)total_reps << std::endl;
        std::cout << std::scientific << std::setprecision(3) << "-- removed " << (double)total_reps - (double)allpoints.size() << " (" << std::fixed << std::setprecision(1) << ((double)total_reps - (double)allpoints.size()) / total_reps * 100 <<"%) --" << std::endl;
        std::cout << std::scientific << std::setprecision(3) << "-> Dowsing " << (double)allpoints.size() << std::endl;
        
        PRINT_LOG(WriteWaterPDB(allpoints, output_file + "_all_internals_before_protein_overlap.pdb"));


    // ---------- remove overlaps with Protein atoms ----------

        std::vector<Atom> watervector = {};
        int dbug = 0;
        std::cout << "\033[?25l";
        int total_points = allpoints.size();

        for (int i = 0; i < total_points; i++){
            

            double percent = 0.0;
            if (i % 1000 == 0 || i == total_points - 1) {
                double percent = ((double)(i + 1) / total_points) * 100.0;
                std::cout << "\r** Iteration: " << i + 1 << " of " << total_points << " **\033[K\n"
                          << "   Progress:  " << std::fixed << std::setprecision(1) << percent << "%\033[K" << std::flush;
                std::cout << "\033[1A";
            }

            std::array<double, 3> temp_array = {allpoints[i].x, allpoints[i].y ,allpoints[i].z};
            if(getOverlap_cluster(map, atomvector, temp_array, hash_spacing, water_diameter, cutoff_distance)){
                Atom newatom("HOH", "O", temp_array);
                watervector.push_back(newatom);
            }
            dbug += 1;

        }
        std::cout << "\n\n\033[?25h";

        std::cout << "\n** There are " << watervector.size() << " waters **" <<std::endl;


        std::cout << "-> Entering vectortopdb" << std::endl;
        
        vectortopdb(watervector, output_file + "_all_internal_gridpoints.pdb");

    

    // ---------- run categorize_water ----------
    // ---------- this further separates gridpoints into internal/external based on proximity to surface ----------

    std::cout << "-> Entering categorize_water" << std::endl;

    std::string categorize_water_python = "python categorize_water/categorize_water.py -f " + output_file + "_all_internal_gridpoints.pdb" 
                                                                        " -s " + vert_file + 
                                                                        " -o " + output_file +
                                                                        " -r " + r_value;
    int result = std::system(categorize_water_python.c_str());

    if (result == 0) {
    } else {
        std::cerr << "categorize_water.py failed." << std::endl;
        return 1;
    }

    std::string reformat_python = "python scripts/reformat.py -o "+ output_file + "_reformatted.pdb " + output_file + "_internal.pdb " + output_file + "_surface.pdb";

    result = std::system(reformat_python.c_str());

    if (result == 0) {
    } else {
        std::cerr << "reformat.py failed." << std::endl;
        return 1;
    }

    input_file = output_file + "_reformatted.pdb";
    } 
    //---------- begin clustering ----------

    std::tuple<std::vector<Atom>, double, double, double, double, double, double> cluster_tuple = pdbtovector(input_file);

    std::vector<std::string> remarks = extractRemarks(input_file);
    std::string grid_spacing_remark = "grid spacing = " + std::to_string(grid_spacing);
    remarks.push_back(grid_spacing_remark);
    if(!only_cluster){
        std::string radius = "surface +- " + r_value;
        std::string water_diameter_remark = "water diameter = " + std::to_string(water_diameter);
        std::string vert_file_remark = "vert file = " + vert_file;
        remarks.push_back(radius);  remarks.push_back(water_diameter_remark); remarks.push_back(vert_file_remark);
    }

    std::vector<Atom> allAtoms = std::get<0>(cluster_tuple);

    std::cout << "-> Clustering " << allAtoms.size() << " points" << std::flush;
    std::vector<std::vector<Atom>> clusters = clusterAtoms(allAtoms, grid_spacing, hash_spacing);
    std::cout << "\n** Found " << clusters.size() << " clusters **" << std::endl;
    
    size_t totalClusteredAtoms = 0;
    for (const auto& cluster : clusters) {
        totalClusteredAtoms += cluster.size();
    }

    std::cout << "-> Total Clustered Atoms: " << totalClusteredAtoms << "  "
              << (totalClusteredAtoms == allAtoms.size() ? "PASS" : "FAIL") << std::endl;

    std::cout << "-> Writing to " << output_file << ".pdb" << std::flush;
    writeClusteredPDB(clusters, output_file, remarks);

    std::string internals_file = output_file + "_all_internal_gridpoints.pdb";
    std::string internals_categorized_file = output_file + "_internal.pdb";
    std::string surface_categorized_file = output_file + "_surface.pdb";
    std::string reformatted_file = output_file + "_reformatted.pdb";



    if(std::filesystem::exists(internals_file)) {
        std::filesystem::remove(internals_file);
    }
    if(std::filesystem::exists(internals_categorized_file)) {
        std::filesystem::remove(internals_categorized_file);
    }
    if(std::filesystem::exists(surface_categorized_file)) {
        std::filesystem::remove(surface_categorized_file);
    }
    if(std::filesystem::exists(reformatted_file)) {
        std::filesystem::remove(reformatted_file);
    }


    std::cout << "\n-> Success" << std::endl;

    //--------- timer ----------
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

    return 1;
}