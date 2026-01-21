#include "AtomicRadii.h"
#include "pdbtovector.h"
#include "map.h"
#include "internals.h"

#include <iostream>
#include <iomanip>


double hash_spacing = 3;
double grid_spacing = .25;
double water_diameter = 2.5;
double cutoff_distance = 5;

std::string input_file = "pdbfiles/L-sub_met.pdb";
// std::string input_file = "pdbfiles/test_sphere.pdb";

std::string vert_file = "vert_files/L-sub.vert";
// std::string vert_file = "vert_files/test_sphere.vert";

std::string output_file = "results/L-sub";
// std::string output_file = "results/test_sphere";



int main(int argc, char* argv[]){

    auto start_time = std::chrono::high_resolution_clock::now();


    // AtomParams params = getParams("ALA", "N");
    // std::cout << params.radius_ua << std::endl;

    std::cout << "-> Entering pdbtovector" << std::endl;

    std::tuple<std::vector<Atom>, double, double, double, double, double, double> newtuple = pdbtovector(input_file);

    std::vector<Atom> atomvector = std::get<0>(newtuple);

    double minx = std::get<1>(newtuple);
    double maxx = std::get<2>(newtuple);
    double miny = std::get<3>(newtuple);
    double maxy = std::get<4>(newtuple);
    double minz = std::get<5>(newtuple);
    double maxz = std::get<6>(newtuple);

    // std::cout << minx << maxx << miny << maxy << minz << maxz << std::endl;
    // std::cout << "atomvector size : " << atomvector.size() << std::endl;

    std::cout << "-> Building hashmap" << std::endl;

    std::unordered_map<GridKey, std::vector<int>> map = buildSpatialGrid(atomvector, hash_spacing);

    // printSpatialGrid(map);

    // std::cout << "-> Entering toString" << std::endl;

    // std::cout <<     atomvector[0].toString();

    // double start_x = std::floor(minx - 5);
    // double end_x   = std::ceil(maxx + 5);

    // double start_y = std::floor(miny - 5);
    // double end_y   = std::ceil(maxy + 5);

    // double start_z = std::floor(minz - 5);
    // double end_z   = std::ceil(maxz + 5);

    double start_x = std::floor(minx);
    double end_x   = std::ceil(maxx);

    double start_y = std::floor(miny);
    double end_y   = std::ceil(maxy);

    double start_z = std::floor(minz);
    double end_z   = std::ceil(maxz);

    int nx = (int)((end_x - start_x) / grid_spacing + 1e-5) + 1;
    int ny = (int)((end_y - start_y) / grid_spacing + 1e-5) + 1;
    int nz = (int)((end_z - start_z) / grid_spacing + 1e-5) + 1;


    long total_reps = (long)nx * ny * nz;

    std::cout << "-> Processing Vertices" << std::endl;


    std::vector<Vertex> mySurface = vert_to_vector(vert_file);

    std::vector<Vec3> insidePoints;
    std::vector<Vec3> outsidePoints;
    
    Vec3 minB = {(float)start_x - 5, (float)start_y - 5, (float)start_z - 5};
    Vec3 maxB = {(float)end_x + 5, (float)end_y + 5, (float)end_z + 5};
    
    SeparateGridPoints(mySurface, minB, maxB, 0.25f, 2.0f, insidePoints, outsidePoints);

    WriteWaterPDB(insidePoints, output_file + "_in.pdb");
    WriteWaterPDB(outsidePoints, output_file + "_out.pdb");

    std::cout << insidePoints.size() << "  " << outsidePoints.size() << std::endl;
    /*
    std::cout << "-> Filling Void" << std::endl;

    
    std::cout  << "-> total gridpoints = " << std::scientific << std::setprecision(3) << (double)total_reps << std::endl;
    std::cout << std::scientific << std::setprecision(3) << "-- removed " << (double)total_reps - (double)finalVolume.size() << " --" << std::endl;
    std::cout << std::scientific << std::setprecision(3) << "-> Dowsing " << (double)finalVolume.size() << std::endl;


    std::vector<Atom> watervector = {};
    int dbug = 0;

    // for(double x = std::floor(minx - 5); x <= std::ceil(maxx + 5); x += grid_spacing) {
    //     for (double y = std::floor(miny - 5); y <= std::ceil(maxy + 5); y+= grid_spacing){
    //         for( double z = std::floor(minz - 5); z <= std::ceil(maxz + 5); z+= grid_spacing) {

    //             std::cout << "\r** Iteration: " << dbug+1 << " of " << total_reps << " **" << std::flush;

    //             // std::cout << "x = " << x << "   y = " << y << "   z = " << z << std::endl;
    //             std::array<double, 3> temp_array = {x, y ,z};
    //             if(getOverlap_cluster(map, atomvector, temp_array, hash_spacing, water_diameter, cutoff_distance)){
    //                 Atom newatom("HOH", "O", temp_array);
    //                 watervector.push_back(newatom);
    //             }
    //             dbug += 1;
    //         }
    //     }
    // }
    for (int i = 0; i < finalVolume.size(); i++){
        std::cout << "\r** Iteration: " << dbug+1 << " of " << finalVolume.size() << " **" << std::flush;

        std::array<double, 3> temp_array = {finalVolume[i].position.x, finalVolume[i].position.y ,finalVolume[i].position.z};
        if(getOverlap_cluster(map, atomvector, temp_array, hash_spacing, water_diameter, cutoff_distance)){
            Atom newatom("HOH", "O", temp_array);
            watervector.push_back(newatom);
        }
        dbug += 1;

    }

    std::cout << "\n** There are " << watervector.size() << " waters **" <<std::endl;


    std::cout << "-> Entering vectortopdb" << std::endl;
    
    vectortopdb(watervector, output_file);

    std::cout << "-> Success" << std::endl;
    
    */
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Time: " << elapsed.count() << " s" << std::endl;

    return 1;
}