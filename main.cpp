#include "AtomicRadii.h"
#include "pdbtovector.h"
#include "map.h"

#include <iostream>


double hash_spacing = 3;
double grid_spacing = .5;
double water_diameter = 2.5;
double cutoff_distance = 3.5;


int main(int argc, char* argv[]){

    // AtomParams params = getParams("ALA", "N");
    // std::cout << params.radius_ua << std::endl;

    std::string input_file = "8OM1_structure_met.pdb";

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

    double start_x = std::floor(minx - 5);
    double end_x   = std::ceil(maxx + 5);

    double start_y = std::floor(miny - 5);
    double end_y   = std::ceil(maxy + 5);

    double start_z = std::floor(minz - 5);
    double end_z   = std::ceil(maxz + 5);

    int nx = (int)((end_x - start_x) / grid_spacing + 1e-5) + 1;
    int ny = (int)((end_y - start_y) / grid_spacing + 1e-5) + 1;
    int nz = (int)((end_z - start_z) / grid_spacing + 1e-5) + 1;

    long total_reps = (long)nx * ny * nz;

    std::cout << "-> Dowsing " << total_reps << std::endl;


    std::vector<Atom> watervector = {};
    int dbug = 0;

    for(double x = std::floor(minx - 5); x <= std::ceil(maxx + 5); x += grid_spacing) {
        for (double y = std::floor(miny - 5); y <= std::ceil(maxy + 5); y+= grid_spacing){
            for( double z = std::floor(minz - 5); z <= std::ceil(maxz + 5); z+= grid_spacing) {

                std::cout << "\r** Iteration: " << dbug+1 << " of " << total_reps << " **" << std::flush;

                // std::cout << "x = " << x << "   y = " << y << "   z = " << z << std::endl;
                std::array<double, 3> temp_array = {x, y ,z};
                if(getOverlap_cluster(map, atomvector, temp_array, hash_spacing, water_diameter, cutoff_distance)){
                    Atom newatom("HOH", "O", temp_array);
                    watervector.push_back(newatom);
                }
                dbug += 1;
            }
        }
    }

    std::cout << "\n** There are " << watervector.size() << " waters **" <<std::endl;


    std::cout << "-> Entering vectortopdb" << std::endl;
    
    vectortopdb(watervector, "waters-05-8OM1.pdb");

    std::cout << "-> Success" << std::endl;



    return 1;
}