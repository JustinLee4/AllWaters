#include "map.h"

#include "AtomicRadii.h"
#include "common.h"


GridKey getGridKey(const Atom& input_water, double gridCellSize) {
    std::array<double, 3> pos = input_water.getCoords();
    return {
        (int)std::floor(pos[0] / gridCellSize),
        (int)std::floor(pos[1] / gridCellSize),
        (int)std::floor(pos[2] / gridCellSize)
    };
}

GridKey getGridKey_pos(const std::array<double, 3> pos, double gridCellSize) {
    return {
        (int)std::floor(pos[0] / gridCellSize),
        (int)std::floor(pos[1] / gridCellSize),
        (int)std::floor(pos[2] / gridCellSize)
    };
}

std::unordered_map<GridKey, std::vector<int>>
buildSpatialGrid(std::vector<Atom> objects, double gridCellSize) {

    std::unordered_map<GridKey, std::vector<int>> grid;

    for (int i = 0; i < objects.size(); ++i) {
        GridKey key = getGridKey(objects[i], gridCellSize);
        grid[key].push_back(i);
    }
    return grid;
}

void printSpatialGrid(const std::unordered_map<GridKey, std::vector<int>>& grid) {
    
    std::cout << "--- Spatial Grid Contents ---" << std::endl;
    // grid.size() tells you how many unique cells (keys) were created.
    std::cout << "Total unique cells (keys): " << grid.size() << std::endl;
    std::cout << "-------------------------------" << std::endl;

    // We loop through the map. 'pair' will be one {key, value} item.
    for (const auto& pair : grid) {
        // 'pair.first' is the GridKey
        const GridKey& key = pair.first;
        
        // 'pair.second' is the std::vector<int>
        const std::vector<int>& objectIndices = pair.second;
        
        // Print the Key
        std::cout << "Key (Cell): (" << key.x << ", " << key.y << ", " << key.z << ")" << std::endl;
        
        // Print the Value (the list of indices)
        std::cout << "  Indices : { ";
        for (int index : objectIndices) {
            std::cout << index << " ";
        }
        std::cout << "}" << std::endl;
    }
    std::cout << "-------------------------------" << std::endl;
}


//find overlaps and also update nearest neighbors
bool getOverlap_cluster(const std::unordered_map<GridKey, std::vector<int>>& grid, std::vector<Atom>& Atomvector, std::array<double,3> target, double gridCellSize, double diameter, double cutoff_dist) {

    double cutoff_dist_sq = (cutoff_dist * cutoff_dist);
    double targetRadius = diameter / 2.0; 


    // 1. Get the center coordinates of the target particle
    GridKey centerKey = getGridKey_pos(target, gridCellSize);

    bool at_least_one_neighbor = false;

    // 2. Loop through X, Y, Z offsets (-1, 0, 1)
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                
                // 3. Construct the neighbor key
                // Note: Assuming GridKey is a struct with members x, y, z
                GridKey neighborKey;
                neighborKey.x = centerKey.x + dx;
                neighborKey.y = centerKey.y + dy;
                neighborKey.z = centerKey.z + dz;

                // 4. Try to find this neighbor box in the map
                auto it = grid.find(neighborKey);

                // 5. If the box exists (contains particles)
                if (it != grid.end()) {

                    const std::vector<int>& indices = it->second;
                    
                    for (int neighborIndex : indices) {  

                        std::array<double,3> origCoords = target;
                        std::array<double,3> targetCoords = Atomvector[neighborIndex].getCoords();

                        double distance = 0;
                        double dX = origCoords[0] - targetCoords[0];
                        double dY = origCoords[1] - targetCoords[1];
                        double dZ = origCoords[2] - targetCoords[2];

                        distance = (dX*dX) + (dY*dY) + (dZ*dZ);

                        AtomParams params = getParams(Atomvector[neighborIndex].get_resname(), Atomvector[neighborIndex].get_atomname());
                        double atomradius =  params.radius_aa;

                        double collisionThreshold = targetRadius + atomradius;

                        if(distance <= (collisionThreshold * collisionThreshold)) {
                            return false;
                        }
                        if (distance <= cutoff_dist_sq) {
                            at_least_one_neighbor = true;
                        }
                    }
                }
            } 
        }
    }
    return at_least_one_neighbor;
}



// // void testGrid() {
//     // 1. Define a simple test case
//     double myGridCellSize = 10.0;
    
//     std::vector<Atom> testObjects;
    
//     // Manually create some Atom objects.
//     // Assume Atom has a constructor or method to set its coords.
//     Atom obj0 = Atom(0, 0,std::array<double, 3> {5.0, 0.0, 0.0}); // obj0 at (5.0, 0.0, 0.0)
    
//     Atom obj1 = Atom(0, 0, std::array<double, 3> {12.0, 0.0, 0.0}); // obj1 at (12.0, 0.0, 0.0)

//     Atom obj2 = Atom(0, 0, std::array<double, 3> {15.0, 0.0, 0.0}); // obj2 at (15.0, 0.0, 0.0)

//     Atom obj3 = Atom(0, 0, std::array<double, 3> {-8.0, 0.0, 0.0}); // obj3 at (-8.0, 0.0, 0.0)

//     testObjects.push_back(obj0); // Index 0
//     testObjects.push_back(obj1); // Index 1
//     testObjects.push_back(obj2); // Index 2
//     testObjects.push_back(obj3); // Index 3

//     // 2. Manually calculate your expected result:
//     // With grid size 10.0:
//     // obj0 (5.0) -> key (0, 0, 0)
//     // obj1 (12.0) -> key (1, 0, 0)
//     // obj2 (15.0) -> key (1, 0, 0)
//     // obj3 (-8.0) -> key (-1, 0, 0)

//     // 3. Build the grid
//     std::unordered_map<GridKey, std::vector<int>> myGrid = 
//         buildSpatialGrid(testObjects, myGridCellSize);

//     // 4. Print and verify
//     printSpatialGrid(myGrid);
// }

