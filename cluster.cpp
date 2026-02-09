#include "cluster.h"
#include <cmath>

std::vector<std::vector<Atom>> clusterAtoms(const std::vector<Atom>& atoms, double spacing) {
    //Build the spatial grid
    std::unordered_map<GridKey, std::vector<int>> grid = buildSpatialGrid(atoms, spacing);
    
    std::vector<std::vector<Atom>> clusters;
    std::vector<bool> visited(atoms.size(), false);

    // Define the maximum distance squared for adjacency (including diagonals)
    // We add a tiny epsilon (1.05 factor) to handle floating point inaccuracies.
    double maxDistSq = 3.0 * spacing * spacing * 1.1; 

    for (int i = 0; i < (int)atoms.size(); ++i) {
        if (visited[i]) continue;

        std::vector<Atom> currentCluster;
        std::queue<int> q;

        visited[i] = true;
        q.push(i);

        while (!q.empty()) {
            int currIdx = q.front();
            q.pop();
            currentCluster.push_back(atoms[currIdx]);

            std::array<double, 3> posA = atoms[currIdx].getCoords();
            GridKey k = getGridKey(atoms[currIdx], spacing);

            // Check 3x3x3 neighborhood
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        GridKey neighborKey = {k.x + dx, k.y + dy, k.z + dz};
                        
                        auto it = grid.find(neighborKey);
                        if (it != grid.end()) {
                            for (int neighborIdx : it->second) {
                                
                                if (visited[neighborIdx]) continue;

                                std::array<double, 3> posB = atoms[neighborIdx].getCoords();

                                double dX = posA[0] - posB[0];
                                double dY = posA[1] - posB[1];
                                double dZ = posA[2] - posB[2];
                                double distSq = dX*dX + dY*dY + dZ*dZ;

                                if (distSq <= maxDistSq) {
                                    visited[neighborIdx] = true;
                                    q.push(neighborIdx);
                                }
                            }
                        }
                    }
                }
            }
        }
        clusters.push_back(currentCluster);
    }
    return clusters;
}


void writeClusteredPDB(const std::vector<std::vector<Atom>>& clusters, const std::string& filename) {
    FILE* pFile = fopen((filename + ".pdb").c_str(), "w");
    if (pFile == NULL) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return;
    }

    int atomGlobalCount = 1;
    for (int i = 0; i < (int)clusters.size(); ++i) {
        // PDB residue IDs wrap at 9999
        int resID = (i % 9999) + 1;
        
        for (const auto& atom : clusters[i]) {
            std::array<double, 3> c = atom.getCoords();
            
            // Format: HETATM, ID, AtomName, ResName, Chain, ResID, X, Y, Z, Occ, B, Element
            fprintf(pFile, "HETATM%5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s\n",
                    atomGlobalCount % 100000, 
                    atom.get_atomname().c_str(),
                    atom.get_resname().c_str(),
                    resID,
                    c[0], c[1], c[2],
                    atom.get_atomname().substr(0,1).c_str());
            
            atomGlobalCount++;
        }
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}