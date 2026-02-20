#include "cluster.h"

#include "common.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

std::vector<std::vector<Atom>> clusterAtoms(const std::vector<Atom>& atoms, double grid_spacing, double map_spacing) {
    //Build the spatial grid
    std::unordered_map<GridKey, std::vector<int>> grid = buildSpatialGrid(atoms, map_spacing);
    
    std::vector<std::vector<Atom>> clusters;
    std::vector<bool> visited(atoms.size(), false);

    // Define the maximum distance squared for adjacency (including diagonals)
    // We add a tiny epsilon (1.05 factor) to handle floating point inaccuracies.
    double maxDistSq = 3.0 * grid_spacing * grid_spacing * 1.05; 

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
            GridKey k = getGridKey(atoms[currIdx], map_spacing);

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
    std::sort(clusters.begin(), clusters.end(), 
        [](const std::vector<Atom>& a, const std::vector<Atom>& b) {
            return a.size() > b.size();
        }    
    );

    return clusters;
}

void writeClusteredPDB(const std::vector<std::vector<Atom>>& clusters, 
                       const std::string& filename, 
                       const std::vector<std::string>& remarks) {
    // 1. Open File
    FILE* pFile = fopen((filename + ".pdb").c_str(), "w");
    
    if (pFile == NULL) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return;
    }

    // --- NEW: Write Comments (REMARK) at the top ---
    for (const auto& remark : remarks) {
        // PDB standard usually has at least one space after REMARK. 
        // Adding 4 spaces keeps it visually aligned with typical PDB headers.
        fprintf(pFile, "REMARK    %s\n", remark.c_str());
    }
    fprintf(pFile, "REMARK    --------------------------------\n");

    // 2. Write Atom Data
    int atomGlobalCount = 1;
    for (int i = 0; i < (int)clusters.size(); ++i) {
        // PDB residue IDs wrap at 9999
        int resID = (i % 9999) + 1;
        
        for (const auto& atom : clusters[i]) {
            std::array<double, 3> c = atom.getCoords();
            
            fprintf(pFile, "HETATM%5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %s\n",
                    atomGlobalCount % 100000, 
                    atom.get_atomname().c_str(),
                    atom.get_resname().c_str(),
                    resID,
                    c[0], c[1], c[2],
                    atom.get_bfactor(), 
                    atom.get_atomname().substr(0,1).c_str());
            
            atomGlobalCount++;
        }
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}

std::vector<std::string> extractRemarks(const std::string& filename) {
    std::vector<std::string> remarks;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open PDB file for reading remarks: " << filename << std::endl;
        return remarks; // Returns an empty vector safely
    }

    std::string line;
    while (std::getline(file, line)) {
        // 1. Check if the line starts with "REMARK"
        if (line.rfind("REMARK", 0) == 0) { // rfind at index 0 is a fast "starts_with" check
            
            // 2. Extract everything after the word "REMARK" (which is 6 characters long)
            std::string content = line.substr(6);
            
            // 3. Trim any leading whitespace so it aligns perfectly in your output file
            size_t first_char = content.find_first_not_of(" \t");
            if (first_char != std::string::npos) {
                content = content.substr(first_char);
            } else {
                content = ""; // Handles completely blank "REMARK" lines safely
            }
            
            remarks.push_back(content);
        }
        // 4. Stop reading early once we hit the actual atom data to save time!
        else if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) {
            break; 
        }
    }

    file.close();
    return remarks;
}