#include "internals.h" 

#include "common.h"


std::vector<Vertex> vert_to_vector(std::string vert_file) {
    std::vector<Vertex> output;
    std::ifstream infile(vert_file);
    std::string line_string;
    
    int z = 0;
    while(std::getline(infile, line_string)) {
        if(!line_string.empty() & (z > 2)){

            std::stringstream ss(line_string);
            std::istream_iterator<std::string> start(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> tokens(start, end);

            Vertex v = {{std::stof(tokens[0]),
                std::stof(tokens[1]), 
                std::stof(tokens[2])},
                {std::stof(tokens[3]),
                std::stof(tokens[4]), 
                std::stof(tokens[5])}};
            output.push_back(v);


        }
        z++;
    }
    return output;
}

int vectoIndex(Vec3 vector, int dimX, int dimY, int dimZ, Vec3 minBound, double spacing) {

    double grid_x = (vector.x - minBound.x) * (1/spacing);
    double grid_y = (vector.y - minBound.y) * (1/spacing);
    double grid_z = (vector.z - minBound.z) * (1/spacing);

    int idx = grid_z * (dimX * dimY) + grid_y * dimX + grid_x;
    return idx;
}
Vec3 indextoVec3(int index, int dimX, int dimY, int dimZ, Vec3 minBound, double spacing) {
    double z = index / (dimX * dimY);
    int remainder = index % (dimX * dimY);
    double y = remainder / dimX;
    double x = remainder % dimX;

    Vec3 vector;
    vector.x = (minBound.x) + (x * spacing);
    vector.y = (minBound.y) + (y * spacing);
    vector.z = (minBound.z) + (z * spacing);

    return vector;
}

void SeparateGridPoints(
    const std::vector<Vertex>& surfaceVertices,
    Vec3 minBound,
    Vec3 maxBound,
    float spacing,
    float searchRadius,
    std::vector<Vec3>& outInside, 
    std::vector<Vec3>& outOutside
) {
    
    int dimX = static_cast<int>(std::ceil((maxBound.x - minBound.x) / spacing));
    int dimY = static_cast<int>(std::ceil((maxBound.y - minBound.y) / spacing));
    int dimZ = static_cast<int>(std::ceil((maxBound.z - minBound.z) / spacing));

    size_t totalCells = (size_t)dimX * dimY * dimZ;
    
    std::vector<GridCellInfo> grid(totalCells);

    float searchRadiusSq = searchRadius * searchRadius;
    int searchRadius_cells = static_cast<int>(std::ceil(searchRadius / spacing));

    for(int i = 0; i < surfaceVertices.size(); i++){
        const Vertex& vert = surfaceVertices[i];
        
        int cx = static_cast<int>(std::floor((vert.position.x - minBound.x) / spacing));
        int cy = static_cast<int>(std::floor((vert.position.y - minBound.y) / spacing));
        int cz = static_cast<int>(std::floor((vert.position.z - minBound.z) / spacing));

        for (int z = cz - searchRadius_cells; z <= cz + searchRadius_cells; z++) {
            for (int y = cy - searchRadius_cells; y <= cy + searchRadius_cells; y++) {
                for (int x = cx - searchRadius_cells; x <= cx + searchRadius_cells; x++) {   
                    
                    if(x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ) {
                        continue;
                    }

                    Vec3 gridPos;
                    gridPos.x = minBound.x + (x * spacing);
                    gridPos.y = minBound.y + (y * spacing);
                    gridPos.z = minBound.z + (z * spacing);

                    double distanceSq = (gridPos - vert.position).lengthSq();

                    if(distanceSq <= searchRadiusSq) {
                        
                        size_t idx = (size_t)z * (dimX * dimY) + (size_t)y * dimX + (size_t)x;

                        if (distanceSq < grid[idx].minDistSq) {
                            grid[idx].closestVertexIndex = i;
                            grid[idx].minDistSq = distanceSq;
                        }
                    }
                }
            }
        }
    }

    // --- PHASE 2: CLASSIFY ---
    for(int z = 0; z < dimZ; z++) {
        for(int y = 0; y < dimY; y++) {
            for(int x = 0; x < dimX; x++) {
                
                size_t i = (size_t)z * (dimX * dimY) + (size_t)y * dimX + (size_t)x;

                if(grid[i].closestVertexIndex != -1){
                    
                    Vec3 test_point;
                    test_point.x = minBound.x + (x * spacing);
                    test_point.y = minBound.y + (y * spacing);
                    test_point.z = minBound.z + (z * spacing);

                    const Vertex& closestVert = surfaceVertices[grid[i].closestVertexIndex];
                    Vec3 dir = test_point - closestVert.position;
                    
                    // Dot Product check
                    float dotproduct = closestVert.normal.dot(dir);

                    if (dotproduct < 0) {
                        outInside.push_back(test_point);
                    } else {
                        outOutside.push_back(test_point);
                    }
                }
            }
        }
    }
}

//weighted voting
// void SeparateGridPoints(
//     const std::vector<Vertex>& surfaceVertices,
//     Vec3 minBound,
//     Vec3 maxBound,
//     float spacing,
//     float searchRadius,
//     std::vector<Vec3>& outInside, 
//     std::vector<Vec3>& outOutside
// ) {
//     // --- SETUP GRID DIMENSIONS ---
//     int dimX = static_cast<int>(std::ceil((maxBound.x - minBound.x) / spacing));
//     int dimY = static_cast<int>(std::ceil((maxBound.y - minBound.y) / spacing));
//     int dimZ = static_cast<int>(std::ceil((maxBound.z - minBound.z) / spacing));
    
//     size_t totalCells = (size_t)dimX * dimY * dimZ;

//     // 1. ALLOCATE: One large contiguous block
//     // 50M cells * ~30 bytes = ~1.5 GB RAM. Safe for modern machines.
//     std::vector<CellVotes> grid(totalCells);

//     float searchRadiusSq = searchRadius * searchRadius;
//     int searchRadius_cells = static_cast<int>(std::ceil(searchRadius / spacing));

//     // --- PHASE 1: COLLECT CANDIDATES (SCATTER) ---
//     // Iterate over surface vertices and "scatter" them into nearby grid cells
//     for(int i = 0; i < surfaceVertices.size(); i++){
//         const Vertex& vert = surfaceVertices[i];
        
//         // Find which cell this vertex falls into
//         int cx = static_cast<int>(std::floor((vert.position.x - minBound.x) / spacing));
//         int cy = static_cast<int>(std::floor((vert.position.y - minBound.y) / spacing));
//         int cz = static_cast<int>(std::floor((vert.position.z - minBound.z) / spacing));

//         // Search neighboring cells within radius
//         for (int z = cz - searchRadius_cells; z <= cz + searchRadius_cells; z++) {
//             for (int y = cy - searchRadius_cells; y <= cy + searchRadius_cells; y++) {
//                 for (int x = cx - searchRadius_cells; x <= cx + searchRadius_cells; x++) {   
                    
//                     // Bounds Check
//                     if(x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ) continue;

//                     // Calculate Grid Point Position
//                     Vec3 gridPos;
//                     gridPos.x = minBound.x + (x * spacing);
//                     gridPos.y = minBound.y + (y * spacing);
//                     gridPos.z = minBound.z + (z * spacing);

//                     // Quick Box Check (Avoid sqrt)
//                     if (std::abs(gridPos.x - vert.position.x) > searchRadius) continue;
//                     if (std::abs(gridPos.y - vert.position.y) > searchRadius) continue;
//                     if (std::abs(gridPos.z - vert.position.z) > searchRadius) continue;

//                     // Precise Distance Check
//                     double d2 = (gridPos - vert.position).lengthSq();

//                     if(d2 <= searchRadiusSq) {
//                         size_t idx = (size_t)z * (dimX * dimY) + (size_t)y * dimX + (size_t)x;
                        
//                         // INSERT into fixed-size priority queue
//                         grid[idx].insert(i, (float)d2);
//                     }
//                 }
//             }
//         }
//     }

    // --- PHASE 2: WEIGHTED VOTING ---
    // Iterate over grid cells and classify them based on their collected candidates
//     for(int z = 0; z < dimZ; z++) {
//         for(int y = 0; y < dimY; y++) {
//             for(int x = 0; x < dimX; x++) {
                
//                 size_t idx = (size_t)z * (dimX * dimY) + (size_t)y * dimX + (size_t)x;
//                 CellVotes& cell = grid[idx];

//                 // Only process cells that found at least one surface vertex
//                 if(cell.count > 0){
//                     Vec3 test_point;
//                     test_point.x = minBound.x + (x * spacing);
//                     test_point.y = minBound.y + (y * spacing);
//                     test_point.z = minBound.z + (z * spacing);

//                     double scoreInside = 0.0;
//                     double scoreOutside = 0.0;
                    
//                     // Iterate through top 3 (or fewer) candidates
//                     for(int j = 0; j < cell.count; j++) {
//                         const Vertex& v = surfaceVertices[cell.indices[j]];
//                         Vec3 dir = test_point - v.position;
                        
//                         // WEIGHT CALCULATION: Inverse Distance Weighting
                        
//                         // 1. Calculate Weight (1 / distance)
//                         // Add epsilon to avoid divide-by-zero
//                         float weight = 1.0f / (cell.distSq[j]  + 1e-6f);

//                         // 2. Dot Product Check
//                         if (v.normal.dot(dir) < 0) {
//                             scoreInside += weight;
//                         } else {
//                             scoreOutside += weight;
//                         }
//                     }

//                     // DECISION: Whichever score is higher wins
//                     if (scoreInside > scoreOutside) {
//                         outInside.push_back(test_point);
//                     } else {
//                         outOutside.push_back(test_point);
//                     }
//                 }
//             }
//         }
//     }
// }

void FillInternalVoid(
    const std::vector<Vec3>& shellPoints,
    const std::vector<Vec3>& initialInsidePoints,
    Vec3 minBound,
    Vec3 maxBound,
    float spacing,
    std::vector<Vec3>& allPoints,
    std::vector<std::vector<Vec3>>& layers
) {
    // Clear outputs to ensure fresh start
    allPoints.clear();
    layers.clear();

    // 1. Setup Grid Dimensions
    // Added +1 buffer to ensure coverage
    int dimX = static_cast<int>(std::ceil((maxBound.x - minBound.x) / spacing)) + 1;
    int dimY = static_cast<int>(std::ceil((maxBound.y - minBound.y) / spacing)) + 1;
    int dimZ = static_cast<int>(std::ceil((maxBound.z - minBound.z) / spacing)) + 1;
    size_t totalCells = (size_t)dimX * dimY * dimZ;

    // 2. The "Marked" Map (visited/wall tracker)
    std::vector<bool> marked(totalCells, false);

    // 3. Mark the Shells (The Walls)
    for (const auto& p : shellPoints) {
        int idx = vectoIndex(p, dimX, dimY, dimZ, minBound, spacing);
        if (idx != -1) {
            marked[idx] = true; 
        }
    }

    // 4. Initialize Queue with Seeds
    std::vector<int> currentLayerIndices;
    std::vector<Vec3> layer0;
    
    for (const auto& p : initialInsidePoints) {
        int idx = vectoIndex(p, dimX, dimY, dimZ, minBound, spacing);
        
        // Only process if valid and NOT already marked (shell or duplicate)
        if (idx != -1 && !marked[idx]) {
            marked[idx] = true; 
            currentLayerIndices.push_back(idx);
            
            // Reconstruct point to ensure grid alignment consistency
            Vec3 alignedP = indextoVec3(idx, dimX, dimY, dimZ, minBound, spacing);
            
            allPoints.push_back(alignedP);
            layer0.push_back(alignedP);
        }
    }
    
    if (!layer0.empty()) {
        layers.push_back(layer0);
    }

    // 5. BFS Loop (Flood Fill)
    // Neighbors: +X, -X, +Y, -Y, +Z, -Z
    const int dx[6] = {1, -1, 0, 0, 0, 0};
    const int dy[6] = {0, 0, 1, -1, 0, 0};
    const int dz[6] = {0, 0, 0, 0, 1, -1};

    while (!currentLayerIndices.empty()) {
        std::vector<int> nextLayerIndices;
        std::vector<Vec3> nextLayerVecs;

        for (int currentIdx : currentLayerIndices) {
            
            // Manual index unpacking for neighbor calculation
            int slice = dimX * dimY;
            int cz = currentIdx / slice;
            int rem = currentIdx % slice;
            int cy = rem / dimX;
            int cx = rem % dimX;

            for (int i = 0; i < 6; i++) {
                int nx = cx + dx[i];
                int ny = cy + dy[i];
                int nz = cz + dz[i];

                // Bounds Check
                if (nx >= 0 && nx < dimX && ny >= 0 && ny < dimY && nz >= 0 && nz < dimZ) {
                    
                    int nIdx = nz * slice + ny * dimX + nx;

                    if (!marked[nIdx]) {
                        marked[nIdx] = true; // Mark visited immediately
                        nextLayerIndices.push_back(nIdx);
                        
                        Vec3 nPos = indextoVec3(nIdx, dimX, dimY, dimZ, minBound, spacing);
                        
                        nextLayerVecs.push_back(nPos);
                        allPoints.push_back(nPos);
                    }
                }
            }
        }

        if (!nextLayerVecs.empty()) {
            layers.push_back(nextLayerVecs);
        }

        currentLayerIndices = nextLayerIndices;
    }
}


void WriteWaterPDB(const std::vector<Vec3>& waterPositions, const std::string& filename) {
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename.c_str());
        return;
    }

    // PDB Header
    fprintf(file, "REMARK    GENERATED BY SOLVATION TOOL\n");
    fprintf(file, "REMARK    THIS FILE CONTAINS %zu WATER MOLECULES\n", waterPositions.size());

    int serial = 1;  
    int resSeq = 1;  

    for (const auto& pos : waterPositions) {
        // Safety Check: Skip NaNs or Infinite numbers (these corrupt PDBs)
        if (std::isnan(pos.x) || std::isnan(pos.y) || std::isnan(pos.z)) continue;

        // Safety Check: Wrap indices BEFORE they print to prevent column shifting.
        // PDB columns are strict: Serial is 5 chars max, ResSeq is 4 chars max.
        int printSerial = serial > 99999 ? serial % 100000 : serial;
        int printResSeq = resSeq > 9999 ? resSeq % 10000 : resSeq;
        
        // Handle the "0" case for wrapped indices (optional, keeps IDs clean)
        if (printSerial == 0) printSerial = 1;
        if (printResSeq == 0) printResSeq = 1;

        // PDB ATOM Record Format (Strict Column alignment)
        // 1-4:   "ATOM"
        // 7-11:  Serial (5d)
        // 13-16: Name " OW " (Unique alignment for water)
        // 18-20: Residue "HOH"
        // 22:    Chain 'A'
        // 23-26: ResSeq (4d)
        // 31-54: Coords (8.3f)
        fprintf(file, "ATOM  %5d  OW  HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           O  \n",
            printSerial,
            printResSeq,
            pos.x,
            pos.y,
            pos.z
        );

        serial++;
        resSeq++;
    }

    fprintf(file, "END\n");
    fclose(file);
    
    printf("Successfully wrote %.3e waters to %s\n", (double)waterPositions.size(), filename.c_str());
}