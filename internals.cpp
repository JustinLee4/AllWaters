#include "internals.h" 

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
    double grid_y = (vector.y - minBound.x) * (1/spacing);
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

void GenerateShell_Pull(
    const std::vector<Vertex>& surfaceVertices,
    Vec3 minBound,
    Vec3 maxBound,
    float spacing,      // Fine grid spacing (e.g., 0.25)
    float searchRadius, // Max distance to search (e.g., 1.5)
    std::vector<Vec3>& outInside, 
    std::vector<Vec3>& outOutside
) {
    // ---------------------------------------------------------
    // PHASE 1: Build Acceleration Grid (Spatial Hashing)
    // ---------------------------------------------------------
    
    // The bin size should effectively match the search radius.
    // This guarantees we only check 3x3x3 bins.
    float binSize = searchRadius; 

    int binDimX = static_cast<int>(std::ceil((maxBound.x - minBound.x) / binSize)) + 1;
    int binDimY = static_cast<int>(std::ceil((maxBound.y - minBound.y) / binSize)) + 1;
    int binDimZ = static_cast<int>(std::ceil((maxBound.z - minBound.z) / binSize)) + 1;

    // A flat vector of vectors. 
    // storage[idx] contains a list of vertex indices that live in that bin.
    std::vector<std::vector<int>> bins(binDimX * binDimY * binDimZ);

    for (int i = 0; i < surfaceVertices.size(); i++) {
        const Vec3& p = surfaceVertices[i].position;
        
        // Map vertex to Bin Coordinate
        int bx = static_cast<int>((p.x - minBound.x) / binSize);
        int by = static_cast<int>((p.y - minBound.y) / binSize);
        int bz = static_cast<int>((p.z - minBound.z) / binSize);

        // Bounds safety check
        if (bx >= 0 && bx < binDimX && by >= 0 && by < binDimY && bz >= 0 && bz < binDimZ) {
            size_t binIdx = (size_t)bz * (binDimX * binDimY) + (size_t)by * binDimX + bx;
            bins[binIdx].push_back(i);
        }
    }

    // ---------------------------------------------------------
    // PHASE 2: The "Pull" Loop (Iterate Voxels)
    // ---------------------------------------------------------

    // Dimensions of the FINE grid (your output resolution)
    int gridDimX = static_cast<int>(std::ceil((maxBound.x - minBound.x) / spacing)) + 1;
    int gridDimY = static_cast<int>(std::ceil((maxBound.y - minBound.y) / spacing)) + 1;
    int gridDimZ = static_cast<int>(std::ceil((maxBound.z - minBound.z) / spacing)) + 1;
    float searchRadiusSq = searchRadius * searchRadius;

    // Iterate over every single target voxel
    for (int z = 0; z < gridDimZ; z++) {
        for (int y = 0; y < gridDimY; y++) {
            for (int x = 0; x < gridDimX; x++) {
                
                // 1. Calculate World Position of this voxel
                Vec3 voxelPos;
                voxelPos.x = minBound.x + (x * spacing);
                voxelPos.y = minBound.y + (y * spacing);
                voxelPos.z = minBound.z + (z * spacing);

                // 2. Which Bin is this voxel inside?
                int bx = static_cast<int>((voxelPos.x - minBound.x) / binSize);
                int by = static_cast<int>((voxelPos.y - minBound.y) / binSize);
                int bz = static_cast<int>((voxelPos.z - minBound.z) / binSize);

                int closestVertIdx = -1;
                float minD = std::numeric_limits<float>::max();

                // 3. Search Neighboring Bins (3x3x3 Neighborhood)
                // Since BinSize == SearchRadius, the closest point MUST be in this bin 
                // or one of the 26 immediate neighbors.
                for (int nz = bz - 1; nz <= bz + 1; nz++) {
                    for (int ny = by - 1; ny <= by + 1; ny++) {
                        for (int nx = bx - 1; nx <= bx + 1; nx++) {
                            
                            if (nx < 0 || nx >= binDimX || ny < 0 || ny >= binDimY || nz < 0 || nz >= binDimZ)
                                continue;

                            size_t neighborBinIdx = (size_t)nz * (binDimX * binDimY) + (size_t)ny * binDimX + nx;

                            // Check all vertices in this bin
                            for (int vIdx : bins[neighborBinIdx]) {
                                const Vertex& v = surfaceVertices[vIdx];
                                float d2 = (voxelPos - v.position).lengthSq();

                                if (d2 < minD) {
                                    minD = d2;
                                    closestVertIdx = vIdx;
                                }
                            }
                        }
                    }
                }

                // 4. Result Processing
                // If we found a valid vertex within range, classify the point
                if (closestVertIdx != -1 && minD <= searchRadiusSq) {
                    const Vertex& closest = surfaceVertices[closestVertIdx];
                    Vec3 dir = voxelPos - closest.position;
                    
                    // Simple dot product check
                    // Note: For perfect shells, you might need a more robust sign check,
                    // but this matches your current logic.
                    if (closest.normal.x * dir.x + closest.normal.y * dir.y + closest.normal.z * dir.z < 0) {
                        outInside.push_back(voxelPos);
                    } else {
                        outOutside.push_back(voxelPos);
                    }
                }
            }
        }
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
    
    printf("Successfully wrote %zu waters to %s\n", waterPositions.size(), filename.c_str());
}