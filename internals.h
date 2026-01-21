#ifndef GRID_PROCESSOR_H
#define GRID_PROCESSOR_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <cstdio>
#include <cfloat>




#include <chrono>


// --- Data Structures ---

struct Vec3 {
    float x, y, z;

    // Inline operators allow this struct to be defined in the header
    // without causing linker errors.
    inline Vec3 operator+(const Vec3& b) const { return {x + b.x, y + b.y, z + b.z}; }
    inline Vec3 operator-(const Vec3& b) const { return {x - b.x, y - b.y, z - b.z}; }
    inline Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    
    inline float dot(const Vec3& b) const { return x * b.x + y * b.y + z * b.z; }
    
    // Length Squared is faster than length() because it avoids sqrt()
    inline float lengthSq() const { return x * x + y * y + z * z; }
};

struct Vertex {
    Vec3 position;
    Vec3 normal;
};

struct GridCellInfo {
    float minDistSq = std::numeric_limits<float>::max();
    int closestVertexIndex = -1;
};

// struct CellVotes {
//     int indices[3];
//     float dists[3];
//     int count; // How many valid candidates we have found so far (0, 1, 2, or 3)

//     CellVotes() {
//         count = 0;
//         for(int i=0; i<3; i++) {
//             indices[i] = -1;
//             dists[i] = FLT_MAX;
//         }
//     }

//     // Smart Insert: Keeps the array sorted (Smallest distance at [0])
//     // This eliminates the need for std::sort later
//     void insert(int index, float dist) {
//         // If dist is worse than our worst candidate, ignore it
//         if (count == 3 && dist >= dists[2]) return;

//         // Find insertion point
//         int i;
//         for (i = count - 1; i >= 0; i--) {
//             if (dists[i] > dist) {
//                 if (i + 1 < 3) { // Shift down
//                     dists[i + 1] = dists[i];
//                     indices[i + 1] = indices[i];
//                 }
//             } else {
//                 break;
//             }
//         }
        
//         // Insert at the correct spot
//         if (i + 1 < 3) {
//             dists[i + 1] = dist;
//             indices[i + 1] = index;
//             if (count < 3) count++;
//         }
//     }
// };

struct CellVotes {
    int indices[3];    // Index into surfaceVertices
    float distSq[3];   // Squared distance to that vertex
    int count;         // How many valid candidates found (0-3)

    CellVotes() {
        count = 0;
        for(int i=0; i<3; i++) {
            indices[i] = -1;
            distSq[i] = FLT_MAX;
        }
    }

    // Smart Insert: Keeps the array sorted (Smallest distance at [0])
    void insert(int index, float d2) {
        // Optimization: If full and new point is worse than worst, skip immediately
        if (count == 3 && d2 >= distSq[2]) return;

        // Find insertion point (Linear scan is fastest for N=3)
        int i;
        for (i = count - 1; i >= 0; i--) {
            if (distSq[i] > d2) {
                if (i + 1 < 3) { // Shift down
                    distSq[i + 1] = distSq[i];
                    indices[i + 1] = indices[i];
                }
            } else {
                break;
            }
        }
        
        // Insert at the correct spot
        if (i + 1 < 3) {
            distSq[i + 1] = d2;
            indices[i + 1] = index;
            if (count < 3) count++;
        }
    }
};

// --- Function Declarations ---

std::vector<Vertex> vert_to_vector(std::string vert_file);

void SeparateGridPoints(
    const std::vector<Vertex>& surfaceVertices,
    Vec3 minBound,
    Vec3 maxBound,
    float spacing,
    float searchRadius,
    std::vector<Vec3>& outInside,  // Output Vector 1
    std::vector<Vec3>& outOutside  // Output Vector 2
);

void FillInternalVoid(
    const std::vector<Vec3>& shellPoints,
    const std::vector<Vec3>& initialInsidePoints,
    Vec3 minBound,
    Vec3 maxBound,
    float spacing,
    std::vector<Vec3>& allPoints,
    std::vector<std::vector<Vec3>>& layers
);

void WriteWaterPDB(const std::vector<Vec3>& waterPositions, const std::string& filename);



#endif // GRID_PROCESSOR_H