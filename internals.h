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
#include <cstdio> // For FILE, fprintf




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