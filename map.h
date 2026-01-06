#ifndef MAP
#define MAP

#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <iostream>
#include <array>

#include "atom.h"


struct GridKey {
    int x, y, z;

    // We must tell C++ how to check if two keys are equal
    bool operator==(const GridKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

namespace std {
    template <>
    struct hash<GridKey> {
        size_t operator()(const GridKey& k) const {
            // A common way to combine hashes: use xor (^) and bit-shifting.
            // This creates a good "hash" value from the three integers.
            return (
                (std::hash<int>()(k.x) ^
                 (std::hash<int>()(k.y) << 1)) >> 1
            ) ^ (std::hash<int>()(k.z) << 1);
        }
    };
}

GridKey getGridKey(const Atom& input_Atom, double gridCellSize);

GridKey getGridKey_pos(const std::array<double, 3> pos, double gridCellSize);


std::unordered_map<GridKey, std::vector<int>>
buildSpatialGrid(std::vector<Atom> objects, double gridCellSize);

void printSpatialGrid(const std::unordered_map<GridKey, std::vector<int>>& grid);

bool getOverlap_cluster(const std::unordered_map<GridKey, std::vector<int>>& grid, std::vector<Atom>& atomvector, std::array<double,3> target, double gridCellSize, double diameter, double cutoff_dist);

void testGrid();



#endif