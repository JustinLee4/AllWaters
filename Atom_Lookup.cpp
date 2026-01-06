#include "AtomicRadii.h"


AtomParams getParams(std::string resName, std::string atomName) {

    // TRY 1: Look for exact match (e.g., "MET", "CA")
    // We create a pair {res, atom} to query the map
    auto it = ATOMIC_RADII.find({resName, atomName});
    
    if (it != ATOMIC_RADII.end()) {
        return it->second; // Found it!
    }

    // TRY 2: Look for generic/fallback match (e.g., "", "CA")
    // In Python this was (None, "CA"), in C++ we use empty string ""
    auto fallbackIt = ATOMIC_RADII.find({"", atomName});
    
    if (fallbackIt != ATOMIC_RADII.end()) {
        return fallbackIt->second; // Found generic version
    }

    // FAIL: Atom not found
    std::cerr << "WARNING: No parameters found for residue [" 
              << resName << "] atom [" << atomName << "]" << std::endl;

    // Return a zeroed-out struct
    return {0.0, 0.0, 0.0, 0.0, 0.0, 0};
}