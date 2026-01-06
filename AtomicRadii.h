#ifndef ATOM_RADII
#define ATOM_RADII

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

// 1. Define the structure for the values
struct AtomParams {
    double radius_ua;      // Column 0: Adjusted VDW (United Atom)
    double radius_aa;      // Column 1: All Atom Radius
    double hc_1986;        // Column 2: Hydrophobicity 1986
    double hc_1989;        // Column 3: Hydrophobicity 1989
    double hc_1998;        // Column 4: Hydrophobicity 1998
    int type_id;           // Column 5: Integer flag

};

using AtomKey = std::pair<std::string, std::string>;

extern const std::map<AtomKey, AtomParams> ATOMIC_RADII;

AtomParams getParams(std::string resName, std::string atomName);

#endif