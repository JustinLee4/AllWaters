#ifndef CLUSTER
#define CLUSTER

#include "atom.h"
#include "map.h"

#include <queue>
#include <fstream>
#include <iostream>




std::vector<std::vector<Atom>> clusterAtoms(const std::vector<Atom>& atoms, double grid_spacing, double map_spacing);

void writeClusteredPDB(const std::vector<std::vector<Atom>>& clusters, 
                       const std::string& filename, 
                       const std::vector<std::string>& remarks = {});

std::vector<std::string> extractRemarks(const std::string& filename);

#endif