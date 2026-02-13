#ifndef CLUSTER
#define CLUSTER

#include "atom.h"
#include "map.h"

#include <queue>
#include <fstream>
#include <iostream>




std::vector<std::vector<Atom>> clusterAtoms(const std::vector<Atom>& atoms, double spacing);

void writeClusteredPDB(const std::vector<std::vector<Atom>>& clusters, const std::string& filename);

#endif