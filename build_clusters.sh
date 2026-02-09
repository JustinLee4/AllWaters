#!/bin/bash

g++ -O3 -std=c++17 find_clusters.cpp Atom_Lookup.cpp AtomicRadii_Map.cpp cluster.cpp pdbtovector.cpp atom.cpp map.cpp -o find_clusters
