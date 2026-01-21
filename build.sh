#!/bin/bash

#g++ -std=c++17 main.cpp Atom_Lookup.cpp AtomicRadii_Map.cpp pdbtovector.cpp atom.cpp map.cpp -o main
g++ -O3 -std=c++17 *.cpp -o main