#!/bin/bash

g++ -O3 -std=c++17 main.cpp Atom_Lookup.cpp AtomicRadii_Map.cpp atom.cpp internals.cpp map.cpp pdbtovector.cpp -o main
# g++ -O3 -std=c++17 *.cpp -o main