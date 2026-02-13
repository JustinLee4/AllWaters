#ifndef PDB_TO_VECTOR
#define PDB_TO_VECTOR

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "atom.h"


std::array<double, 3> get_coords(const std::string& input);


std::tuple<std::string, std::string> get_data(const std::string& input);


//turn pdbfile to a vector of atoms
std::tuple<std::vector<Atom>, double, double, double, double, double, double> pdbtovector(std::string filename);

//---------------------------------------------------
//converting vector back into pdb file

void vectortopdb(const std::vector<Atom> &atomvector, std::string output_filename);

//-----------------------------
//take a.pdb, append b.pdb, result c.pdb
bool append_pdb_files(const std::string& filepath_1, const std::string& filepath_2, const std::string& output_path, int targetLine);

#endif