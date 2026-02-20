#ifndef PYMOL_H
#define PYMOL_H

#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
#include <thread>
#include <chrono>


void createPyMOLSession(const std::string& pdb_file, const std::string& session_name, const std::string& structure_file = "");

#endif