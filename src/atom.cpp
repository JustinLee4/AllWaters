#include "atom.h"

#include "common.h"

std::string Atom::toString() {
    std::stringstream ss;

            //Write the formatted data into the stream

            ss << "  Resname : " << resname;
            ss << "  Atomname : " << atomname;
            ss << "  Position (x, y, z): [";
            
            // Print the array elements
            ss << position[0] << ", " 
            << position[1] << ", " 
            << position[2] << "]" << "\n";
            
            return ss.str();
}

std::array<double, 3> Atom::getCoords() const {
    return position;
}

void Atom::set_radius(double new_radius) {
    radius = new_radius;
}

double Atom::get_radius() {
    return radius;
}


std::string Atom::get_resname()const{
    return resname;
}

std::string Atom::get_atomname() const {
    return atomname;
}

double Atom::get_bfactor() const {
    return b_factor;
}

