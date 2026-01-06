#ifndef WATER
#define WATER

#include <iostream>
#include <sstream>
#include <array>
#include <vector>

class Atom {
    std::string resname;
    std::string atomname;
    double radius;
    std::array<double, 3> position;
    
    public:
    Atom(std::string init_resname, std::string init_atomname, std::array<double,3> init_position):
        resname(init_resname),
        atomname(init_atomname),
        radius(0.0)
        {
            if (init_position.size() == 3) {
                std::copy(init_position.begin(), init_position.end(), position.begin());
            } else {
            // Handle error or default initialization if necessary
            std::cerr << "Warning: Position list must have exactly 3 elements." << std::endl;
            }
        }
        
    // print binary value and coordinates as a string
    std::string toString();

    //get coordinate information as an array of doubles
    std::array<double, 3> getCoords() const;

    void set_radius(double new_radius);

    std::string get_resname() const;

    std::string get_atomname() const;

    
    ~Atom() {
    }
};

#endif