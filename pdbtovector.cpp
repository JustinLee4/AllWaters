#include "pdbtovector.h"
//this is bugged when not separated by whitespace
// //get coords for Atom class
// std::array<double, 3> get_coords(std::string input) {
//     std::array<double,3> output;
//     std::string s;
//     std::stringstream ss(input);
   
//     // 1. Create the "start" iterator, initialized with the stringstream
//     std::istream_iterator<std::string> start(ss);
    
//     // 2. Create the "end" (default) iterator
//     std::istream_iterator<std::string> end;
    
//     // 3. Create the vector. This is now completely unambiguous.
//     std::vector<std::string> tokens(start, end);
//     // for(int i = 0; i < tokens.size(); i++){
//     //     std::cout << tokens[i] << ", i= " << i << std::endl;
//     // }
//     for(int i = 0; i < 3; i++){
//         output[i] = std::stod(tokens[i + 6]);
//     }
//     return output;
// }

std::array<double, 3> get_coords(const std::string& input) {
    // Initialize with 0.0 or a sentinel value (e.g., infinity)
    std::array<double, 3> output = {0.0, 0.0, 0.0};

    // 2. Extract Substrings based on PDB Standard (0-based indexing)
    // X: Columns 31-38 -> Index 30, Length 8
    // Y: Columns 39-46 -> Index 38, Length 8
    // Z: Columns 47-54 -> Index 46, Length 8
    try {
        std::string x_str = input.substr(30, 8);
        std::string y_str = input.substr(38, 8);
        std::string z_str = input.substr(46, 8);

        // 3. Convert to Double
        // std::stod automatically handles leading/trailing whitespace in the substring
        output[0] = std::stod(x_str);
        output[1] = std::stod(y_str);
        output[2] = std::stod(z_str);
    } 
    catch (...) {
        // Handles cases where columns are empty, contain "*******",
        // or contain non-numeric data.
        return {0.0, 0.0, 0.0};
    }

    return output;
}

// std::tuple<std::string, std::string> get_data(std::string input) {
//     std::string resname;
//     std::string atomname;
//     std::tuple<std::string, std::string> output;

//     std::stringstream ss(input);

//     std::istream_iterator<std::string> start(ss);
//     std::istream_iterator<std::string> end;

//     std::vector<std::string> tokens(start, end);

//     std::get<0>(output) = tokens[3];
//     std::get<1>(output) = tokens[2];

//     return output;

// }

// Helper function to trim whitespace from the result
std::string trim(const std::string& str) {
    auto first = str.find_first_not_of(' ');
    if (std::string::npos == first) {
        return str;
    }
    auto last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

std::tuple<std::string, std::string> get_data(const std::string& input) {
    // Default return values (Residue Name, Atom Name)
    std::string resName = "UNK";
    std::string atomName = " X  ";
    
    // PDB Fixed Column Widths (0-indexed):
    // Atom Name:    Columns 13-16 -> Index 12, Length 4
    // Residue Name: Columns 18-20 -> Index 17, Length 3
    
    try {
        if (input.length() > 20) {
            // Extract raw substrings
            std::string atom_raw = input.substr(12, 4);
            std::string res_raw = input.substr(17, 3);
            
            // Clean up whitespace (turn " O  " into "O")
            atomName = trim(atom_raw);
            resName = trim(res_raw);
        }
    } 
    catch (...) {
        // Fallback for empty or malformed lines
        // Returns defaults set above
    }

    // Return format matching your original tuple order: 
    // <0> Residue Name, <1> Atom Name
    return std::make_tuple(resName, atomName);
}


std::tuple<std::vector<Atom>, double, double, double, double, double, double> pdbtovector(std::string filename) {
    std::vector<Atom> output;
    std::ifstream infile(filename);
    double minx = INFINITY, maxx = -INFINITY, miny = INFINITY, maxy = -INFINITY, minz= INFINITY, maxz = -INFINITY;
    // std::ofstream out_file;
    // out_file.open("output.pdb");
    // out_file << "entered " << filename << std::endl;
    std::string line_string;
    int i = 0;
    while (std::getline(infile, line_string))
    {   
        i = i+1;
        // std::cout << i << std::endl;
        //Choosing HETATM and ATOM rows
        if(line_string.substr(0,6) == "HETATM") {
            //printing HETATM rows            
            // std::cout << "ans = " << line_string << std::endl;
            if(!line_string.empty()) {
                // out_file << line_string << "\n";
                std::array<double, 3> this_array = get_coords(line_string);
                std::tuple<std::string, std::string> this_data = get_data(line_string);
                Atom this_Atom = Atom(std::get<0>(this_data), std::get<1>(this_data), this_array);
                output.push_back(this_Atom);

                if(this_array[0] < minx) {
                    minx = this_array[0];
                }if(this_array[0] > maxx) {
                    maxx = this_array[0];
                }if(this_array[1] < miny) {
                    miny = this_array[1];
                }if(this_array[1] > maxy) {
                    maxy = this_array[1];
                }if(this_array[2] < minz) {
                    minz = this_array[2];
                }if(this_array[2] > maxz) {
                    maxz = this_array[2];
                }
            }  
        }
        if(line_string.substr(0,4) == "ATOM"){

            // std::cout << "** Entered Atom row" << std::endl;
            //Print ATOM rows
            if(!(line_string.empty())) {
                // out_file << line_string << "\n";
                std::array<double, 3> this_array = get_coords(line_string);
                std::tuple<std::string, std::string> this_data = get_data(line_string);
                Atom this_Atom = Atom(std::get<0>(this_data), std::get<1>(this_data), this_array);
                output.push_back(this_Atom);

                if(this_array[0] < minx) {
                    minx = this_array[0];
                }if(this_array[0] > maxx) {
                    maxx = this_array[0];
                }if(this_array[1] < miny) {
                    miny = this_array[1];
                }if(this_array[1] > maxy) {
                    maxy = this_array[1];
                }if(this_array[2] < minz) {
                    minz = this_array[2];
                }if(this_array[2] > maxz) {
                    maxz = this_array[2];
                }
            }
        }        
    }
    // out_file.close();

    std::tuple<std::vector<Atom>, double, double, double, double, double, double> ans {output, minx, maxx, miny, maxy, minz, maxz};
    return ans;
}

void test() {
    ;
}

//-------------------------------------------

void vectortopdb(const std::vector<Atom> &atomvector, std::string output_filename) {
    std::vector<Atom> output;
    std::ofstream out_file;
    out_file.open(output_filename);
    out_file << std::fixed;
    out_file << std::setprecision(3);

    // out_file << "REMARK  total number of initial water molecules = " << N << "\n"
            //  << "REMARK  total remaining water molecules = " << k << "\n";

    for(int i = 0; i < atomvector.size(); i++) {
        std::array<double, 3> pos = atomvector[i].getCoords();
        int z = i + 1;
        if( z > 9999) {
            z = 9999;
        }

        // out_file << "HETATM"                                   // [ 1- 6] Record
        //         << std::setw(5) << std::right << z         // [ 7-11] Serial
        //         << " "                                        // [12]    Blank
        //         << " "                                        // [13]    Space (standard for Element O)
        //         << std::setw(3) << std::left << atomvector[i].get_atomname() // [14-16] Name "O  " (padded to 3)
        //         << std::setw(3) << std::right << atomvector[i].get_resname() // [18-20] ResName
        //         << "  "                                       // [21-22] Chain ID (A) + Blank
        //         << std::setw(4) << std::right << z       // [23-26] ResSeq
        //         << "    "                                     // [27-30] Code + Blanks
        //         << std::setw(8) << std::right << pos[0]       // [31-38] X
        //         << std::setw(8) << std::right << pos[1]       // [39-46] Y
        //         << std::setw(8) << std::right << pos[2]       // [47-54] Z
        //         << std::setw(6) << std::right << "1.00"       // [55-60] Occupancy
        //         << std::setw(6) << std::right << "0.00"       // [61-66] Temp Factor
        //         << "           "                              // [67-77] Spacing
        //         << "O"                                        // [78]    Element Symbol
        //         << "\n";
        out_file << "HETATM"
         << std::setw(5) << std::right << z        // Col 7-11
         << " "                                         // Col 12
         << " O  "                                      // Col 13-16 (Atom Name)
         << " "                                         // Col 17
         << "HOH"                                       // Col 18-20 (ResName)
         << " "                                         // Col 21
         << "A"                                         // Col 22 (Chain)
         << std::setw(4) << std::right << z        // Col 23-26 (ResSeq)
         << "    "                                      // Col 27-30
         << std::setw(8) << std::fixed << std::right << pos[0] // X
         << std::setw(8) << std::fixed << std::right << pos[1] // Y
         << std::setw(8) << std::fixed << std::right << pos[2] // Z
         << std::setw(6) << "1.00"                      // Occ
         << std::setw(6) << "0.00"                      // Temp
         << "          "                                // Spacing
         << " O" << "\n";                               // Element
    }   
    out_file.close();
}


bool append_pdb_files(const std::string& filepath_1, const std::string& filepath_2, const std::string& output_path, int targetLine) {
    std::ifstream file1(filepath_1);
    std::ifstream file2(filepath_2);
    std::ofstream outFile(output_path);

    // Check if files opened successfully
    if (!file1.is_open()) {
        std::cerr << "Error: Could not open input file: " << filepath_1 << std::endl;
        return false;
    }
    if (!file2.is_open()) {
        std::cerr << "Error: Could not open input file: " << filepath_2 << std::endl;
        return false;
    }
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not create output file: " << output_path << std::endl;
        return false;
    }

    std::string line;
    int currentLine = 0;

    // 1. Read File 1 up to the target line
    // We use a loop to read line by line so we can count them
    while (currentLine < targetLine && std::getline(file1, line)) {
        outFile << line << "\n";
        currentLine++;
    }

    // 2. Inject the entirety of File 2
    // rdbuf() is a highly efficient way to copy the entire buffer of a file
    outFile << file2.rdbuf();

    // Ensure File 2 ends with a newline before resuming File 1 (optional but recommended)
    // Un-comment the next line if File 2 might be missing a terminal newline
    // outFile << "\n"; 

    // 3. Write the remainder of File 1
    // rdbuf() will dump the rest of the file starting from where the previous loop left off
    outFile << file1.rdbuf();

    // Close streams (optional, destructors will handle this automatically)
    file1.close();
    file2.close();
    outFile.close();

    return true;
}
