#include "pymol.h"

void createPyMOLSession(const std::string& pdb_file, const std::string& session_name, const std::string& structure_file) {

    if (std::filesystem::exists("temp/render.pml")) {
        std::filesystem::remove("temp/render.pml");
    }
    

    std::string input_file = pdb_file + ".pdb"; 
    std::string input_file_stem = std::filesystem::path(pdb_file).stem().string();
    std::string script_name = "temp/render.pml";

    std::ofstream pml(script_name);
    
    if(!pml) {
        std::cout << "ERROR: Could not create .pml file! Check folder permissions or disk space." << std::endl;
        return;
    }

    if (pml.is_open()) {
        pml << "load " << input_file << "\n";
        pml << "hide all\n";
        pml << "show spheres\n";
        pml << "set sphere_scale, 0.2\n";
        
        pml << "color red, b < .5 and " << input_file_stem << "\n";
        pml << "color white, b > 0.5 and " << input_file_stem << "\n";

        pml << "set seq_view, 1, " << input_file_stem << "\n";


        if(!structure_file.empty()){
            pml << "load " << structure_file << "\n";
            pml << "set seq_view, 0, " << std::filesystem::path(structure_file).stem().string() << "\n";
        }

        pml << "show seq_view\n";
        pml << "refresh\n";

        pml << "save " << session_name << ".pse\n";
        
        pml << "quit\n";
        pml.close();

        #ifdef _WIN32
            std::string pymol_path = "\"C:\\ProgramData\\pymol\\PyMOLWin.exe\"";            
            std::string cmd = pymol_path + "-cq -r " + script_name;        
        #elif __APPLE__
            std::string cmd = "/Applications/PyMOL.app/Contents/MacOS/PyMOL -cq -r " + script_name;
        #else
            std::string cmd = "pymol -cq -r " + script_name + " &";
        #endif
        std::system(cmd.c_str());
    } 


}