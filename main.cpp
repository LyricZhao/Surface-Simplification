# include <iostream>

# include "mesh.hpp"

int main(int argc, char **argv) {
    if(argc != 4) {
        std:: cout << "Usage: simplifier <input> <output> <ratio>." << std:: endl;
        return 0;
    }
    std:: string input(argv[1]), output(argv[2]); double ratio = atof(argv[3]);
    
    Mesh mesh;
    mesh.load_from_file(input);
    mesh.simplify(ratio);
    mesh.write_into_file(output);
    return 0;
}