#include <bits/stdc++.h>
#include <filesystem>
#include "test_final_energy.h"
#include "../src/read_file.h"
#include "../src/utils.h"
#include "../src/hf.h"

using namespace std;

void testFinalEnergy(){
    // Test read input file with sample data
    string inp_dir_path = "input/";
    string inp2_dir_path = "output/";
    string out_dir_path = "output_test/";
    for(const auto & inp_file_path: filesystem::directory_iterator(inp_dir_path)){
        string inp_file_s = inp_file_path.path().u8string();
        cout << "Reading " << inp_file_s << endl;
        ifstream inp_file(inp_file_s);
        ifstream inp2_file(inp2_dir_path + inp_file_s.substr(6, inp_file_s.size()-19) + ".input.out");
        string out_path = out_dir_path + inp_file_s.substr(6, inp_file_s.size()-19) + ".input.out";
        assert(inp_file.is_open() && "Cannot open input file");
        assert(inp2_file.is_open() && "Cannot open input file");
        Molecule mol = readMolecule(inp_file, false);
        inp_file.close();
        double e_tot_inp = hartreeFock(mol, 1, out_path);
        double e_tot_out = readOutputEnergy(inp2_file);
        inp2_file.close();
        assert(isEqual(e_tot_inp*0.1, e_tot_out*0.1) && "Incorrect final energy"); // threshold = 0.00001
    }
    cout << "Test final energy passed" << endl;
}