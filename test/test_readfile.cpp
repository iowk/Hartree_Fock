#include <bits/stdc++.h>
#include <filesystem>
#include "test_readfile.h"
#include "../src/read_file.h"

using namespace std;

void testReadInput(){
    // Test read input file with sample data
    string inp_dir_path = "input/";
    for(const auto & inp_file_path: filesystem::directory_iterator(inp_dir_path)){
        cout << "Reading " << inp_file_path.path() << endl;
        ifstream inp_file(inp_file_path.path());
        assert(inp_file.is_open() && "Cannot open file");
        Molecule mol_implicit = readMolecule(inp_file, false);        
        assert(mol_implicit.s_itg.size() == mol_implicit.n_basis && "Incorrect size for s_itg (mol_implicit)");
        assert(mol_implicit.t_itg.size() == mol_implicit.n_basis && "Incorrect size for t_itg (mol_implicit)");
        assert(mol_implicit.v_itg.size() == mol_implicit.n_basis && "Incorrect size for v_itg (mol_implicit)");
        assert(mol_implicit.e_itg.size() == mol_implicit.n_basis && "Incorrect size for e_itg (mol_implicit)");

        inp_file.clear();
        inp_file.seekg(0, ios::beg);        
        Molecule mol_explicit = readMolecule(inp_file, true);
        inp_file.close();
        assert(mol_explicit.s_itg.size() == mol_explicit.n_basis && "Incorrect size for s_itg (mol_explicit)");
        assert(mol_explicit.t_itg.size() == mol_explicit.n_basis && "Incorrect size for t_itg (mol_explicit)");
        assert(mol_explicit.v_itg.size() == mol_explicit.n_basis && "Incorrect size for v_itg (mol_explicit)");
        assert(mol_explicit.e_itg.size() == mol_explicit.n_basis && "Incorrect size for e_itg (mol_explicit)");
    }
    cout << "Test read input passed" << endl;
}
void testReadOutput(){
    // Test read output file with sample data
    string inp_dir_path = "output/";
    for(const auto & inp_file_path: filesystem::directory_iterator(inp_dir_path)){
        cout << "Reading " << inp_file_path.path() << endl;
        ifstream inp_file(inp_file_path.path());
        assert(inp_file.is_open() && "Cannot open file");
        double e_tot = readOutputEnergy(inp_file);
    }
    cout << "Test read output passed" << endl;
}