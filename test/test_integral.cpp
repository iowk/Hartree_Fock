#include <bits/stdc++.h>
#include <filesystem>
#include "test_integral.h"
#include "../src/read_file.h"
#include "../src/utils.h"

using namespace std;

void testIntegrals(){
    // Test read input file with sample data
    string inp_dir_path = "input/";
    for(const auto & inp_file_path: filesystem::directory_iterator(inp_dir_path)){
        cout << "Reading " << inp_file_path.path() << endl;
        ifstream inp_file(inp_file_path.path());
        assert(inp_file.is_open() && "Cannot open file");
        Molecule mol_implicit = readMolecule(inp_file, false);
        inp_file.clear();
        inp_file.seekg(0, ios::beg);
        Molecule mol_explicit = readMolecule(inp_file, true);
        inp_file.close();
        assert(isItg2idxEqual(mol_implicit.s_itg, mol_explicit.s_itg) && "Incorrect s_itg");
        assert(isItg2idxEqual(mol_implicit.t_itg, mol_explicit.t_itg) && "Incorrect t_itg");
        assert(isItg2idxEqual(mol_implicit.v_itg, mol_explicit.v_itg) && "Incorrect v_itg");
        assert(isEItgEqual(mol_implicit.e_itg, mol_explicit.e_itg) && "Incorrect e_itg"); 
    }
    cout << "Test integrals passed" << endl;
}