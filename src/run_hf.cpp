#include <bits/stdc++.h>
#include "read_file.h"
#include "hf.h"

using namespace std;

int main(){
    int verbose = 1; // Verbose, 0: write final energy to outputfile, 1: write final matrices to output file
    string inp_path;
    cout << "Enter input file path:" << endl;
    getline(cin, inp_path);
    ifstream inp_file(inp_path);
    assert(inp_file.is_open() && "Cannot open file");
    string s_read_itg;
    bool is_read_itg = false;
    cout << "Read integrals from input file (Y/n)?" << endl;
    getline(cin, s_read_itg);
    if(s_read_itg.size() == 0 || s_read_itg[0] == 'Y' || s_read_itg[0] == 'y'){
        is_read_itg = true;
        cout << "Integrals will be read from input file explicitly" << endl;
    }
    else{
        cout << "Integrals will be calculated implicitly" << endl;
    }
    Molecule mol = readMolecule(inp_file, is_read_itg);
    cout << "Read file complete" << endl;
    cout << "Starting Hartree Fock" << endl;
    string out_path = inp_path + ".out";
    double final_energy = hartreeFock(mol, verbose, out_path);
    return 0;
}