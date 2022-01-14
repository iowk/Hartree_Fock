#include <bits/stdc++.h>
#include <boost/algorithm/string.hpp>
#include "read_input.h"
#include "utils.h"

using namespace std;
using namespace boost::algorithm;

void setWordMap(map<string, string>& word_map){
    /*
    Set dictionary for input keywords
    Input:
        Variable         Data type                                  Description
        word_map         map<string, string>                        word map to be set
    */
    word_map["number of atoms"] = "n_atom";
    word_map["Atom labels, atom number, coords (Angstrom)"] = "atom_coord";
    word_map["Overall charge"] = "charge";
    word_map["Number of basis funcs"] = "n_basis";
    word_map["Maximum number of primitives"] = "skip";
    word_map["Basis set: Func no, At label, Z, Atom no //  nprim // (zeta cjk)"] = "basis";
    word_map["Overlap integrals"] = "s_itg";
    word_map["Kinetic integrals"] = "t_itg";
    word_map["Nuclear Attraction integrals"] = "v_itg";
    word_map["Two-Electron integrals"] = "e_itg";
    return;
}

Molecule readMolecule(ifstream& inp_file, const bool is_read_itg){
    /*
    Read molecule object from input file
    Input:
        Variable         Data type                                  Description
        is_read_itg      bool                                       read integrals from input file
        inp_file         ifstream                                   input file
    Return:
        Variable         Data type                                  Description
        val              Molecule                                   molecule object
    */
    int charge;
    int n_atom;
    int n_basis;
    vector<Atom> atoms;
    vector<Basis> basis_set;
    map<string, string> word_map;
    setWordMap(word_map);
    string line;
    bool is_n_atom_read = false;
    bool is_n_basis_read = false;
    bool is_atom_coord_read = false;
    bool is_charge_read = false;
    bool is_basis_read = false;
    while(getline(inp_file, line)){
        trim(line);
        if(word_map.find(line) != word_map.end()){
            if(word_map[line] == "skip") continue;
            else if(word_map[line] == "n_atom"){
                cout << "Reading number of atoms" << endl; 
                n_atom = readInt(inp_file);
                is_n_atom_read = true;
            }
            else if(word_map[line] == "n_basis"){
                cout << "Reading number of basis functions" << endl; 
                n_basis = readInt(inp_file);
                is_n_basis_read = true;
            }
            else if(word_map[line] == "charge"){
                cout << "Reading molecular charge" << endl; 
                charge = readInt(inp_file);
                is_charge_read = true;
            }
            else if(word_map[line] == "atom_coord"){ 
                assert(is_n_atom_read && "Number of atoms should be read before atomic coordinates");
                cout << "Reading atomic coordinates" << endl;               
                readAtomCoord(inp_file, atoms, n_atom);
                is_atom_coord_read = true;
            }            
            else if(word_map[line] == "basis"){
                cout << "Reading basis set" << endl; 
                assert(is_atom_coord_read && "Atomic coordinates should be read before basis set");
                assert(is_n_basis_read && "Number of basis functions should be read before basis set");
                readBasis(inp_file, basis_set, atoms, n_basis);
                is_basis_read = true;
            }
            else if(is_read_itg){
                if(word_map[line].substr(1,-1) == "_itg"){
                    assert(false && "Integrals should be put at the end of the input file");
                }
            }
            if(is_atom_coord_read && is_charge_read && is_basis_read) break; 
        }
    }
    assert(is_atom_coord_read && "Atomic coordinates not read");
    assert(is_charge_read && "Molecular charge not read");
    assert(is_basis_read && "Basis set not read");
    //cout << atoms.size() << ' ' << charge << ' ' << basis_set.size() << endl;
    Molecule mol(charge, atoms, basis_set, is_read_itg);

    if(is_read_itg){
        bool is_s_itg_read = false;
        bool is_t_itg_read = false;
        bool is_v_itg_read = false;
        bool is_e_itg_read = false;
        while(getline(inp_file, line)){
            trim(line);
            if(word_map[line] == "s_itg"){
                cout << "Reading overlap integrals" << endl;
                readItg2idx(inp_file, mol.s_itg, n_basis);
                is_s_itg_read = true;
            }
            else if(word_map[line] == "t_itg"){
                cout << "Reading kinetic energy integrals" << endl;
                readItg2idx(inp_file, mol.t_itg, n_basis);
                is_t_itg_read = true;
            }
            else if(word_map[line] == "v_itg"){
                cout << "Reading nucleus-electron coulomb integrals" << endl;
                readItg2idx(inp_file, mol.v_itg, n_basis);
                is_v_itg_read = true;
            }            
            else if(word_map[line] == "e_itg"){
                cout << "Reading 2-electron integrals" << endl;
                readEItg(inp_file, mol.e_itg, n_basis);
                is_e_itg_read = true;
            }            
        }
        assert(is_s_itg_read && "Overlap integrals not read");
        assert(is_t_itg_read && "Kinetic energy integrals not read");
        assert(is_v_itg_read && "Nucleus-electron coulomb integrals not read");
        assert(is_e_itg_read && "2-electron integrals not read");
    }
    /*cout << "S[0][0]: " << mol.s_itg[0][0] << endl;
    cout << "T[0][0]: " << mol.t_itg[0][0] << endl;
    cout << "V[0][0]: " << mol.v_itg[0][0] << endl;
    cout << "E[0][0][0][0]: " << mol.e_itg[0][0][0][0] << endl;
    cout << "S[1][0]: " << mol.s_itg[1][0] << endl;
    cout << "T[1][0]: " << mol.t_itg[1][0] << endl;
    cout << "V[1][0]: " << mol.v_itg[1][0] << endl;
    cout << "E[1][0][1][0]: " << mol.e_itg[1][0][1][0] << endl;*/ 
    atoms.clear();
    basis_set.clear();
    word_map.clear();
    cout << "Read file complete" << endl;
    return mol;
}
int readInt(ifstream& inp_file){
    /*
    Read charge from input file
    Input:
        Variable         Data type                                  Description        
        inp_file         ifstream                                   input file
    Return:
        Variable         Data type                                  Description
        number           int                                        integer to be read
    */
    string line;
    getline(inp_file, line);
    trim(line);
    return stoi(line);
}
void readAtomCoord(ifstream& inp_file, vector<Atom>& atoms, const int n_atom){
    /*
    Read atomic coordinates from input file
    Input:
        Variable         Data type                                  Description        
        inp_file         ifstream                                   input file
        atoms            vector<Atom>                               atoms to be read
        n_atom           int                                        number of atoms
    */
    string line;
    int cnt_atom = 0;
    while(getline(inp_file, line)){
        vector<string> ls;
        trim(line);
        split(ls, line, is_any_of("\t "), token_compress_on);
        assert(ls.size()==5 && "There should be 5 columns in atomic coordinates");
        vector<double> coord{angstromToBohr(stod(ls[2])), angstromToBohr(stod(ls[3])), angstromToBohr(stod(ls[4]))};
        atoms.push_back(Atom(coord, stoi(ls[1])));        
        coord.clear();
        ls.clear();    
        ++cnt_atom;
        if(cnt_atom == n_atom) break;    
    }
    return;
}
void readBasis(ifstream& inp_file, vector<Basis>& basis_set, vector<Atom> const& atoms, const int n_basis){
    /*
    Read atomic coordinates from input file
    Input:
        Variable         Data type                                  Description        
        inp_file         ifstream                                   input file
        basis_set        vector<Basis>                              basis set to be read
        atoms            vector<Atom>                               atoms
        n_basis          int                                        number of basis functions
    */
    string line;
    vector<Primitive> primitives;
    int cur_atom_idx;
    int cnt_basis = 0;
    int cnt_prim = 0;
    int n_prim;
    while(getline(inp_file, line)){
        vector<string> ls;
        trim(line);
        split(ls, line, is_any_of("\t "), token_compress_on);
        if(ls.size() == 4){
            cur_atom_idx = stoi(ls[3]) - 1;
            assert(cur_atom_idx <= atoms.size() && "Incorrect atom index");
        }            
        else if(ls.size() == 1) n_prim = stoi(ls[0]);
        else if(ls.size() == 2){
            Primitive prim_tmp(stod(ls[0]), stod(ls[1]));
            primitives.push_back(prim_tmp);
            ++cnt_prim;
            if(cnt_prim == n_prim){
                cnt_prim = 0;
                Basis basis_tmp(primitives, atoms[cur_atom_idx]);
                primitives.clear();
                basis_set.push_back(basis_tmp);
                ++cnt_basis;
            }
        }
        else assert(false && "Incorrect number of columns");        
        ls.clear();
        if(cnt_basis == n_basis) break;
    }
    return;
}
void readItg2idx(ifstream& inp_file, vector<vector<double>>& mat, const int n_basis){
    /*
    Read atomic coordinates from input file
    Input:
        Variable         Data type                                  Description        
        inp_file         ifstream                                   input file
        mat              vector<vector<double>>                     matrix to be read
        n_basis          int                                        number of basis functions
    */
    string line;
    while(getline(inp_file, line)){
        vector<string> ls;
        trim(line);
        split(ls, line, is_any_of("\t "), token_compress_on);
        mat[stoi(ls[0])-1][stoi(ls[1])-1] = stod(ls[2]);
        assert(ls.size() == 3 && "There should be 3 columns");
        assert(stoi(ls[0]) >= stoi(ls[1]) && "Incorrect integral index");
        if(stoi(ls[0]) == n_basis && stoi(ls[1]) == n_basis) break;
        ls.clear();
    }
    return;
}
void readEItg(ifstream& inp_file, vector<vector<vector<vector<double>>>>& mat, const int n_basis){
    /*
    Read atomic coordinates from input file
    Input:
        Variable         Data type                                  Description        
        inp_file         ifstream                                   input file
        mat              vector<vector<double>>                     matrix to be read
        n_basis          int                                        number of basis functions
    */
    string line;
    while(getline(inp_file, line)){
        vector<string> ls;
        trim(line);
        split(ls, line, is_any_of("\t "), token_compress_on);
        mat[stoi(ls[0])-1][stoi(ls[1])-1][stoi(ls[2])-1][stoi(ls[3])-1] = stod(ls[4]);
        assert(ls.size() == 5 && "There should be 5 columns");
        assert((stoi(ls[0]) >= stoi(ls[1]) && stoi(ls[2]) >= stoi(ls[3]) && stoi(ls[0]) >= stoi(ls[2]) && (stoi(ls[0]) != stoi(ls[2]) || stoi(ls[1]) >= stoi(ls[3]))) && "Incorrect integral index");
        if(stoi(ls[0]) == n_basis && stoi(ls[1]) == n_basis) break;
        ls.clear();
    }
    return;
}