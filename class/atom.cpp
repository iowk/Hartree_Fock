#include <bits/stdc++.h>

using namespace std;
class Atom
{
    /* class Atom:
    Atom object with basis sets
    Members:
        Variable         Data type                                  Description
        coord            vector<double>(3)                          atomic coordinate in Bohr
        n_basis          int                                        number of basis functions
        basis_set        vector<Basis>(n_basis)                     basis set
    */  
    public:
    vector<double> coord = vector<double>(3);
    int n_basis;
    vector<Basis> basis_set;

    Atom(vector<double> const& , vector<Basis> const& );
    !Atom();
    double operator - (Atom const& );
};

Atom::Atom(vector<double> const& coord_inp, vector<Basis> const& basis_set_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description
        coord_inp        vector<double>(3)                          atomic coordinate in Bohr
        basis_set_inp    vector<Basis>(n_basis)                     basis set
    */
    coord = coord_inp;
    basis_set = basis_set_inp;
    n_basis = basis_set.size();
}
Atom::~Atom(){
    // Destructor
    coord.clear();
    basis_set.clear();
}
double Atom::operator - (Atom const& atom2){
    // Get distance between two atoms
    double sum = 0.0;
    for(int xyz = 0 ; xyz < 3 ; ++xyz){
        sum += (atom2.coord[xyz] - coord[xyz]) * (atom2.coord[xyz] - coord[xyz]);
    }
    return sqrt(sum);
}