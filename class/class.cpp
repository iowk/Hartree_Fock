#include <bits/stdc++.h>
#include "class.h"

using namespace std;

Primitive::Primitive(double const& alpha_inp, double const& c_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description        
        alpha_inp        double                                     exponent of the primitive gaussian
        c_inp            double                                     coefficient of the primitive gaussian
    */    
    alpha = alpha_inp;
    c = c_inp;
}
Primitive::~Primitive(){
    // Destructor
}
double getZeta(Primitive p1, Primitive p2){
    /*
    Get zeta value from two primitives  
    Input:
        Variable         Data type                                  Description
        p1               Primitive                                  primitive gaussian 1
        p2               Primitive                                  primitive gaussian 2
    Return:
        Variable         Data type                                  Description
        zeta             double                                     zeta = p1.alpha + p2.alpha;
    */ 
    double zeta = p1.alpha + p2.alpha;
    return zeta;
}
double getXi(Primitive p1, Primitive p2){
    /*
    Get xi value from two primitives  
    Input:
        Variable         Data type                                  Description
        p1               Primitive                                  primitive gaussian 1
        p2               Primitive                                  primitive gaussian 2
    Return:
        Variable         Data type                                  Description
        xi               double                                     xi = p1.alpha + p2.alpha;
    */ 
    double xi = p1.alpha*p2.alpha / (p1.alpha+p2.alpha);
    return xi;
}

Basis::Basis(vector<Primitive> const& primitives_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description
        primitives_inp   vector<Primitive>(n_primitive)             primitive gaussians
    */
    primitives = primitives_inp;
    n_primitive = primitives.size();
}
Basis::~Basis(){
    // Destructor
    primitives.clear();
}

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
Molecule::Molecule(const int charge_inp,
vector<Atom> const& atoms_inp, 
vector<vector<double>> const& s_itg_inp, 
vector<vector<double>> const& t_itg_inp, 
vector<vector<double>> const& v_itg_inp, 
vector<vector<vector<vector<double>>>> const& e_itg_inp){
    /* Constructor
    Read integral from input file explicitly
    Input:
        Variable         Data type                                  Description
        charge_inp       int                                        charge
        atoms_inp        vector<Atom>(n_atom)                       atoms
        s_itg_inp        vector<vector<double>>                     overlap integrals (optional)
        t_itg_inp        vector<vector<double>>                     kinetic integrals (optional)
        v_itg_inp        vector<vector<double>>                     nucleus-electron columb integrals (optional)
        e_itg_inp        vector<vector<vector<vector<double>>>>     two-electron integrals (optional)
    */ 
    _initialize(charge_inp, atoms_inp);
    s_itg = s_itg_inp;
    t_itg = t_itg_inp;
    v_itg = v_itg_inp;
    e_itg = e_itg_inp;
}    
Molecule::Molecule(const int charge_inp, vector<Atom> const& atoms_inp){
    /* Constructor
    Compute integrals implicitly
    Input:
        Variable         Data type                                  Description
        charge_inp       int                                        charge
        atoms_inp        vector<Atom>(n_atom)                       atoms
    */ 
    _initialize(charge_inp, atoms_inp);
    _calcSItg4idx();
    _calcSItg();
    _calcTItg();
    _calcVItg();
    _calcEItg();
}
Molecule::~Molecule(){
    // Destructor
    atoms.clear();
    s_itg.clear();
    t_itg.clear();
    v_itg.clear();
    e_itg.clear();
    s_itg_4idx.clear();
    k_factor.clear();
}
double Molecule::getSItg(int idx1, int idx2){
    /*
    get overlap integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2);
    double val = s_itg[idx1][idx2];
    return val;
}
double Molecule::getTItg(int idx1, int idx2){
    /*
    get kinetic energy integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2);
    double val = s_itg[idx1][idx2];
    return val;
}
double Molecule::getVItg(int idx1, int idx2){
    /*
    get nucleus-electron coulomb integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2);
    double val = s_itg[idx1][idx2];
    return val;
}    
double Molecule::getEItg(int idx1, int idx2, int idx3, int idx4){
    /*
    get two-electron integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description            
        idx1             int                                        the first index
        idx2             int                                        the second index
        idx3             int                                        the third index
        idx4             int                                        the forth index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2);
    if(idx4 > idx3) swap(idx3, idx4);
    double val = e_itg[idx1][idx2][idx3][idx4];
    return val;
}
void Molecule::_initialize(const int charge_inp, vector<Atom> const& atoms_inp){
    // common initialization for explicit and implicit integral
    atoms = atoms_inp;
    n_atom = atoms.size();        
    assert(n_atom > 0 && "Number of atoms should be positive");
    charge = charge_inp;
    _calcNBasis();
    return;
}
void Molecule::_calcNBasis(){
    // calculate total number of basis functions
    n_basis = 0;
    for(const auto& atom: atoms){
        n_basis += atom.n_basis;
    }
    return;
}
void Molecule::_calcSItg4idx(){
    // calculate overlap integrals s_itg and 4-indexed overlap integral s_itg_4idx
    is_s_itg_4idx_calculated = true; // set s_itg_4idx as calculated
    for(int atom1_idx = 0 ; atom1_idx < n_atom ; ++atom1_idx){        
        Atom atom1 = atoms[atom1_idx];
        for(const auto& basis1: atom1.basis_set){
            for(int atom2_idx = 0 ; atom2_idx <= atom1_idx ; ++atom2_idx){
                Atom atom2 = atoms[atom2_idx];
                for(const auto& basis2: atom2.basis_set){
                    // TODO
                    break;
                }
            }
        }
    }
    return;
}
void Molecule::_checkSItg4idxCalculated(){
    assert(is_s_itg_4idx_calculated && "s_itg_4idx is not calculated, should call _calc_s_itg_4idx() first");
    return;
}
void Molecule::_calcSItg(){
    // calculate overlap integrals s_itg
    _checkSItg4idxCalculated();
    return;
}
void Molecule::_calcTItg(){
    _checkSItg4idxCalculated();
    return;
}
void Molecule::_calcVItg(){
    _checkSItg4idxCalculated();
    return;
}
void Molecule::_calcKFactor(){
    return;
}
void Molecule::_calcEItg(){
    _calcKFactor();
    return;
}