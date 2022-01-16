#include <bits/stdc++.h>
#include "class.h"
#include "../src/utils.h"
#include "../src/integral.h"

using namespace std;

Primitive::Primitive(double const& alpha_inp, double const& c_inp) : alpha(alpha_inp), c(c_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description        
        alpha_inp        double                                     exponent of the primitive gaussian
        c_inp            double                                     coefficient of the primitive gaussian
    */    
    alpha = alpha_inp;
    assert(alpha > 0.0 && "exponent of the primitive gaussian should be positive");
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
Atom::Atom(vector<double> const& coord_inp, const int z_inp) : coord(coord_inp), z(z_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description
        coord_inp        vector<double>(3)                          atomic coordinate in Bohr
        z_inp            int                                        nuclear charge
    */
    assert(z > 0 && "nuclear charge should be positive");
}
Atom::~Atom(){
    // Destructor
    coord.clear();
}
double Atom::operator - (Atom const& atom2){
    // Get distance between two atoms
    return calcDis(coord, atom2.coord);
}
Basis::Basis(vector<Primitive> const& primitives_inp, const Atom& atom_inp) : primitives(primitives_inp), atom(atom_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description
        primitives_inp   vector<Primitive>(n_primitive)             primitive gaussians
        atom_inp         Atom                                       atom which the basis is located on
    */
    n_primitive = primitives.size();
    assert(n_primitive > 0 && "Number of primitive gaussians should be positive");
}
Basis::~Basis(){
    // Destructor
    primitives.clear();
}
Molecule::Molecule(const int charge_inp, vector<Atom> const& atoms_inp, vector<Basis> const& basis_set_inp, const bool is_read_itg){
    /* Constructor
    Compute integrals implicitly
    Input:
        Variable         Data type                                  Description
        charge_inp       int                                        charge
        atoms_inp        vector<Atom>(n_atom)                       atoms
        basis_set_inp    vector<Basis>(n_basis)                     basis set of all basis functions
        is_read_itg      bool                                       read integral from file explicitly
    */
    _initialConstruct(charge_inp, atoms_inp, basis_set_inp);
    initializeMat2idx(s_itg, n_basis);
    initializeMat2idx(t_itg, n_basis);
    initializeMat2idx(v_itg, n_basis);
    initializeEItg(e_itg, n_basis);
    if(!is_read_itg){
        initializeMat2idx(dis_mat, n_basis);        
        calcDisMat(dis_mat, basis_set);
        initializeMat4idx(s_itg_4idx, basis_set);
        calcSItg4idx(s_itg_4idx, basis_set, dis_mat);
        calcSItg(s_itg, s_itg_4idx, basis_set);
        calcTItg(t_itg, s_itg_4idx, basis_set, dis_mat);
        calcVItg(v_itg, s_itg_4idx, basis_set, atoms);
        calcEItg(e_itg, basis_set, atoms, dis_mat);
    }   
}
Molecule::~Molecule(){
    // Destructor
    atoms.clear();
    basis_set.clear();
    dis_mat.clear();
    s_itg.clear();
    t_itg.clear();
    v_itg.clear();
    e_itg.clear();
    s_itg_4idx.clear();
}
double Molecule::getDisMat(int idx1, int idx2){
    /*
    get distance between atoms that basis idx1 and idx2 locate on
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     distance
    */
    if(idx2 > idx1) swap(idx1, idx2);
    double val = dis_mat[idx1][idx2];
    return val;
}
double Molecule::getSItg(int idx1, int idx2){
    /*
    get overlap integral value between basis idx1 and idx2
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
    get kinetic energy integral value between basis idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2);
    double val = t_itg[idx1][idx2];
    return val;
}
double Molecule::getVItg(int idx1, int idx2){
    /*
    get nucleus-electron coulomb integral value between basis idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2);
    double val = v_itg[idx1][idx2];
    return val;
}    
double Molecule::getEItg(int idx1, int idx2, int idx3, int idx4){
    /*
    get two-electron integral value between basis idx1 ~ idx4
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
    if(idx2 > idx1) swap(idx2, idx1);
    if(idx4 > idx3) swap(idx4, idx3);
    if(idx3 > idx1){
        swap(idx1, idx3);
        swap(idx2, idx4);
    }
    else if(idx3 == idx1 && idx4 > idx2) swap(idx2, idx4);
    double val = e_itg[idx1][idx2][idx3][idx4];
    return val;
}
void Molecule::_initialConstruct(const int charge_inp, vector<Atom> const& atoms_inp, vector<Basis> const& basis_set_inp){
    // common initialization for explicit and implicit integral
    charge = charge_inp;
    atoms = atoms_inp;
    n_atom = atoms.size();        
    assert(n_atom > 0 && "Number of atoms should be positive");
    basis_set = basis_set_inp;
    n_basis = basis_set.size();
    assert(n_basis > 0 && "Number of basis functions should be positive");
    _calcNElectron();
    assert(n_electron >= 0 && "Number of electrons should be non-negative");
    assert(n_electron <= n_basis*2 && "Number of electrons is large than n_basis*2");
    return;
}
void Molecule::_calcNElectron(){
    n_electron = 0;
    for(const auto& atom: atoms){
        n_electron += int(atom.z);
    }
    n_electron -= charge;
    return;
}
