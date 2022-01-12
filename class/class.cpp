#include <bits/stdc++.h>
#include "class.h"
#include "../src/utils.h"

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
    primitives.shrink_to_fit();
}

Atom::Atom(vector<double> const& coord_inp, vector<Basis> const& basis_set_inp, const double mass_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description
        coord_inp        vector<double>(3)                          atomic coordinate in Bohr
        basis_set_inp    vector<Basis>(n_basis)                     basis set
        mass_inp         double                                     nuclear mass
    */
    coord = coord_inp;
    basis_set = basis_set_inp;
    n_basis = basis_set.size();
    mass = mass_inp;
    assert(mass > 0.0 && "mass should be positive");
}
Atom::~Atom(){
    // Destructor
    coord.shrink_to_fit();
    basis_set.shrink_to_fit();
}
double Atom::operator - (Atom const& atom2){
    // Get distance between two atoms
    return calcDis(coord, atom2.coord);
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
    _initialConstruct(charge_inp, atoms_inp);
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
    _initialConstruct(charge_inp, atoms_inp);
    _calcDisMat();
    _calcSItg4idx();
    _calcSItg();
    _calcTItg();
    _calcVItg();
    _calcEItg();
}
Molecule::~Molecule(){
    // Destructor
    atoms.shrink_to_fit();
    dis_mat.shrink_to_fit();
    s_itg.shrink_to_fit();
    t_itg.shrink_to_fit();
    v_itg.shrink_to_fit();
    e_itg.shrink_to_fit();
    s_itg_4idx.shrink_to_fit();
    k_factor.shrink_to_fit();
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
    double val = s_itg[idx1][idx2];
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
    double val = s_itg[idx1][idx2];
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
    if(idx2 > idx1) swap(idx1, idx2);
    if(idx4 > idx3) swap(idx3, idx4);
    double val = e_itg[idx1][idx2][idx3][idx4];
    return val;
}
void Molecule::_initialConstruct(const int charge_inp, vector<Atom> const& atoms_inp){
    // common initialization for explicit and implicit integral
    atoms = atoms_inp;
    n_atom = atoms.size();        
    assert(n_atom > 0 && "Number of atoms should be positive");
    charge = charge_inp;
    _calcNBasis();
    return;
}
void Molecule::_initializeMat2idx(vector<vector<double>>& mat){
    /*
    Initialize 2-indexed matrices: dis_mat, s_itg, t_itg, v_itg
    Input:
        Variable         Data type                                  Description
        mat              vector<vector<double>>                     matrix to be initialized
    */
    for(int i = 0 ; i < n_basis ; ++i){
        vector<double> v1;
        for(int j = 0 ; j <= i ; ++j){
            v1.push_back(0.0);
        }
        mat.push_back(v1);
        v1.shrink_to_fit();
    }
    return;
}
void Molecule::_initializeMat4idx(vector<vector<vector<vector<double>>>>& mat){
    /*
    Initialize 4-indexed matrices: s_itg_4idx, k_factor
    Input:
        Variable         Data type                                  Description
        mat              vector<vector<vector<vector<double>>>>     matrix to be initialized
    */
    int basis1_idx = 0;
    for(int atom1_idx = 0 ; atom1_idx < n_atom ; ++atom1_idx){        
        Atom atom1 = atoms[atom1_idx];        
        for(const auto& basis1: atom1.basis_set){
            vector<vector<vector<double>>> v3;
            for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
                vector<vector<double>> v2;
                int basis2_idx = 0;
                for(int atom2_idx = 0 ; atom2_idx < n_atom ; ++atom2_idx){        
                    Atom atom2 = atoms[atom2_idx];                
                    for(const auto& basis2: atom2.basis_set){                        
                        vector<double> v1;
                        for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                            v1.push_back(0);
                        }
                        v2.push_back(v1);
                        v1.shrink_to_fit();
                        ++basis2_idx;
                        if(basis2_idx > basis1_idx) break;
                    }                    
                    if(basis2_idx > basis1_idx) break;
                }
                v3.push_back(v2);
                v2.shrink_to_fit();
            }
            mat.push_back(v3);
            v3.shrink_to_fit();
            ++basis1_idx;
        }
    }
    return;
}
void Molecule::_initializeEItg(vector<vector<vector<vector<double>>>>& mat){
    /*
    Initialize e_itg
    Input:
        Variable         Data type                                  Description
        mat              vector<vector<vector<vector<double>>>>     matrix to be initialized
    */
    for(int i1 = 0 ; i1 < n_basis ; ++i1){
        vector<vector<vector<double>>> v3;
        for(int j1 = 0 ; j1 <= i1 ; ++j1){
            vector<vector<double>> v2;
            for(int i2 = 0 ; i2 <= i1 ; ++i2){
                vector<double> v1;
                for(int j2 = 0 ; j2 <= i2 ; ++j2){
                    v1.push_back(0.0);
                }
                v2.push_back(v1);
                v1.shrink_to_fit();
            }
            v3.push_back(v2);
            v2.shrink_to_fit();
        }
        mat.push_back(v3);
        v3.shrink_to_fit();
    }
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
void Molecule::_calcDisMat(){
    // calculate distance matrix
    is_dis_calculated = true;
    _initializeMat2idx(dis_mat);
    int basis1_idx = 0;
    for(int atom1_idx = 0 ; atom1_idx < n_atom ; ++atom1_idx){        
        Atom atom1 = atoms[atom1_idx];        
        for(const auto& basis1: atom1.basis_set){
            int basis2_idx = 0;
            for(int atom2_idx = 0 ; atom2_idx <= atom1_idx ; ++atom2_idx){
                Atom atom2 = atoms[atom2_idx];
                double cur_dis = atom2 - atom1;                
                for(const auto& basis2: atom2.basis_set){
                    dis_mat[basis1_idx][basis2_idx] = cur_dis;                    
                    ++basis2_idx;
                    if(basis2_idx > basis1_idx) break;
                }
                if(basis2_idx > basis1_idx) break;
            }
            ++basis1_idx;
        }        
    }
    return;
}
void Molecule::_checkDisMatCalculated(){
    assert(is_dis_calculated && "dis is not calculated, should call _calcDis() first");
    return;
}
void Molecule::_calcSItg4idx(){
    // calculate 4-indexed overlap integrals s_itg_4idx
    _checkDisMatCalculated();
    is_s_itg_4idx_calculated = true;
    _initializeMat4idx(s_itg_4idx);
    int basis1_idx = 0;
    for(int atom1_idx = 0 ; atom1_idx < n_atom ; ++atom1_idx){        
        Atom atom1 = atoms[atom1_idx];        
        for(const auto& basis1: atom1.basis_set){
            for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
                int basis2_idx = 0;
                for(int atom2_idx = 0 ; atom2_idx < n_atom ; ++atom2_idx){        
                    Atom atom2 = atoms[atom2_idx];                
                    for(const auto& basis2: atom2.basis_set){                        
                        for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                            s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx] = ;
                        }
                        ++basis2_idx;
                        if(basis2_idx > basis1_idx) break;
                    }                    
                    if(basis2_idx > basis1_idx) break;
                }
            }
            ++basis1_idx;
        }
    }
    return;
}
void Molecule::_checkSItg4idxCalculated(){
    assert(is_s_itg_4idx_calculated && "s_itg_4idx is not calculated, should call _calcSItg4idx() first");
    return;
}
void Molecule::_calcSItg(){
    // calculate overlap integrals s_itg    
    _checkDisMatCalculated();
    _checkSItg4idxCalculated();
    _initializeMat2idx(s_itg);
    return;
}
void Molecule::_calcTItg(){
    // calculate kinetic energy integrals t_itg    
    _checkDisMatCalculated();
    _checkSItg4idxCalculated();
    _initializeMat2idx(t_itg);
    return;
}
void Molecule::_calcVItg(){
    // calculate nucleus-elctron coulomb integrals t_itg       
    _checkSItg4idxCalculated();
    _initializeMat2idx(v_itg);
    return;
}
void Molecule::_calcKFactor(){
    // calculate prefactors k for two-electron integrals
    _checkDisMatCalculated();
    _initializeMat4idx(k_factor);
    return;
}
void Molecule::_calcEItg(){
    // calculate two-electron integrals
    _initializeEItg(e_itg);
    _calcKFactor();
    return;
}