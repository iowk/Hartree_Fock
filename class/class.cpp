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
    assert(c > 0.0 && "coefficient of the primitive gaussian should be positive");
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
Atom::Atom(vector<double> const& coord_inp, const double mass_inp) : coord(coord_inp), mass(mass_inp){
    /* Constructor
    Input:
        Variable         Data type                                  Description
        coord_inp        vector<double>(3)                          atomic coordinate in Bohr
        mass_inp         double                                     nuclear mass
    */
    assert(mass > 0.0 && "mass should be positive");
}
Atom::~Atom(){
    // Destructor
    coord.shrink_to_fit();
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
}
Basis::~Basis(){
    // Destructor
    primitives.shrink_to_fit();
}
Molecule::Molecule(const int charge_inp,
vector<Atom> const& atoms_inp, 
vector<Basis> const& basis_set_inp,
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
        basis_set_inp    vector<Basis>(n_basis)                     basis set of all basis functions
        s_itg_inp        vector<vector<double>>                     overlap integrals (optional)
        t_itg_inp        vector<vector<double>>                     kinetic integrals (optional)
        v_itg_inp        vector<vector<double>>                     nucleus-electron columb integrals (optional)
        e_itg_inp        vector<vector<vector<vector<double>>>>     two-electron integrals (optional)
    */ 
    _initialConstruct(charge_inp, atoms_inp, basis_set_inp);
    s_itg = s_itg_inp;
    t_itg = t_itg_inp;
    v_itg = v_itg_inp;
    e_itg = e_itg_inp;
}    
Molecule::Molecule(const int charge_inp, vector<Atom> const& atoms_inp, vector<Basis> const& basis_set_inp){
    /* Constructor
    Compute integrals implicitly
    Input:
        Variable         Data type                                  Description
        charge_inp       int                                        charge
        atoms_inp        vector<Atom>(n_atom)                       atoms
        basis_set_inp    vector<Basis>(n_basis)                     basis set of all basis functions
    */ 
    _initialConstruct(charge_inp, atoms_inp, basis_set_inp);
    cout << "Calculating distance matrix" << endl;
    _calcDisMat();
    cout << "Calculating 4-indexed overlap integral" << endl;
    _calcSItg4idx();
    cout << "Calculating 2-indexed overlap integral" << endl;
    _calcSItg();
    cout << "Calculating kinetic energy integral" << endl;
    _calcTItg();
    cout << "Calculating nucleus-elctron coulomb integral" << endl;
    _calcVItg();
    cout << "Calculating 2-electron integral" << endl;
    _calcEItg();
}
Molecule::~Molecule(){
    // Destructor
    atoms.shrink_to_fit();
    basis_set.shrink_to_fit();
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
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        vector<vector<vector<double>>> v3;
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            vector<vector<double>> v2;               
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){                       
                vector<double> v1;
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    v1.push_back(0);
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
                for(int j2 = 0 ; j2 <= (i1==i2 ? j1 : i2) ; ++j2){
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
void Molecule::_calcDisMat(){
    // calculate distance matrix
    is_dis_calculated = true;
    _initializeMat2idx(dis_mat);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
            Basis basis2 = basis_set[basis2_idx];
            dis_mat[basis1_idx][basis2_idx] = calcDis(basis1.atom.coord, basis2.atom.coord);
        }        
    }
    return;
}
void Molecule::_checkDisMatCalculated(){
    assert(is_dis_calculated && "dis is not calculated, should call _calcDis() first");
    return;
}
void Molecule::_calcSItg4idx(){
    // calculate 4-indexed overlap integrals s_itg_4idx with equation 24-1
    _checkDisMatCalculated();
    is_s_itg_4idx_calculated = true;    
    _initializeMat4idx(s_itg_4idx);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    double s_element = calcSItg4idxElement(getZeta(prim1, prim2), getXi(prim1, prim2), getDisMat(basis1_idx, basis2_idx));
                    s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx] = s_element;
                }
            }
        }
    }
    return;
}
void Molecule::_checkSItg4idxCalculated(){
    assert(is_s_itg_4idx_calculated && "s_itg_4idx is not calculated, should call _calcSItg() first");
    return;
}
void Molecule::_calcSItg(){
    // calculate kinetic energy integrals t_itg    
    _checkDisMatCalculated();
    _checkSItg4idxCalculated();
    _initializeMat2idx(s_itg);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    s_itg[basis1_idx][basis2_idx] += prim1.c*prim2.c*s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx];
                }
            }
        }
    }
    return;
}
void Molecule::_calcTItg(){
    // calculate kinetic energy integrals t_itg    
    _checkDisMatCalculated();
    _checkSItg4idxCalculated();
    _initializeMat2idx(t_itg);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    double t_element = calcTElementWithS(s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx], getXi(prim1, prim2), getDisMat(basis1_idx, basis2_idx));
                    t_itg[basis1_idx][basis2_idx] += prim1.c*prim2.c*t_element;
                }
            }
        }
    }
    return;
}
void Molecule::_calcVItg(){
    // calculate nucleus-elctron coulomb integrals t_itg       
    _checkSItg4idxCalculated();
    _initializeMat2idx(v_itg);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    vector<double> rp = calcRP(prim1.alpha, prim2.alpha, basis1.atom.coord, basis2.atom.coord);
                    for(const auto& atom_i: atoms){
                        double v_element = calcVElementWithS(s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx], 
                        atom_i.mass, getZeta(prim1, prim2), calcDis(atom_i.coord, rp));
                        v_itg[basis1_idx][basis2_idx] += prim1.c*prim2.c*v_element;
                    }
                    rp.shrink_to_fit();
                }
            }
        }
    }
    return;
}
void Molecule::_calcKFactor(){
    // calculate prefactors k for two-electron integrals
    _checkDisMatCalculated();
    _initializeMat4idx(k_factor);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    k_factor[basis1_idx][prim1_idx][basis2_idx][prim2_idx] = 
                    calcKFactorElement(getZeta(prim1, prim2), getXi(prim1, prim2), getDisMat(basis1_idx, basis2_idx));
                }
            }
        }
    }
    return;
}
void Molecule::_calcEItg(){
    // calculate two-electron integrals
    _calcKFactor();
    _initializeEItg(e_itg);
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
            for(int basis3_idx = 0 ; basis3_idx <= basis1_idx ; ++basis3_idx){
                for(int basis4_idx = 0 ; basis4_idx <= (basis1_idx==basis3_idx ? basis2_idx : basis3_idx) ; ++basis4_idx){
                    e_itg[basis1_idx][basis2_idx][basis3_idx][basis4_idx] = 
                    _calcEItgElement(basis1_idx, basis2_idx, basis3_idx, basis4_idx);
                }
            }
        }
    } 
    return;
}
double Molecule::_calcEItgElement(const int basis1_idx, const int basis2_idx, const int basis3_idx, const int basis4_idx){
    /* 
    Calculate element of the 2-electron integral with equation 36, 37   
    Input:
        Variable         Data type                                  Description
        basis_idx1~4     int                                        indexes of bases
    Return:
        Variable         Data type                                  Description
        val              double                                     2-electron integral element
    */
    double val = 0.0;
    Basis basis1 = basis_set[basis1_idx];
    Basis basis2 = basis_set[basis2_idx];
    Basis basis3 = basis_set[basis3_idx];
    Basis basis4 = basis_set[basis4_idx];
    for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
        Primitive prim1 = basis1.primitives[prim1_idx];
        for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
            Primitive prim2 = basis2.primitives[prim2_idx];
            for(int prim3_idx = 0 ; prim3_idx < basis3.n_primitive ; ++prim3_idx){
                Primitive prim3 = basis3.primitives[prim3_idx];
                for(int prim4_idx = 0 ; prim4_idx < basis4.n_primitive ; ++prim4_idx){
                    Primitive prim4 = basis4.primitives[prim4_idx];
                    vector<double> rp1, rp2;
                    rp1 = calcRP(prim1.alpha, prim2.alpha, basis1.atom.coord, basis2.atom.coord);
                    rp2 = calcRP(prim3.alpha, prim4.alpha, basis3.atom.coord, basis4.atom.coord);
                    double rpq = calcDis(rp1, rp2);
                    double zeta1 = getZeta(prim1, prim2);
                    double zeta2 = getZeta(prim3, prim4);
                    double k1 = k_factor[basis1_idx][prim1_idx][basis2_idx][prim2_idx];
                    double k2 = k_factor[basis3_idx][prim3_idx][basis4_idx][prim4_idx];
                    val += prim1.c*prim2.c*prim3.c*prim4.c*calcEItgElementWithKs(zeta1, zeta2, k1, k2, rpq);
                    rp1.shrink_to_fit();
                    rp2.shrink_to_fit();
                }
            }
        }
    }
    return val;
}