#include <bits/stdc++.h>
#include "../class/class.h"
#include "../src/utils.h"
#define EPSILON 0.00001

using namespace std;

double _boysFunction(const double x){
    /* 
    Boys function (equation 32) 
    Input:
        Variable         Data type                                  Description
        x                double                                     input variable
    Return:
        Variable         Data type                                  Description
        val              double                                     return value
    */
    if(abs(x) < EPSILON) return 1.0;
    assert(x >= 0.0 && "input of Boys function should be non-negative");
    double val = 0.5 * sqrt(M_PI/x) * erf(sqrt(x));
    return val;
}
vector<double> _calcRP(const double alpha1, const double alpha2, vector<double> const& coord1, vector<double> const& coord2){
    /* 
    Calculate rp (equation 23-3)
    Input:
        Variable         Data type                                  Description
        alpha1           double                                     exponent of gaussian 1
        alpha2           double                                     exponent of gaussian 2
        coord1           vector<double>(3)                          coordinate of the atom gaussian 1 is located on
        coord2           vector<double>(3)                          coordinate of the atom gaussian 2 is located on
    Return:
        Variable         Data type                                  Description
        rp               vector<double>(3)                          rp
    */
    vector<double> rp(3);
    assert(coord1.size()==3 && "Incorrect shape for coord1");
    assert(coord2.size()==3 && "Incorrect shape for coord2");
    for(int xyz = 0 ; xyz < 3 ; ++xyz){
        rp[xyz] = (alpha1*coord1[xyz] + alpha2*coord2[xyz]) / (alpha1 + alpha2);
    }
    return rp;
}
double _calcSItg4idxElement(const double zeta, const double xi, const double dis){
    /* 
    Calculate 4-indexed overlap integral with equation 24-1   
    Input:
        Variable         Data type                                  Description
        zeta             double                                     sum of gaussian exponents
        xi               double                                     reduced mass of gaussian exponents
        dis              double                                     distance between atoms
    Return:
        Variable         Data type                                  Description
        val              double                                     4-indexed overlap integral
    */
    double val = exp(-xi*dis*dis) * pow((M_PI/zeta), (3.0/2.0));
    return val;
}
double _calcTElementWithS(const double s_element, const double xi, const double dis){
    /* 
    Calculate 4-indexed kinetic energy integral with equation 26  
    Input:
        Variable         Data type                                  Description
        s_element        double                                     element of s_itg_4idx
        xi               double                                     reduced mass of gaussian exponents
        dis              double                                     distance between atoms
    Return:
        Variable         Data type                                  Description
        val              double                                     4-indexed kinetic energy integral
    */
    double val = xi * (3.0 - 2.0*xi*dis*dis) * s_element;
    return val;
}
double _calcVElementWithS(const double s_element, const double charge, const double zeta, const double dis){
    /* 
    Calculate 5-indexed nucleus-electron coulomb integral with equation 31 
    Input:
        Variable         Data type                                  Description
        s_element        double                                     element of s_itg_4idx
        charge           double                                     nuclear charge
        zeta             double                                     sum of gaussian exponents
        dis              double                                     rIP
    Return:
        Variable         Data type                                  Description
        val              double                                     4-indexed kinetic energy integral
    */
    double val = -2.0 * charge * sqrt(zeta/M_PI) * _boysFunction(zeta*dis*dis) * s_element;
    return val;
}
double _calcKFactorElement(const double zeta, const double xi, const double dis){
    /* 
    Calculate prefactor k of the 2-electron integral with equation 35   
    Input:
        Variable         Data type                                  Description
        zeta             double                                     sum of gaussian exponents
        xi               double                                     reduced mass of gaussian exponents
        dis              double                                     distance between atoms
    Return:
        Variable         Data type                                  Description
        val              double                                     prefactor k
    */
    double val = sqrt(2.0) * pow(M_PI,5.0/4.0) / zeta * exp(-xi*dis*dis);
    return val;
}
double _calcEItgElementWithKs(const double zeta1, const double zeta2, const double k1, const double k2, const double rpq){
    /* 
    Calculate 2-electron integral element with prefactor ks (equation 37)  
    Input:
        Variable         Data type                                  Description
        zeta1            double                                     sum of gaussian exponents 1 and 2
        zeta2            double                                     sum of gaussian exponents 3 and 4
        k1               double                                     prefactor from primitives 1 and 2
        k2               double                                     prefactor from primitives 3 and 4
        rpq              double                                     rpq
    Return:
        Variable         Data type                                  Description
        val              double                                     2-electron integral element
    */
    double ro = zeta1*zeta2 / (zeta1+zeta2);
    double val = 1.0 / sqrt(zeta1+zeta2) * k1 * k2 * _boysFunction(ro*rpq*rpq);
    return val;
}
double _calcEItgElement(vector<vector<vector<vector<double>>>>& k_factor, vector<Basis> const& basis_set, const int basis1_idx, const int basis2_idx, const int basis3_idx, const int basis4_idx){
    /* 
    Calculate element of the 2-electron integral with equation 36, 37   
    Input:
        Variable         Data type                                  Description
        basis_set        vector<Basis>                              basis set
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
                    rp1 = _calcRP(prim1.alpha, prim2.alpha, basis1.atom.coord, basis2.atom.coord);
                    rp2 = _calcRP(prim3.alpha, prim4.alpha, basis3.atom.coord, basis4.atom.coord);
                    double rpq = calcDis(rp1, rp2);
                    double zeta1 = getZeta(prim1, prim2);
                    double zeta2 = getZeta(prim3, prim4);
                    double k1 = k_factor[basis1_idx][prim1_idx][basis2_idx][prim2_idx];
                    double k2 = k_factor[basis3_idx][prim3_idx][basis4_idx][prim4_idx];
                    val += prim1.c*prim2.c*prim3.c*prim4.c*_calcEItgElementWithKs(zeta1, zeta2, k1, k2, rpq);
                    rp1.clear();
                    rp2.clear();
                }
            }
        }
    }
    return val;
}
void initializeMat2idx(vector<vector<double>>& mat, const int n_basis){
    /*
    Initialize 2-indexed matrices: dis_mat, s_itg, t_itg, v_itg
    Input:
        Variable         Data type                                  Description
        mat              vector<vector<double>>                     matrix to be initialized
        n_basis          int                                        basis set size
    */
    for(int i = 0 ; i < n_basis ; ++i){
        vector<double> v1;
        for(int j = 0 ; j <= i ; ++j){
            v1.push_back(0.0);
        }
        mat.push_back(v1);
        v1.clear();
    }
    return;
}
void initializeMat4idx(vector<vector<vector<vector<double>>>>& mat, vector<Basis> const& basis_set){
    /*
    Initialize 4-indexed matrices: s_itg_4idx, k_factor
    Input:
        Variable         Data type                                  Description
        mat              vector<vector<vector<vector<double>>>>     matrix to be initialized
        basis_set        vector<Basis>                              basis set
    */     
    int n_basis = basis_set.size();
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
                v1.clear();
            }
            v3.push_back(v2);
            v2.clear();
        }
        mat.push_back(v3);
        v3.clear();
    }
    return;
}
void initializeEItg(vector<vector<vector<vector<double>>>>& mat, const int n_basis){
    /*
    Initialize e_itg
    Input:
        Variable         Data type                                  Description
        mat              vector<vector<vector<vector<double>>>>     matrix to be initialized
        n_basis          int                                        basis set size
    */
    for(int i1 = 0 ; i1 < n_basis ; ++i1){
        vector<vector<vector<double>>> v3;
        for(int j1 = 0 ; j1 <= i1 ; ++j1){
            vector<vector<double>> v2;
            for(int i2 = 0 ; i2 <= i1 ; ++i2){
                vector<double> v1;
                for(int j2 = 0 ; j2 <= (i1==i2 ? j1 : i2) ; ++j2){
                    //cout << i1 << ' ' << j1 << ' ' << i2 << ' ' << j2 << endl;
                    v1.push_back(0.0);
                }
                v2.push_back(v1);
                v1.clear();
            }
            v3.push_back(v2);
            v2.clear();
        }
        mat.push_back(v3);
        v3.clear();
    }
    return;
}
void calcDisMat(vector<vector<double>>& dis_mat, vector<Basis> const& basis_set){ 
    /*
    Calculate distance matrix
    Input:
        Variable         Data type                                  Description
        dis_mat          vector<vector<double>>                     distance matrix to be calculated
        basis_set        vector<Basis>                              basis set
    */
    int n_basis = basis_set.size();
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
            Basis basis2 = basis_set[basis2_idx];
            dis_mat[basis1_idx][basis2_idx] = calcDis(basis1.atom.coord, basis2.atom.coord);
        }        
    }
    return;
}
void calcSItg4idx(vector<vector<vector<vector<double>>>>& s_itg_4idx, vector<Basis> const& basis_set, vector<vector<double>> const& dis_mat){
    /*
    Calculate 4-indexed overlap integrals s_itg_4idx with equation 24-1
    Input:
        Variable         Data type                                  Description
        s_itg_4idx       vector<vector<vector<vector<double>>>>     4-indexed overlap integrals to be calculated
        basis_set        vector<Basis>                              basis set
        dis_mat          vector<vector<double>>                     distance matrix
    */
    int n_basis = basis_set.size();
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    double s_element = _calcSItg4idxElement(getZeta(prim1, prim2), getXi(prim1, prim2), dis_mat[basis1_idx][basis2_idx]);
                    s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx] = s_element;
                }
            }
        }
    }
    return;
}
void calcSItg(vector<vector<double>>& s_itg, vector<vector<vector<vector<double>>>> const& s_itg_4idx, vector<Basis> const& basis_set){
    /*
    Calculate overlap integrals s_itg
    Input:
        Variable         Data type                                  Description
        s_itg            vector<vector<double>>                     overlap integrals to be calculated
        s_itg_4idx       vector<vector<vector<vector<double>>>>     4-indexed overlap integrals
        basis_set        vector<Basis>                              basis set
    */  
    int n_basis = basis_set.size();
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
void calcTItg(vector<vector<double>>& t_itg, vector<vector<vector<vector<double>>>> const& s_itg_4idx, vector<Basis> const& basis_set, vector<vector<double>> const& dis_mat){
    /*
    Calculate kinetic energy integrals t_itg
    Input:
        Variable         Data type                                  Description
        t_itg            vector<vector<double>>                     kinetic energy integrals to be calculated
        s_itg_4idx       vector<vector<vector<vector<double>>>>     4-indexed overlap integrals
        basis_set        vector<Basis>                              basis set
        dis_mat          vector<vector<double>>                     distance matrix
    */
    int n_basis = basis_set.size();
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    double t_element = _calcTElementWithS(s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx], getXi(prim1, prim2), dis_mat[basis1_idx][basis2_idx]);
                    t_itg[basis1_idx][basis2_idx] += prim1.c*prim2.c*t_element;
                }
            }
        }
    }
    return;
}
void calcVItg(vector<vector<double>>& v_itg, vector<vector<vector<vector<double>>>> const& s_itg_4idx, vector<Basis> const& basis_set, vector<Atom> const& atoms){
    /*
    Calculate kinetic energy integrals v_itg
    Input:
        Variable         Data type                                  Description
        v_itg            vector<vector<double>>                     nucleus-electron coulomb integrals to be calculated
        s_itg_4idx       vector<vector<vector<vector<double>>>>     4-indexed overlap integrals
        basis_set        vector<Basis>                              basis set
        atoms            vector<Atom>                               atoms
    */
    int n_basis = basis_set.size();
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    vector<double> rp = _calcRP(prim1.alpha, prim2.alpha, basis1.atom.coord, basis2.atom.coord);
                    for(const auto& atom_i: atoms){
                        double v_element = _calcVElementWithS(s_itg_4idx[basis1_idx][prim1_idx][basis2_idx][prim2_idx], 
                        double(atom_i.z), getZeta(prim1, prim2), calcDis(atom_i.coord, rp));
                        v_itg[basis1_idx][basis2_idx] += prim1.c*prim2.c*v_element;
                    }
                    rp.clear();
                }
            }
        }
    }
    return;
}
void _calcKFactor(vector<vector<vector<vector<double>>>>& k_factor, vector<Basis> const& basis_set, vector<Atom> const& atoms, vector<vector<double>> const& dis_mat){
    /*
    Calculate prefactors k for 2-electron integrals
    Input:
        Variable         Data type                                  Description
        k_factor         vector<vector<vector<vector<double>>>>     prefactors k to be calculated for 2-electron integrals
        basis_set        vector<Basis>                              basis set
        dis_mat          vector<vector<double>>                     distance matrix
    */
    int n_basis = basis_set.size();
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        Basis basis1 = basis_set[basis1_idx];
        for(int prim1_idx = 0 ; prim1_idx < basis1.n_primitive ; ++prim1_idx){
            Primitive prim1 = basis1.primitives[prim1_idx];              
            for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
                Basis basis2 = basis_set[basis2_idx];
                for(int prim2_idx = 0 ; prim2_idx < basis2.n_primitive ; ++prim2_idx){
                    //cout << basis1_idx << ' ' << prim1_idx << ' ' << basis2_idx << ' ' << prim2_idx << endl;
                    Primitive prim2 = basis2.primitives[prim2_idx];
                    k_factor[basis1_idx][prim1_idx][basis2_idx][prim2_idx] = 
                    _calcKFactorElement(getZeta(prim1, prim2), getXi(prim1, prim2), dis_mat[basis1_idx][basis2_idx]);
                }
            }
        }
    }
    return;
}
void calcEItg(vector<vector<vector<vector<double>>>>& e_itg, vector<Basis> const& basis_set, vector<Atom> const& atoms, vector<vector<double>> const& dis_mat){
    /*
    Calculate 2-electron integrals
    Input:
        Variable         Data type                                  Description
        e_itg            vector<vector<vector<vector<double>>>>     2-electron integrals
        basis_set        vector<Basis>                              basis set
    */
    vector<vector<vector<vector<double>>>> k_factor;
    initializeMat4idx(k_factor, basis_set);
    _calcKFactor(k_factor, basis_set, atoms, dis_mat);
    int n_basis = basis_set.size();
    for(int basis1_idx = 0 ; basis1_idx < n_basis ; ++basis1_idx){
        for(int basis2_idx = 0 ; basis2_idx <= basis1_idx ; ++basis2_idx){
            for(int basis3_idx = 0 ; basis3_idx <= basis1_idx ; ++basis3_idx){
                for(int basis4_idx = 0 ; basis4_idx <= (basis1_idx==basis3_idx ? basis2_idx : basis3_idx) ; ++basis4_idx){
                    //cout << basis1_idx << ' ' << basis2_idx << ' ' << basis3_idx << ' ' << basis4_idx << endl;
                    e_itg[basis1_idx][basis2_idx][basis3_idx][basis4_idx] = 
                    _calcEItgElement(k_factor, basis_set, basis1_idx, basis2_idx, basis3_idx, basis4_idx);
                }
            }
        }
    }
    k_factor.clear();
    return;
}
