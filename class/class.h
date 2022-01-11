#include <bits/stdc++.h>

using namespace std;

class Primitive
{
    /* class Primitive:
    Primitive gaussian
    Members:
        Variable         Data type                                  Description
        alpha            double                                     exponent of the primitive gaussian
        c                double                                     coefficient of the primitive gaussian
    */  
    public:        
        double alpha;
        double c;

        Primitive(double const& , double const& );
        ~Primitive();
};
double getZeta(Primitive , Primitive );
double getXi(Primitive , Primitive );
class Basis
{
    /* class Basis:
    Basis function   
    
    Members:
        Variable         Data type                                  Description
        n_primitive      int                                        number of primitive gaussians
        primitives       vector<Primitive>(n_primitive)             primitive gaussians
    */  
    public:
    vector<Primitive> primitives;
    int n_primitive;

    Basis(vector<Primitive> const& );
    ~Basis();
};
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
    ~Atom();
    double operator - (Atom const& );
};
class Molecule
{
    /* class Molecule:
    Molecule object with Atoms
    Members:
        Variable         Data type                                  Description
        n_atom           int                                        number of atoms
        atoms            vector<Atom>(n_atom)                       atoms
        charge           int                                        charge
        n_basis          int                                        number of basis functions
        s_itg            vector<vector<double>>(n_basis, n_basis)   overlap integrals
        t_itg            vector<vector<double>>(n_basis, n_basis)   kinetic integrals
        v_itg            vector<vector<double>>(n_basis, n_basis)   nucleus-electron coulomb integrals
        e_itg            vector<vector<vector<vector<double>>>>(n_basis, n_basis, n_basis, n_basis)   two-electron integrals
        s_itg_4idx       vector<vector<vector<vector<double>>>>(n_basis, n_prim, n_basis, n_prim)   overlap integrals (4-indexed)
        k_factor         vector<vector<vector<vector<double>>>>(n_basis, n_prim, n_basis, n_prim)   prefactors k for two-electron integrals
        is_s_itg_4idx_calculated    bool                            Record if s_itg_4idx calculated            
    */  
    public:
    int n_atom;
    vector<Atom> atoms;
    int charge;
    int n_basis;
    vector<vector<double>> s_itg;
    vector<vector<double>> t_itg;
    vector<vector<double>> v_itg;    
    vector<vector<vector<vector<double>>>> e_itg;
    vector<vector<vector<vector<double>>>> s_itg_4idx;
    vector<vector<vector<vector<double>>>> k_factor;
    bool is_s_itg_4idx_calculated = false;

    Molecule(const int ,
    vector<Atom> const& , 
    vector<vector<double>> const& , 
    vector<vector<double>> const& , 
    vector<vector<double>> const& , 
    vector<vector<vector<vector<double>>>> const& );
    Molecule(const int , vector<Atom> const& );
    ~Molecule();
    double getSItg(int , int );
    double getTItg(int , int );
    double getVItg(int , int );
    double getEItg(int , int , int , int );

    private:
    void _initialize(const int, vector<Atom> const& );
    void _calcNBasis();
    void _calcSItg4idx();
    void _checkSItg4idxCalculated();
    void _calcSItg();
    void _calcTItg();
    void _calcVItg();
    void _calcKFactor();
    void _calcEItg();
};