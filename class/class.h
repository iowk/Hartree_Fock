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
class Atom
{
    /* class Atom:
    Atom object with basis sets
    Members:
        Variable         Data type                                  Description
        coord            vector<double>(3)                          atomic coordinate in Bohr
        mass             double                                     nuclear mass
    */  
    public:
        vector<double> coord = vector<double>(3);
        double mass;

        Atom(vector<double> const& , const double );
        ~Atom();
        double operator - (Atom const& );
};
class Basis
{
    /* class Basis:
    Basis function   
    
    Members:
        Variable         Data type                                  Description
        n_primitive      int                                        number of primitive gaussians
        primitives       vector<Primitive>(n_primitive)             primitive gaussians
        atom             Atom                                       atom which the basis is located on
    */  
    public:
        vector<Primitive> primitives;
        int n_primitive;
        Atom atom;

        Basis(vector<Primitive> const& , const Atom& );
        ~Basis();
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
        basis_set        vector<Basis>(n_basis)                     basis set of all basis functions
        s_itg            vector<vector<double>>(n_basis, n_basis)   overlap integrals
        t_itg            vector<vector<double>>(n_basis, n_basis)   kinetic integrals
        v_itg            vector<vector<double>>(n_basis, n_basis)   nucleus-electron coulomb integrals
        e_itg            vector<vector<vector<vector<double>>>>(n_basis, n_basis, n_basis, n_basis)   two-electron integrals
        dis_mat          vector<vector<double>>(n_basis, n_basis)   distance matrix
        s_itg_4idx       vector<vector<vector<vector<double>>>>(n_basis, n_prim, n_basis, n_prim)   overlap integrals (4-indexed)
        k_factor         vector<vector<vector<vector<double>>>>(n_basis, n_prim, n_basis, n_prim)   prefactors k for two-electron integrals
        is_s_itg_4idx_calculated    bool                            Record if s_itg_4idx calculated            
    */  
    public:
        int n_atom;
        vector<Atom> atoms;
        int charge;
        vector<Basis> basis_set;
        int n_basis;
        vector<vector<double>> s_itg;
        vector<vector<double>> t_itg;
        vector<vector<double>> v_itg;   
        vector<vector<vector<vector<double>>>> e_itg;
        vector<vector<double>> dis_mat; 
        vector<vector<vector<vector<double>>>> s_itg_4idx;
        vector<vector<vector<vector<double>>>> k_factor;
        bool is_dis_calculated = false;
        bool is_s_itg_4idx_calculated = false;
        
        Molecule(const int ,
        vector<Atom> const& , 
        vector<Basis> const& ,
        vector<vector<double>> const& , 
        vector<vector<double>> const& , 
        vector<vector<double>> const& , 
        vector<vector<vector<vector<double>>>> const& );
        Molecule(const int , vector<Atom> const& , vector<Basis> const& );
        ~Molecule();
        double getDisMat(int, int );
        double getSItg(int , int );
        double getTItg(int , int );
        double getVItg(int , int );
        double getEItg(int , int , int , int );

    private:
        void _initialConstruct(const int, vector<Atom> const& , vector<Basis> const& );
        void _initializeMat2idx(vector<vector<double>>& );
        void _initializeMat4idx(vector<vector<vector<vector<double>>>>& );
        void _initializeEItg(vector<vector<vector<vector<double>>>>& );
        void _calcNBasis();
        void _calcDisMat();
        void _checkDisMatCalculated();    
        void _calcSItg4idx();
        void _checkSItg4idxCalculated();
        void _calcSItg();    
        void _calcTItg();
        void _calcVItg();
        void _calcKFactor();
        void _calcEItg();
        double _calcEItgElement(const int , const int , const int , const int );
};