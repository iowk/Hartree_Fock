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
        z                double                                     nuclear charge
    */  
    public:
        vector<double> coord = vector<double>(3);
        double z;

        Atom(vector<double> const& , const int );
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
        is_s_itg_4idx_calculated    bool                            Record if s_itg_4idx calculated            
    */  
    public:
        int n_atom;
        vector<Atom> atoms;
        int charge;
        vector<Basis> basis_set;
        int n_basis;
        int n_electron;
        vector<vector<double>> s_itg;
        vector<vector<double>> t_itg;
        vector<vector<double>> v_itg;   
        vector<vector<vector<vector<double>>>> e_itg;
        vector<vector<double>> dis_mat; 
        vector<vector<vector<vector<double>>>> s_itg_4idx;
        
        Molecule(const int , vector<Atom> const& , vector<Basis> const& , const bool );
        ~Molecule();
        double getDisMat(int, int );
        double getSItg(int , int );
        double getTItg(int , int );
        double getVItg(int , int );
        double getEItg(int , int , int , int );        

    private:
        void _initialConstruct(const int, vector<Atom> const& , vector<Basis> const& );
        void _calcNElectron();
};