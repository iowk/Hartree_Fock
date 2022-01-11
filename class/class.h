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
    !Atom();
    double operator - (Atom const& );
};

