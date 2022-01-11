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
