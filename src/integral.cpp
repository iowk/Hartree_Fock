#include <bits/stdc++.h>
#define EPSILON 0.000000001

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
    double val = 0.5 * sqrt(M_PI/x) * erf(x);
    return val;
}
vector<double> calcRP(const double alpha1, const double alpha2, vector<double> const& coord1, vector<double> const& coord2){
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
double calcSItg4idxElement(const double zeta, const double xi, const double dis){
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
double calcTElementWithS(const double s_element, const double xi, const double dis){
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
double calcVElementWithS(const double s_element, const double mass, const double zeta, const double dis){
    /* 
    Calculate 4-indexed nucleus-electron coulomb integral with equation 31 
    Input:
        Variable         Data type                                  Description
        s_element        double                                     element of s_itg_4idx
        mass             double                                     nuclear mass
        zeta             double                                     sum of gaussian exponents
        dis              double                                     rIP
    Return:
        Variable         Data type                                  Description
        val              double                                     4-indexed kinetic energy integral
    */
    double val = -2.0 * mass * sqrt(zeta/M_PI) * _boysFunction(zeta*dis*dis) * s_element;
    return val;
}
double calcKFactorElement(const double zeta, const double xi, const double dis){
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
double calcEItgElementWithKs(const double zeta1, const double zeta2, const double k1, const double k2, const double rpq){
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