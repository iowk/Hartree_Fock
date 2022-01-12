#include <bits/stdc++.h>
#include "utils.h"
#define EPSILON 0.000000001

using namespace std;

double isEqual(const double n1, const double n2){
    /*    
    Input:
        Variable         Data type                                  Description
        n1               double                                     double number 1
        n2               double                                     double number 2
    Return:
        Variable         Data type                                  Description
        is_equal         bool                                       check if abs(n1-n2) < epsilon
    */  
    bool is_equal = abs(n1-n2) < EPSILON;
    return is_equal;
}
double angstromToBohr(const double r_angstrom){
    /*    
    Input:
        Variable         Data type                                  Description
        r_angstrom       double                                     variable in Angstrom
    Return:
        Variable         Data type                                  Description
        r_bohr           double                                     variable in Bohr
    */  
    double r_bohr = 1.8897259886 * r_angstrom;
    return r_angstrom;
}
double calcDis(vector<double> const& coord1, vector<double> const& coord2){
    /*    
    Input:
        Variable         Data type                                  Description
        coord1           vector<double>(3)                          atomic coordinates 1
        coord2           vector<double>(3)                          atomic coordinates 2
    Return:
        Variable         Data type                                  Description
        dis              double                                     distance between two coordinates
    */  
    assert(coord1.size()==3 && "Incorrect shape for coord1");
    assert(coord2.size()==3 && "Incorrect shape for coord2");
    double sum = 0.0;
    for(int xyz = 0 ; xyz < 3 ; ++xyz){
        sum += (coord1[xyz]-coord2[xyz]) * (coord1[xyz]-coord2[xyz]);
    }
    double dis = sqrt(sum);
    return dis;
}