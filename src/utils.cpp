#include <bits/stdc++.h>
#include "utils.h"
#define EPSILON 0.000000001

using namespace std;

double isEqual(double n1, double n2){
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
double angstromToBohr(double r_angstrom){
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