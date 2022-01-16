#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "utils.h"
#include "../class/class.h"
#define ERROR_THRES 0.000001

using namespace std;
using namespace Eigen;

double isEqual(const double n1, const double n2){
    /*    
    Input:
        Variable         Data type                                  Description
        n1               double                                     double number 1
        n2               double                                     double number 2
    Return:
        Variable         Data type                                  Description
                         bool                                       check if abs(n1-n2) < ERROR_THRES
    */  
    return abs(n1-n2) < ERROR_THRES;
}
bool isItg2idxEqual(vector<vector<double>> const& v1, vector<vector<double>> const& v2){
    /*    
    Input:
        Variable         Data type                                  Description
        v1               vector<vector<double>>                     integral 1
        v2               vector<vector<double>>                     integral 2
    Return:
        Variable         Data type                                  Description
                         bool                                       check all elements of integral < ERROR_THRES
    */
    int n_basis = v1.size();
    for(int i = 0 ; i < n_basis ; ++i){
        for(int j = 0 ; j <= i ; ++j){
            if(!isEqual(v1[i][j], v2[i][j])){
                /*cout << i << ' ' << j << endl;
                cout << v1[i][j] << ' ' << v2[i][j] << endl;*/
                return false;
            }
        }
    }
    return true;
}
bool isMatEqual(MatrixXd const& v1, MatrixXd const& v2){
    /*    
    Input:
        Variable         Data type                                  Description
        v1               MatrixXd                                   matrix 1
        v2               MatrixXd                                   matrix 2
    Return:
        Variable         Data type                                  Description
                         bool                                       check all elements of matrix < ERROR_THRES
    */
    for(int i = 0 ; i < v1.rows() ; ++i){
        for(int j = 0 ; j < v1.cols() ; ++j){
            if(!isEqual(v1(i,j), v2(i,j))){
                return false;
            }
        }
    }
    return true;
}
bool isEItgEqual(vector<vector<vector<vector<double>>>> const& v1, vector<vector<vector<vector<double>>>> const& v2){
    /*    
    Input:
        Variable         Data type                                  Description
        v1               vector<vector<vector<vector<double>>>>     integral 1
        v2               vector<vector<vector<vector<double>>>>     integral 2
    Return:
        Variable         Data type                                  Description
                         bool                                       check if abs(n1-n2) < epsilon
    */
    int n_basis = v1.size();
    for(int i1 = 0 ; i1 < n_basis ; ++i1){
        for(int j1 = 0 ; j1 <= i1 ; ++j1){
            for(int i2 = 0 ; i2 <= i1 ; ++i2){
                for(int j2 = 0 ; j2 <= (i1==i2 ? j1 : i2) ; ++j2){
                    if(!isEqual(v1[i1][j1][i2][j2], v2[i1][j1][i2][j2])){
                        return false;
                    }
                }
            }
        }
    }
    return true;
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
    return r_bohr;
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
