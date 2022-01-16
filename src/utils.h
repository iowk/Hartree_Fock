#include <bits/stdc++.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double isEqual(const double , const double );
bool isItg2idxEqual(vector<vector<double>> const& , vector<vector<double>> const& );
bool isMatEqual(MatrixXd const& , MatrixXd const& );
bool isEItgEqual(vector<vector<vector<vector<double>>>> const& , vector<vector<vector<vector<double>>>> const& );
double angstromToBohr(const double );
double calcDis(vector<double> const& , vector<double> const& );