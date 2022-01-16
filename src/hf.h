#include <Eigen/Dense>

using namespace Eigen;

double hartreeFock(Molecule , const int , const string& );
void setSMatInitial(MatrixXd& , Molecule );
void setHCoreInitial(MatrixXd& , Molecule );
void setPMatInitial(MatrixXd& );
void calcFMat(MatrixXd& , MatrixXd const& , MatrixXd const& , Molecule );
void calcPMat(MatrixXd& , MatrixXd const& , const int);
double calcEElec(MatrixXd const& , MatrixXd const& , MatrixXd const& );
double calcVNuc(Molecule );