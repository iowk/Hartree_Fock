#include <bits/stdc++.h>
#include "../class/class.h"
#include "utils.h"
#include "hf.h"
#include <Eigen/Dense>
#define MAX_ITER 100000

using namespace std;
using namespace Eigen;

double hartreeFock(Molecule mol, const int verbose, const string& out_path){
    /*    
    Input:
        Variable         Data type                                  Description
        mol              Molecule                                   molecule object with integrals
        verbose          int                                        Verbose, 0: write final energy to outputfile, 1: write final matrices to output file
        out_path         string                                     output file path
    Return:
        Variable         Data type                                  Description
                         double                                     final total energy
    */
    ofstream out_file(out_path);
    cout << fixed ;
    cout << setprecision(6);
    out_file << fixed ;
    out_file << setprecision(6);
    MatrixXd s_mat(mol.n_basis, mol.n_basis); // Overlap matrix
    MatrixXd sis_mat(mol.n_basis, mol.n_basis); // Inversed square root of the overlap matrix
    MatrixXd hcore_mat(mol.n_basis, mol.n_basis); // Core hamiltonian
    MatrixXd f_mat(mol.n_basis, mol.n_basis); // Fock matrix
    MatrixXd f_mat_orth(mol.n_basis, mol.n_basis); // Orthogonalized Fock matrix
    MatrixXd p_mat(mol.n_basis, mol.n_basis); // Density matrix
    MatrixXd pprev_mat(mol.n_basis, mol.n_basis); // Previous density matrix
    MatrixXd c_mat(mol.n_basis, mol.n_basis); // Coefficient matrix
    setSMatInitial(s_mat, mol);
    setHCoreInitial(hcore_mat, mol);
    f_mat = hcore_mat;
    setPMatInitial(p_mat);
    for(int cur_iter = 1 ; cur_iter <= MAX_ITER ; ++cur_iter){
        pprev_mat = p_mat;
        calcFMat(f_mat, hcore_mat, p_mat, mol);
        SelfAdjointEigenSolver<MatrixXd> sissolver(s_mat);
        sis_mat = sissolver.operatorInverseSqrt();
        f_mat_orth = sis_mat.transpose()*f_mat*sis_mat;
        
        SelfAdjointEigenSolver<MatrixXd> eigensolver(f_mat_orth);
        if(eigensolver.info() != Success) abort();
        c_mat = sis_mat*eigensolver.eigenvectors();
        calcPMat(p_mat, c_mat, mol.n_electron/2);
        if(isMatEqual(p_mat, pprev_mat)){
            cout << "Converged at iteration " << cur_iter << endl;
            out_file << "Converged at iteration " << cur_iter << endl;
            double e_elec = calcEElec(p_mat, hcore_mat, f_mat);
            double v_nuc = calcVNuc(mol);
            if(verbose > 0){
                out_file << "Final electronic energy:" << endl << e_elec << endl;
                out_file << "Final total energy:" << endl << e_elec + v_nuc << endl;
                out_file << "Converged orbital eigenvalues:" << endl << eigensolver.eigenvalues() << endl;
                out_file << "Converged orbital coefficients (orbital by orbital):" << endl << c_mat << endl;
                out_file << "Final Fock matrix:" << endl << f_mat << endl;
                out_file << "Final Fock matrix in orthogonal basis:" << endl << f_mat_orth << endl;
                out_file << "Final density matrix:" << endl << p_mat << endl;
            }
            cout << "Final energy:" << endl << e_elec + v_nuc << endl;
            out_file << "Final energy:" << endl << e_elec + v_nuc << endl;
            out_file.close();
            return e_elec + v_nuc;
        }
        assert(cur_iter < MAX_ITER && "Cannot converge");
    }
    return 0.0;
}
void setSMatInitial(MatrixXd& s_mat, Molecule mol){
    int n_basis = s_mat.rows();
    for(int i = 0 ; i < n_basis ; ++i){
        for(int j = 0 ; j < n_basis ; ++j){
            s_mat(i,j) = mol.getSItg(i, j);
        }
    }
    return;
}
void setHCoreInitial(MatrixXd& hcore_mat, Molecule mol){
    int n_basis = hcore_mat.rows();
    for(int i = 0 ; i < n_basis ; ++i){
        for(int j = 0 ; j < n_basis ; ++j){
            hcore_mat(i,j) = mol.getVItg(i, j) + mol.getTItg(i, j);
        }
    }
    return;
}
void setPMatInitial(MatrixXd& p_mat){
    int n_basis = p_mat.rows();
    for(int i = 0 ; i < n_basis ; ++i){
        for(int j = 0 ; j < n_basis ; ++j){
            p_mat(i,j) = 0.0;
        }
    }
    return;
}
void calcFMat(MatrixXd& f_mat, MatrixXd const& hcore_mat, MatrixXd const& p_mat, Molecule mol){
    int n_basis = p_mat.rows();
    for(int i1 = 0 ; i1 < n_basis ; ++i1){
        for(int i2 = 0 ; i2 < n_basis ; ++i2){
            f_mat(i1,i2) = hcore_mat(i1,i2);
            for(int i3 = 0 ; i3 < n_basis ; ++i3){
                for(int i4 = 0 ; i4 < n_basis ; ++i4){
                    f_mat(i1,i2) += p_mat(i3,i4) * (2*mol.getEItg(i1,i2,i3,i4) - mol.getEItg(i1,i3,i2,i4));
                }
            }
        }
    }
}
void calcPMat(MatrixXd& p_mat, MatrixXd const& c_mat, const int no){
    int n_basis = p_mat.rows();
    for(int i1 = 0 ; i1 < n_basis ; ++i1){
        for(int i2 = 0 ; i2 < n_basis ; ++i2){
            p_mat(i1,i2) = 0.0;
            for(int i3 = 0 ; i3 < no ; ++i3){
                p_mat(i1,i2) += c_mat(i1,i3) * c_mat(i2,i3);
            }
        }
    }
}
double calcEElec(MatrixXd const& p_mat, MatrixXd const& hcore_mat, MatrixXd const& f_mat){
    double val = 0.0;
    int n_basis = p_mat.rows();
    for(int i1 = 0 ; i1 < n_basis ; ++i1){
        for(int i2 = 0 ; i2 < n_basis ; ++i2){
            val += p_mat(i1,i2) * ((hcore_mat(i1,i2)) + f_mat(i1,i2));
        }
    }
    return val;
}
double calcVNuc(Molecule mol){
    double val = 0.0;
    for(int i1 = 0 ; i1 < mol.n_atom-1 ; ++i1){
        for(int i2 = i1+1 ; i2 < mol.n_atom ; ++i2){
            val += mol.atoms[i1].z * mol.atoms[i2].z / calcDis(mol.atoms[i2].coord, mol.atoms[i1].coord);
        }
    }
    return val;
}