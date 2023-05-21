//Pre-processor
//Header files contain some functions and source codes
#include "molecule.h"
#include "hessian.h"
#include <vector>
#include "mass.h"
#include <cmath>
//#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
//Declarations
using namespace std;
//using namespace Eigen;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix;

int main() {
    //Step 1: Read the geometry from the geom.txt file
    Molecule mol("h2o_geom.txt", 0);
    cout << "Number of atoms: " << mol.natom << endl;
    cout << endl;
    cout << "Input Cartesian Coordinates: \n";
    mol.print_geometry();
    cout << endl;

    //Step 2: Read the Hessian from the hessian.txt file
    Hessian hess("h2o_hessian.txt");
    if (mol.natom != hess.natom){
        cout << "Number of atoms is different in both files" << endl;
    }
    cout << "Input Hessian Matrix: \n";
    hess.print_H();
    cout << endl;

    //Step 3: Mass-Weight the Hessian Matrix
    double m[mol.natom];
    for(int i = 0; i < mol.natom; i++)
        m[i] = mass[(int) mol.zvals[i]];

    cout << "Weight of atoms(amu): \n";
    for(int i = 0; i < mol.natom; i++)
        printf("%14.6f \n", m[i]);
    cout << endl;

    for(int i = 0; i < 3 * mol.natom; i++){
        for(int j = 0; j < 3 * mol.natom; j++){
            hess.H[i][j] /= sqrt(m[i/3]*m[j/3]);
        }
    }

    cout << "Mass Weighted Hessian Matrix: \n";
    hess.print_H();
    cout << endl;
    
    //Step 4: Diagonalize the Mass-Weighted Hessian Matrix
    matrix F(3*mol.natom, 3*mol.natom);
    for(int i=0; i<3*mol.natom; i++) {
        for(int j=0; j<3*mol.natom; j++) {
            F(i,j) = hess.H[i][j];
        }
    }
    Eigen::SelfAdjointEigenSolver<matrix>eigenSolver(F);
       matrix evals = eigenSolver.eigenvalues();
       //matrix evecs = eigenSolver.eigenvectors();

    cout << "Eigenvalues of the Hessian Matrix (hartree/amu*bohr^2): \n";
    cout << evals << endl;

    cout << endl;
    int size = evals.size();
    double w[size];

    //Step 5: Compute the Harmonic Vibrational Frequencies

    double conv = 1.0;
    // Energy conversion: 1Joule = 2.29371x10^17 Hartrees
    conv /= 2.293712278e17;

    // Atomic unit conversion: 1amu = 1.661x10^-27kg
    conv /= 1.660539066e-27;

    //Bohr conversion = 5.29x10^-11m
    conv /= pow(5.291772109e-11,2.0);

    double c = 2.99792458e10;
    double pi = 3.141592653;

    for(int i=0; i<size; i++)
        w[i] = (sqrt(conv*evals(i)))/(2*pi*c);

    cout << "Harmonic Vibrational Frequencies (cm^-1): \n";
    for(int i=0; i<size; i++)
        printf("%12.6f \n", w[i]);
    cout << endl;

    return 0;
}