//
// Created by Joshua Atta-Kumi on 5/18/23.
//
// To open and read file we need pre-processor:
#include <iostream>
#include "molecule.h"
#include "mass.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
using namespace std; // just a definition to call the ios function easily

int main(){
    // Read the data in the geometry file
    Molecule mol("acetaldehyde.dat", 0);

    cout << "Number of atoms: " << mol.natom << endl;
    cout << "Input Cartesian coordinates: \n";
    mol.print_geom();

    cout << "Interatomic Distances (bohr): \n";
    for(int i=0; i<mol.natom; i++)
        for(int j=0; j<i; j++)
            printf("%d %d %8.5f\n", i, j, mol.bond(i,j));
    cout << endl;

    cout << "Bond Angles: \n";
    for(int i=0; i<mol.natom; i++){
        for(int j=0; j<i; j++){
            for(int k=0; k<j; k++){
                if(mol.bond(i,j)<4.0 && mol.bond(j,k)<4.0){
                    printf("%2d- %2d- %2d %10.6f \n", i, j, k, mol.angle(i, j, k)*(180/acos(-1)) );
                }
            }
        }
    }
    cout << endl;
    cout << "Out-of-Plane Angles: \n";
    for(int i=0; i<mol.natom; i++){
        for(int k=0; k<mol.natom; k++){
            for(int j=0; j<mol.natom; j++){
                for(int l=0; l<j; l++){
                    if( i!=j && i!=k && i!=l && j!=k && j!=l && k!=l && mol.bond(i,k)<4.0 && mol.bond(k,j)<4.0 && mol.bond(k,l)<4.0 )
                        printf("%2d- %2d- %2d- %2d %10.6f \n", i, j, k, l, mol.oop(i,j,k,l)*(180.0/acos(-1.0)) );
                }
            }
        }
    }
    cout << endl;


    cout << "Torsional/Dihedral Angles: \n";
    for(int i=0; i<mol.natom; i++)
    {
        for(int j=0; j<i; j++)
        {
            for(int k=0; k<j; k++)
            {
                for(int l=0; l<k; l++)
                {
                    if( mol.bond(i,j)<4.0 && mol.bond(j,k)<4.0 && mol.bond(k,l)<4.0 )
                        printf("%2d- %2d- %2d- %2d %10.6f \n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
                }
            }
        }
    }
    cout << endl;



    double M = 0.0;
    for(int i=0; i<mol.natom; i++)
    {
        M += mass[(int) mol.zvals[i]];
    }

    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    double mi;

    for(int i=0; i<mol.natom; i++)
    {
        mi = mass[(int) mol.zvals[i]];
        xcm += mi*mol.geom[i][0];
        ycm += mi*mol.geom[i][1];
        zcm += mi*mol.geom[i][2];
    }

    xcm /= M;
    ycm /= M;
    zcm /= M;
    printf("Molecular Center of Mass: %12.8f %12.8f %12.8f \n", xcm, ycm, zcm);
    cout << endl;

    mol.translate(-xcm, -ycm, -zcm);

    Matrix I(3,3);
    for(int i=0; i<mol.natom; i++)
    {
        mi = mass[(int) mol.zvals[i]];
            I(0,0) += mi * (mol.geom[i][1]*mol.geom[i][1] + mol.geom[i][2]*mol.geom[i][2]);
            I(1,1) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][2]*mol.geom[i][2]);
            I(2,2) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][1]*mol.geom[i][1]);
            I(0,1) += mi * mol.geom[i][0]*mol.geom[i][1];
            I(0,2) += mi * mol.geom[i][0]*mol.geom[i][2];
            I(1,2) += mi * mol.geom[i][1]*mol.geom[i][2];
   }

    I(1,0) = I(0,1);
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);

    cout << "\nMoment of inertia tensor (amu bohr^2):\n";
    cout << I << endl;

    // find the principal moments
    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    cout << "\nPrincipal moments of inertia (amu * bohr^2):\n";
    cout << evals << endl;

    double conv = 0.529177249 * 0.529177249;
    cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
    cout << evals * conv << endl;

    conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
    cout << "\nPrincipal moments of inertia (g * cm^2):\n";
    cout << evals * conv << endl;

    // classify the rotor
    if(mol.natom == 2) cout << "\nMolecule is diatomic.\n";
    else if(evals(0) < 1e-4) cout << "\nMolecule is linear.\n";
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
        cout << "\nMolecule is a spherical top.\n";
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) > 1e-4))
        cout << "\nMolecule is an oblate symmetric top.\n";
    else if((fabs(evals(0) - evals(1)) > 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
        cout << "\nMolecule is a prolate symmetric top.\n";
    else cout << "\nMolecule is an asymmetric top.\n";

    return 0;
}