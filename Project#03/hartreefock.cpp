//
// Created by Joshua Atta-Kumi on 5/21/23.
//

#include "hartreefock.h"
#include "molecule.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

#define BIGNUM 1000
#define INDEX(i,j) (i>j) ? (ioff(i)+j) : (ioff(j)+i)

HF::HF(const char *filename) {
    //Open File here
    std::ifstream input(filename);
    assert(input.good());

    input.seekg(0,std::ios_base::end);
    char ch = ' ';
    while(ch != '\n'){
        input.seekg(-2,std::ios_base::cur);
        if((int)input.tellg() <= 0 ){
            input.seekg(0);
            break;
        }
        input.get(ch);
    }

    input >> norb;
    cout << "Number of atomic orbitals: " << norb << endl;

    S.resize(norb, norb);
    T.resize(norb, norb);
    V.resize(norb, norb);
    H.resize(norb, norb);
    SOM.resize(norb, norb);
    D.resize(norb,norb);

    int M = (norb*(norb+1))/2;
    int N = (M*(M+1))/2;
    ERI.resize(N);
    ioff.resize(BIGNUM);

    //Create Fock Matrix of same dimension as core and oei matrices
    F.resize(norb,norb);

    input.close(); //Close input file
}

void HF::print_matrix(string mat_string, Matrix matrix){
    cout << endl;
    cout << mat_string;
    for(int i=0; i<matrix.rows(); i++) {
        for(int j=0; j<matrix.cols(); j++) {
            printf("%13.7f", matrix(i,j));
        }
        printf("\n");
    }
    cout << endl;
    return;
}
void HF::print_vector(string mat_string, Vector vect){
    cout << endl;
    cout << mat_string;
    for(int i=0; i<vect.size(); i++) {
        printf("%13.7f\n", vect(i));
    }
    cout << endl;
    return;
}

double HF::read_nre(const char *filename){
    std::ifstream nucl(filename);
    assert(nucl.good());
    nucl >> nre;
    cout << endl;
    printf("Nuclear Repulsion Energy: %12.15f \n", nre);
    nucl.close();
    return nre;
}

int HF::read_overlap(HF& hf, const char *filename){
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> hf.S(m-1,n-1) ) {
        hf.S(n-1,m-1) = hf.S(m-1,n-1);
    }
    oei.close(); //Close input file

    print_matrix("Overlap Integral Matrix (s): \n", hf.S);

    return 0;
}

int HF::read_kei(HF& hf, const char *filename){
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> hf.T(m-1,n-1) ) {
        hf.T(n-1,m-1) = hf.T(m-1,n-1);
    }

    oei.close(); //Close input file

    hf.print_matrix("Kinetic Energy Integral Matrix (t): \n", hf.T);

    return 0;
}

int HF::read_nai(HF& hf, const char *filename){
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> hf.V(m-1,n-1) ) {
        hf.V(n-1,m-1) = hf.V(m-1,n-1);
    }

    oei.close(); //Close input file

    hf.print_matrix("Nuclear Attraction Integral Matrix (v): \n", hf.V);

    return 0;

}

int HF::form_core(HF& hf){
    for(int i=0; i<hf.H.rows(); i++) {
        for(int j=0; j<hf.H.cols(); j++) {
            hf.H(i,j) = hf.T(i,j) + hf.V(i,j);
        }
    }
    hf.print_matrix("Core Hamiltonian Matrix (h): \n", hf.H);
    return 0;
}                          //build the core Hamiltonian matrix

int HF::read_eri(HF& hf, const char *filename){
    //Open File here
    std::ifstream eri(filename);
    assert(eri.good());

    //Read in file
    ioff(0) = 0;
    for(int n=1; n<1000; n++) {
        ioff(n) = ioff(n-1) + n;
    }

    int i, j, k, l, ij, kl, ijkl;
    double eri_val;        //Just need something to hold the value read in
    while( eri >> i >> j >> k >> l >> eri_val ) {
        i-=1;
        j-=1;
        k-=1;
        l-=1;

        ij = (i>j) ? (ioff(i) + j) : (ioff(j) + i);
        kl = (k>l) ? (ioff(k) + l) : (ioff(l) + k);
        ijkl = (ij>kl) ? (ioff(ij) + kl) : (ioff(kl) + ij);

        hf.ERI(ijkl) = eri_val;

    }
    eri.close(); //Close input file
    //hf.print_vector("ERI array: \n", hf.ERI);
    return 0;
}    //eri(two-electron repulsion integral)

//Orthogonalization of the Basis set functions
int HF::build_orthog(HF& hf){
    //To diagonalize we need to solve for eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.S);
    Matrix evc = solver.eigenvectors();     //This is a matrix nxn
    Matrix evc_T = evc.transpose();         //This will stay nxn
    Matrix evl = solver.eigenvalues();      //This is a vector nx1

    //Take one over the squareroot of the eigenvalues
    for(int i=0; i<evl.size(); i++)
        evl(i) = 1/(sqrt(evl(i)));

    //Make sure the eigenvalues are a Diagonal Matrix
    Matrix evl_D = evl.asDiagonal();     //This should be nxn
    hf.SOM = evc * evl_D * evc_T;
    hf.print_matrix("Symmetric Orthogonalization Matrix (S^1/2): \n", hf.SOM);
    return 0;
}
//Initial(Guess) Density Matrix function
int HF::build_density(HF& hf, int e_num){
    //Transform Fock matrix (F -> F')
    hf.new_F = hf.SOM.transpose() * hf.F * hf.SOM;            //Now new Fock matrix is used as guess

    //Print F' on first iteration
    if(hf.iter == 0){
        hf.print_matrix("Initial Fock Matrix (F'): \n", hf.new_F);
    }

    //Diagonalize the F'
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.new_F);
    Matrix C_p = solver.eigenvectors();   //The eigenvectors we will use in the transformation
    Matrix E = solver.eigenvalues();     //The eigenvalules - E_0 matrix containing the initial orbital energies

    //print_matrix("Eigenvectors (C' Matrix): \n", C_p0);
    //print_vector("Eigenvalues (Orbital energies): \n", E_0);

    //Transform the e-vectors into the original (non-orthogonal) AO basis
    hf.C = hf.SOM * C_p;

    //Print C Matrix on first iteration
    if(hf.iter == 0){
        hf.print_matrix("Initial Coefficient Matrix (C): \n", hf.C);
    }

    //Build Density
    //We will have calculated total # of electrons in molecule.cc
    int occ = e_num/2;             // occ is the number of doubly-occupied orbitals
    Matrix C_d = hf.C.block(0,0,hf.C.rows(),occ);
    hf.D = C_d * C_d.transpose();      //C_do is 7x5 and C_do_T is 5x7 so my resulting matrix is 7x7

    //Print density on first iteration
    if(hf.iter == 0){
        hf.print_matrix("Initial Density Matrix (D): \n", hf.D);
    }
    return 0;
}           //e_num = electron number

//Self-consistent field iteration
int HF::compute_SCF(HF& hf){
    hf.SCF = 0.0;
    for(int i=0; i<hf.D.rows(); i++) {
        for(int j=0; j<hf.D.cols(); j++) {
            hf.SCF += hf.D(i,j) * (hf.H(i,j) + hf.F(i,j));
        }
    }
    hf.tot_E = hf.SCF + hf.nre;
    return 0;

}
int HF::update_Fock(HF& hf){
    int i, j, k, l, ij, kl, ijkl, ik, jl, ikjl;
    hf.F = hf.H;                    //So the new Fock matrix comes from adding core_H to Density*TEI
    for(i=0; i<hf.F.rows(); i++) {
        for(j=0; j<hf.F.rows(); j++) {
            for(k=0; k<hf.F.rows(); k++) {
                for(l=0; l<hf.F.rows(); l++) {
                    ij = INDEX(i,j);
                    //if(i>j) ij = i*(i+1)/2 + j;
                    //else ij = j*(j+1)/2 + i;

                    kl = INDEX(k,l);
                    ijkl = INDEX(ij,kl);
                    ik = INDEX(i,k);
                    jl = INDEX(j,l);
                    ikjl = INDEX(ik,jl);

                    hf.F(i,j) += hf.D(k,l) * (2.0 * hf.ERI(ijkl) - hf.ERI(ikjl));
                }
            }
        }
    }

    if(hf.iter == 1){
        hf.print_matrix("Fock Matrix (F): \n", hf.F);
    }

    return 0;
}

HF::~HF() {

}