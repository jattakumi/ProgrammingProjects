//
// Created by Joshua Atta-Kumi on 5/21/23.
//

#ifndef PROJECT_03_HARTREEFOCK_H
#define PROJECT_03_HARTREEFOCK_H

#include <string>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

class HF
{
public:
    //variables used in this code
    int norb{};
    Matrix S;   //Overlap matrix
    Matrix T;   //Kinetic matrix
    Matrix V;   //Nuclear attraction integrals
    Matrix H;   //core Hamiltonian matrix
    Matrix SOM;
    Vector ERI;
    Vector ioff;
    Vector ERI_MO;
    Matrix D;
    Matrix F;
    Matrix new_F;
    Matrix C;
    Vector E_p;

    Matrix tmp;
    Matrix X;
    Matrix Y;
    Matrix final_;

    int iter{};
    int iter_max{};
    double SCF{};
    double tot_E{};
    double old_SCF{};
    Matrix old_D;
    double nre{};


    //Matrix and Vector functions;
    void print_matrix(string mat_string, Matrix matrix);
    void print_vector(string mat_string, Vector vect);

    //integral processing functions
    double read_nre(const char *filename);          //nre(nuclear repulsion energy)
    int read_overlap(HF& hf, const char *filename); //overlap integral
    int read_kei(HF& hf, const char *filename);     //kei(kinetic energy integral)
    int read_nai(HF& hf, const char *filename);     //nai(nuclear attraction integral)
    int form_core(HF& hf);                          //build the core Hamiltonian matrix
    int read_eri(HF& hf, const char *filename);     //eri(two-electron repulsion integral)

    //Orthogonalization of the Basis set functions
    int build_orthog(HF& hf);
    //Initial(Guess) Density Matrix function
    int build_density(HF& hf, int e_num);           //e_num = electron number

    //Self-consistent field iteration
    int compute_SCF(HF& hf);
    int update_Fock(HF& hf);

    //MP2 functions
    int transform_AO_2_MO(HF& hf); //This function is to transform from AO to MO basis
    int smart_transform_AO_2_MO(HF& hf); //N^5 smart algorithm
    int MP2_calc(HF& hf, int e_num); //This fxn calculates MP2 energy
    int smart_MP2_calc(HF& hf, int e_num); //This fxn calculates MP2 energy

    HF(const char *filename);                       //constructor
    ~HF();                                          //destructor
};
#endif //PROJECT_03_HARTREEFOCK_H
