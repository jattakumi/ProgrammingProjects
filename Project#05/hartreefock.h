//
// Created by Joshua Atta-Kumi on 5/21/23.
//

#ifndef PROJECT_03_HF_H
#define PROJECT_03_HF_H

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
    int e_count;
    int nmo;
    int no;
    int nv;
    Matrix S;   //Overlap matrix
    Matrix T;   //Kinetic matrix
    Matrix V;   //Nuclear attraction integrals
    Matrix H;   //core Hamiltonian matrix
    Matrix SOM;
    Vector ERI;
    Vector ERI_MO;
    Vector ioff;
    Matrix D;
    Matrix F;
    Matrix new_F;
    Vector E_p;
    Matrix C;
    int iter{};
    int iter_max{};
    double SCF{};
    double tot_E{};
    double old_SCF{};
    Matrix old_D;
    double nre{};
    double Emp2;

    Matrix tmp;
    Matrix X;
    Matrix Y;
    Matrix final_;

    double ****spinERI;
    double ****tau;
    double ****tau_t;
    Vector F_eval;
    Matrix spinFock;
    double spinEmp2;

    Matrix F_ae;
    Matrix F_mi;
    Matrix F_me;
    Matrix T1;
    Matrix old_T1;
    double ****T2;
    double ****old_T2;
    double ****W_mnij;
    double ****W_abef;
    double ****W_mbej;
    Matrix D_ia;
    double ****D_ijab;
    double E_cc;
    double old_E_cc;

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

    int transform_AO_2_MO(HF& hf); //This function is to transform from AO to MO basis
    int smart_transform_AO_2_MO(HF& hf); //N^5 smart algorithm
    int MP2_calc(HF& hf, int e_num); //This fxn calculates MP2 energy
    int smart_MP2_calc(HF& hf, int e_num); //This fxn calculates MP2 energy

    //CCSD Functions
    void transform_to_spin(HF& hf); //Transform MO basis spatial into spin-MO orbital basis
    Matrix create_spinFock(HF& hf, int e_num); //Create Fock in spin-orbital basis
    void build_cluster_amp(HF& hf, int e_num);
    void ccF_intermediates(HF& hf, int e_num);
    void build_tau(HF& hf);
    void ccW_intermediates(HF& hf, int e_num);
    void update_t_ia(HF& hf);
    void update_t_ijab(HF& hf);
    void cc_E(HF& hf);


    HF(const char *filename);                       //constructor
    ~HF();                                          //destructor
};
#endif //PROJECT_03_HF_H
