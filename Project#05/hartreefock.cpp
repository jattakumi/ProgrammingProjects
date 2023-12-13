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
//#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+ (j)) : (((j)*((j)+1/2)+(i))))

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
    nmo = 2*norb;
    no = e_count;
    nv = nmo-no;
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
    ERI_MO.resize(N);

    tmp.resize((norb*(norb+1)/2),(norb*(norb+1)/2));
    X.resize(norb, norb);
    Y.resize(norb, norb);
    ioff.resize(BIGNUM);

    //Create Fock Matrix of same dimension as core and oei matrices
    F.resize(norb,norb);

    //Begin for-loop for creating 4D array (and assign its dimensions)
    spinERI = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        spinERI[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            spinERI[i][j] = new double*[nmo];
            for(int k=0; k<nmo; k++) {
                spinERI[i][j][k] = new double[nmo];
                for(int l=0; l<nmo; l++) {
                    spinERI[i][j][k][l] = 0.0; //Initialize matrix
                }
            }
        }
    }

    int e_uo = nmo - e_count; //unoccupied electrons that would be virtual
    T2 = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        T2[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            T2[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                T2[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    T2[i][j][a][b] = 0.0;
                }
            }
        }
    }

    old_T2 = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        old_T2[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            old_T2[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                old_T2[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    old_T2[i][j][a][b] = 0.0;
                }
            }
        }
    }

    //create effective 2-particle excitation operator tau and tau_p
    tau = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        tau[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            tau[i][j] = new double *[nmo];
            for(int a=0; a<nmo; a++) {
                tau[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    tau[i][j][a][b] = 0.0;
                }
            }
        }
    }

    tau_t = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        tau_t[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            tau_t[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                tau_t[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    tau_t[i][j][a][b] = 0.0;
                }
            }
        }
    }

    T1.resize(nmo, nmo);
    T1.setZero();

    old_T1.resize(nmo, nmo);
    old_T1.setZero();

    D_ia.resize(nmo, nmo);
    D_ia.setZero();

    //F Intermediates
    F_ae.resize(nmo, nmo);
    F_mi.resize(nmo, nmo);
    F_me.resize(nmo, nmo);

    F_ae.setZero();
    F_mi.setZero();
    F_me.setZero();

    //W Intermediates
    W_mnij = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        W_mnij[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            W_mnij[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                W_mnij[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    W_mnij[i][j][a][b] = 0.0;
                }
            }
        }
    }

    W_abef = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        W_abef[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            W_abef[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                W_abef[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    W_abef[i][j][a][b] = 0.0;
                }
            }
        }
    }

    W_mbej = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        W_mbej[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            W_mbej[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                W_mbej[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    W_mbej[i][j][a][b] = 0.0;
                }
            }
        }
    }

    T2 = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        T2[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            T2[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                T2[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    T2[i][j][a][b] = 0.0;
                }
            }
        }
    }

    D_ijab = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        D_ijab[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            D_ijab[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                D_ijab[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    D_ijab[i][j][a][b] = 0.0;
                }
            }
        }
    }

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
    hf.print_vector("ERI array: \n", hf.ERI);
    return 0;
}    //eri(two-electron repulsion integral)

//Orthogonalization of the Basis set functions
int HF::build_orthog(HF& hf){
    // solve for eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.S);
    Matrix evc = solver.eigenvectors();
    Matrix evc_T = evc.transpose();
    Matrix evl = solver.eigenvalues();

    for(int i=0; i<evl.size(); i++)
        evl(i) = 1/(sqrt(evl(i)));


    Matrix evl_D = evl.asDiagonal();
    hf.SOM = evc * evl_D * evc_T;
    hf.print_matrix("Symmetric Orthogonalization Matrix (S^1/2): \n", hf.SOM);
    return 0;
}

//Initial(Guess) Density Matrix function
int HF::build_density(HF& hf, int e_num){
    //Transform Fock matrix (F -> F')
    hf.new_F = hf.SOM.transpose() * hf.F * hf.SOM;

    //Print F' on first iteration
    if(hf.iter == 0){
        hf.print_matrix("Initial Fock Matrix (F'): \n", hf.new_F);
    }

    //Diagonalize the F'
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.new_F);
    Matrix C_p = solver.eigenvectors();   //The eigenvectors we will use in the transformation
    //Matrix E = solver.eigenvalues();     //The eigenvalules - E_0 matrix containing the initial orbital energies
    hf.E_p = solver.eigenvalues();

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
    int occ = e_num/2;
    Matrix C_d = hf.C.block(0,0,hf.C.rows(),occ);
    hf.D = C_d * C_d.transpose();

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
    hf.F = hf.H;
    for(i=0; i<hf.F.rows(); i++) {
        for(j=0; j<hf.F.rows(); j++) {
            for(k=0; k<hf.F.rows(); k++) {
                for(l=0; l<hf.F.rows(); l++) {
                    ij = INDEX(i,j);
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
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.F);
    hf.F_eval = solver.eigenvalues();

    return 0;
}

// Step 3: Transform the 2e integrals from AO to MO basis
int HF::transform_AO_2_MO(HF& hf)
{
    int i, j, k, l, ijkl;
    int p, q, r, s, pq, rs, pqrs;

    //Noddy code:
    for(i=0; i < norb; i++) {
        for(j=0; j <= i; j++) {
            for(k=0; k <= i; k++) {
                for(l=0; l <= (i==k ? j : k); l++,ijkl++) {
                    for(p=0; p < norb; p++) {
                        for(q=0; q < norb; q++) {
                            pq = INDEX(p,q);
                            for(r=0; r < norb; r++) {
                                for(s=0; s < norb; s++) {
                                    rs = INDEX(r,s);
                                    pqrs = INDEX(pq,rs);
                                    ERI_MO(ijkl) += hf.C(p,i) * hf.C(q,j) * hf.C(r,k) * hf.C(s,l) * hf.ERI(pqrs);
                                }
                            }
                        }
                    }

                }
            }
        }
    }
    return 0;
}

//The Smart Algorithm
int HF::smart_transform_AO_2_MO(HF& hf)
{
    int i,j,k,l, ij,kl,ijkl,klij;

    for(i=0, ij=0; i<norb; i++)
        for(j=0; j<=i; j++,ij++) {
            for(k=0, kl=0; k<norb; k++)
                for(l=0; l<=k; l++,kl++) {
                    ijkl=INDEX(ij,kl);
                    hf.X(k,l)=hf.X(l,k)=hf.ERI(ijkl);
                }
            hf.Y.setZero(norb,norb);
            hf.Y = hf.C.transpose()*hf.X;
            hf.X.setZero(norb,norb);
            hf.X = hf.Y*hf.C;

            for(k=0,kl=0; k<norb; k++)
                for(l=0; l<=k; l++, kl++)
                    hf.tmp(kl,ij) = hf.X(k,l);
        }
    hf.ERI.setZero((norb*(norb+1)/2)*((norb*(norb+1)/2)+1)/2);
    for(k=0, kl=0; k<norb; k++)
        for(l=0; l<=k; l++, kl++) {
            hf.X.setZero(norb,norb);
            hf.Y.setZero(norb,norb);
            for(i=0, ij=0; i<norb; i++)
                for(j=0; j<=i; j++,ij++)
                    hf.X(i,j) = hf.X(j,i) = hf.tmp(kl,ij);
            hf.Y.setZero(norb,norb);
            hf.Y = hf.C.transpose()*hf.X;
            hf.X.setZero(norb,norb);
            hf.X = hf.Y*hf.C;
            for(i=0, ij=0; i<norb; i++)
                for(j=0; j<=i; j++, ij++) {
                    klij = INDEX(kl,ij);
                    hf.ERI(klij) = hf.X(i,j);
                }
        }

    return 0;
}


int HF::smart_MP2_calc(HF& hf, int e_num)
{
    hf.Emp2 = 0.0;
    int ndocc = e_num/2;
    int nao = norb;
    int ia, ja, jb, ib, iajb, ibja;
    for(int i=0; i < ndocc; i++) {
        for(int a=ndocc; a<nao; a++) {
            ia = INDEX(i,a);
            for(int j=0; j<ndocc; j++) {
                ja = INDEX(j,a);
                for(int b=ndocc; b<nao; b++) {
                    jb = INDEX(j,b);
                    ib = INDEX(i,b);
                    iajb = INDEX(ia, jb);
                    ibja = INDEX(ib, ja);
                    hf.Emp2 += (hf.ERI[iajb] * (2*hf.ERI[iajb] - hf.ERI[ibja]))/(hf.E_p(i) + hf.E_p(j) - hf.E_p(a) - hf.E_p(b));
                }
            }
        }
    }
    cout << "MP2 Energy: " << hf.Emp2 << std::endl;
    return 0;
}

int HF::MP2_calc(HF& hf, int e_num)
{
    double E_mp2 = 0.0;
    int ndocc = e_num/2;
    int nao = norb;
    int ia, ja, jb, ib, iajb, ibja;

    for(int i=0; i < ndocc; i++) {
        for(int a=ndocc; a<nao; a++) {
            ia = INDEX(i,a);
            for(int j=0; j<ndocc; j++) {
                ja = INDEX(j,a);
                for(int b=ndocc; b<nao; b++) {
                    jb = INDEX(j,b);
                    ib = INDEX(i,b);
                    iajb = INDEX(ia, jb);
                    ibja = INDEX(ib, ja);
                    E_mp2 += (hf.ERI[iajb] * (2*hf.ERI[iajb] - hf.ERI[ibja]))/(hf.E_p(i) + hf.E_p(j) - hf.E_p(a) - hf.E_p(b));
                }
            }
        }
    }
    cout << "MP2 Energy: " << E_mp2 << std::endl;
    return 0;
}


//CCSD functions
void HF::transform_to_spin(HF& hf)
{
    int pr, ps, qr, qs, prqs, psqr;
    double val1, val2;
    for(int p=0; p<2*norb; p++) {
        for(int q=0; q<2*norb; q++) {
            for(int r=0; r<2*norb; r++) {
                for(int s=0; s<2*norb; s++) {
                    //Start Indexing
                    pr = INDEX(p/2, r/2);
                    qs = INDEX(q/2, s/2);
                    prqs = INDEX(pr, qs);
                    val1 = (hf.ERI[prqs])*(p%2==r%2)*(q%2==s%2);
                    ps = INDEX(p/2, s/2);
                    qr = INDEX(q/2, r/2);
                    psqr = INDEX(ps, qr);
                    val2 = (hf.ERI[psqr])*(p%2==s%2)*(q%2==r%2);
                    hf.spinERI[p][q][r][s] = val1 - val2;
                }
            }
        }
    }
    cout<<"The transformation into the spin orbital basis:"<<std::endl<<hf.spinERI<<std::endl;
    cout<<"What is?:"<<std::endl<<norb<<std::endl;
    return;
}

//Create spin-orbital Fock matrix
Matrix HF::create_spinFock(HF& hf, int e_num)
{
    hf.spinFock.setZero(2*norb,2*norb);
    Matrix H_MO(norb,norb);
    H_MO = hf.C.transpose()*hf.H*hf.C;

    //Begin for loops
    for(int p=0; p<2*norb; p++) {
        for(int q=0; q<2*norb; q++) {
            hf.spinFock(p,q) = H_MO(p/2, q/2)*(p%2 == q%2);
            for(int m=0; m<e_num; m++) {
                hf.spinFock(p,q) += hf.spinERI[p][m][q][m];
            }
            if(std::abs(hf.spinFock(p,q)) < std::pow(10,-7)) {
                hf.spinFock(p,q) = 0.0;
            }
        }
    }
    cout<<"The Fock matrix in the spin orbital basis:"<<std::endl<<hf.spinFock<<std::endl;

    return hf.spinFock;
}

//Build initial-guess cluster amplitudes
void HF::build_cluster_amp(HF& hf, int e_n)
{
    int e_uo = (2*norb) - e_n;
    cout << "What is e_n" << e_n<<endl;
    for(int i=0; i<e_n; i++) {
        for(int j=0; j<e_n; j++) {
            for(int a=0; a<e_uo; a++) {
                for(int b=0; b<e_uo; b++) {
                    hf.T2[i][j][a][b] = (hf.spinERI[i][j][a+e_n][b+e_n])/(hf.spinFock(i,i) + hf.spinFock(j,j) - hf.spinFock(a+e_n,a+e_n) - hf.spinFock(b+e_n,b+e_n));
                }
            }
        }
    }
    hf.spinEmp2 = 0.0;
    for(int i=0; i<e_n; i++) {
        for(int j=0; j<e_n; j++) {
            for(int a=0; a<e_uo; a++) {
                for(int b=0; b<e_uo; b++) {
                    hf.spinEmp2 += (hf.T2[i][j][a][b])*(hf.spinERI[i][j][a+e_n][b+e_n]);
                }
            }
        }
    }
    hf.spinEmp2*=0.25;
    cout << "The MP2 energy correction using the spin orbital basis: " << hf.spinEmp2 << "\n" << endl;
    return;
}

void HF::build_tau(HF& hf)
{
    int e_c = e_count;
    int e_uo = (2*norb) - e_c;
    int i,j,a,b;
    for (i=0; i<e_c; i++) {
        for(j=0; j<e_c; j++) {
            for(a=0; a<e_uo; a++) {
                for(b=0; b<e_uo; b++) {
                    hf.tau[i][j][a][b] = hf.T2[i][j][a][b] + hf.T1(i,a)*hf.T1(j,b) - hf.T1(i,b)*hf.T1(j,a);
                    hf.tau_t[i][j][a][b] = hf.T2[i][j][a][b] + 0.5*(hf.T1(i, a)*hf.T1(j,b) - hf.T1(i,b)*hf.T1(j,a));
                }
            }
        }
    }
    return;
}

// Coupled Cluster Intermediates
void HF::ccF_intermediates(HF& hf, int e_c)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;

    for(int a=0; a<nv; a++) {
        for(int e=0; e<nv; e++) {
            for(int m=0; m<no; m++) {
                sum1 += hf.spinFock(m,no+e)*hf.T1(m,a);

                for(int f=0; f<nv; f++) {
                    sum2 += hf.T1(m,f)*hf.spinERI[m][no+a][no+f][no+e];

                    for(int n=0; n<no; n++) {
                        sum3 += hf.tau_t[m][n][a][f] * hf.spinERI[m][n][no+e][no+f];
                    }
                }
            }
            hf.F_ae(a,e) = (1-(a==e))*hf.spinFock(no+a,no+e) - 0.5*sum1 + sum2 - 0.5*sum3;
        }
    }

    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;

    for(int m=0; m<no; m++) {
        for(int i=0; i<no; i++) {
            for(int e=0; e<nv; e++) {
                sum1 += hf.spinFock(m,no+e)*T1(i,e);

                for(int n=0; n<no; n++) {
                    sum2 += hf.T1(n, e)*hf.spinERI[m][n][i][no+e];

                    for(int f=0; f<nv; f++) {
                        sum3 += hf.tau_t[i][n][e][f] * hf.spinERI[m][n][no+e][no+f];
                    }
                }
            }
            hf.F_mi(m,i) = (1-(m==i))*hf.spinFock(m,i) + 0.5*sum1 + sum2 + 0.5*sum3;
        }
    }

    sum1 = 0.0;

    for(int m=0; m<no; m++) {
        for(int e=0; e<nv; e++) {
            for(int n=0; n<no; n++) {
                for(int f=0; f<nv; f++) {
                    sum1 += hf.T1(n,f)*hf.spinERI[m][n][no+e][no+f];
                }
            }
            hf.F_me(m,e) = hf.spinFock(m,no+e) + sum1;
        }
    }
    return;
}

void HF::ccW_intermediates(HF& hf, int e_c)
{
    //Wmnij
    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int m=0; m<no; m++) {
        for(int n=0; n<no; n++) {
            for(int i=0; i<no; i++) {
                for(int j=0; j<no; j++) {
                    for(int e=0; e<nv; e++) {
                        //So before this sum, there is a Permutation operator for ij
                        //I'm guessing this switches ij to ji
                        sum1 += hf.T1(j,e)*hf.spinERI[m][n][i][e+no];
                        sum1 -= hf.T1(i,e)*hf.spinERI[m][n][j][e+no];
                        for(int f=0; f<nv; f++) {
                            sum2 += hf.tau[i][j][e][f]*hf.spinERI[m][n][no+e][no+f];
                        }
                    }
                    hf.W_mnij[m][n][i][j] = hf.spinERI[m][n][i][j] + sum1 + 0.25*sum2;
                }
            }
        }
    }

    //Wabef
    sum1 = 0.0;
    sum2 = 0.0;
    for(int a=0; a<nv; a++) {
        for(int b=0; b<nv; b++) {
            for(int e=0; e<nv; e++) {
                for(int f=0; f<nv; f++) {
                    for(int m=0; m<no; m++) {
                        sum1 += hf.T1(m,b)*hf.spinERI[a+no][m][e+no][f+no];
                        sum1 -= hf.T1(m,a)*hf.spinERI[b+no][m][e+no][f+no];
                        for(int n=0; n<no; n++) {
                            sum2 += hf.tau[m][n][a][b]*hf.spinERI[m][n][no+e][no+f];
                        }
                    }
                     hf.W_abef[a][b][e][f] = hf.spinERI[a+no][b+no][e+no][f+no] - sum1 + 0.25*sum2;
                }
            }
        }
    }

    //Wmbej
    sum1 = 0.0;
    sum2 = 0.0;
    double sum3 = 0.0;
    for(int m=0; m<no; m++) {
        for(int b=0; b<nv; b++) {
            for(int e=0; e<nv; e++) {
                for(int j=0; j<no; j++) {
                    for(int f=0; f<nv; f++) {
                        sum1 += hf.T1(j,f)*hf.spinERI[m][b+no][e+no][f+no];
                        for(int n=0; n<no; n++) {
                            sum2 += hf.T1(n,b)*hf.spinERI[m][n][no+e][j];
                            sum3 += (0.5*hf.T2[j][n][f][b] + hf.T1(j,f)*hf.T1(n,b))*hf.spinERI[m][n][e+no][f+no];
                        }
                    }
                    hf.W_mbej[m][b][e][j] = hf.spinERI[m][b+no][e+no][j] + sum1 - sum2 - sum3;
                }
            }
        }
    }
    return;
}

//Updating CC Terms
void HF::update_t_ia(HF& hf)
{
    for(int i=0; i<no; i++) {
        for(int a=0; a<nv; a++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            double sum6 = 0.0;
            for(int e=0; e<nv; e++) {
                sum1 += T1(i,e)*F_ae(a,e);
            }
            for(int m=0; m<no; m++) {
                sum2 += T1(m,a)*F_mi(m,i);
            }
            for(int e=0; e<nv; e++) {
                for(int m=0; m<no; m++) {
                    sum3 += T2[i][m][a][e]*F_me(m,e);

                    for(int n=0; n<no; n++) {
                        sum6 += T2[n][m][a][e]*hf.spinERI[n][m][e+no][i];
                    }
                }
            }
            for(int f=0; f<nv; f++) {
                for(int n=0; n<no; n++) {
                    sum4 += T1(n,f)*hf.spinERI[n][a+no][i][f+no];
                }
            }
            for(int m=0; m<no; m++) {
                for(int e=0; e<nv; e++) {
                    for(int f=0; f<nv; f++) {
                        sum5 += T2[i][m][e][f]*hf.spinERI[m][a+no][e+no][f+no];
                    }
                }
            }
            D_ia(i,a) = hf.spinFock(i,i) - hf.spinFock(a,a);
            T1(i,a) = hf.spinFock(i,a) + sum1 - sum2 + sum3 - sum4 - 0.5*sum5 - 0.5*sum6;
            T1(i,a) /= D_ia(i,a);
        }
    }
    return;
}

void HF::update_t_ijab(HF& hf)
{
    for(int i=0; i<no; i++) {
        for(int j=0; j<no; j++) {
            for(int a=0; a<nv; a++) {
                for(int b=0; b<nv; b++) {
                    double sum1 = 0.0;
                    double sum2 = 0.0;
                    double sum3 = 0.0;
                    double sum4 = 0.0;
                    double sum5 = 0.0;
                    double sum6 = 0.0;
                    double sum7 = 0.0;
                    double sum8 = 0.0;
                    double sum9 = 0.0;
                    double sum10 = 0.0;

                    for(int e=0; e<nv; e++) {
                        //P_(ab)
                        sum1 += T2[i][j][a][e]*F_ae(b,e);
                        sum1 -= T2[i][j][b][e]*F_ae(a,e);

                        for(int m=0; m<no; m++) {
                            //P_(ab)
                            sum2 += T2[i][j][a][e]*T1(m,b)*F_ae(b,e);
                            sum2 -= T2[i][j][b][e]*T1(m,a)*F_ae(a,e);
                        }
                    }

                    for(int m=0; m<no; m++) {
                        //P_(ij)
                        sum3 += T2[i][m][a][b]*F_me(m,j);
                        sum3 -= T2[j][m][a][b]*F_me(m,i);

                        for(int e=0; e<nv; e++) {
                            //P_(ij)
                            sum4 += T2[i][m][a][b]*(T1(j,e)*F_me(m,e));
                            sum4 -= T2[j][m][a][b]*(T1(i,e)*F_me(m,e));
                        }
                    }

                    for(int m=0; m<no; m++) {
                        for(int n=0; n<no; n++) {
                            for(int e=0; e<nv; e++) {
                                for(int f=0; f<nv; f++) {
                                    //tau & W
                                    sum5 += tau[m][n][a][b]*hf.W_mnij[m][n][i][j];
                                    sum6 += tau[i][j][e][f]*hf.W_abef[a][b][e][f];
                                }
                            }
                        }
                    }

                    for(int m=0; m<no; m++) {
                        for(int e=0; e<nv; e++) {
                            //P_(ij)
                            sum7 += T2[i][m][a][e]*hf.W_mbej[m][b][e][j];
                            sum7 -= T2[j][m][a][e]*hf.W_mbej[m][b][e][i];

                            //P_(ab)
                            sum7 -= T2[i][m][b][e]*hf.W_mbej[m][a][e][j];
                            sum7 += T2[j][m][b][e]*hf.W_mbej[m][a][e][i];
                        }
                    }

                    for(int m=0; m<no; m++) {
                        for(int e=0; e<nv; e++) {
                            //P_(ij)
                            sum8 += T1(i,e)*T1(m,a)*(hf.spinERI[m][b+no][e+no][j]);
                            sum8 -= T1(j,e)*T1(m,a)*(hf.spinERI[m][b+no][e+no][i]);

                            //P_(ab)
                            sum8 -= T1(i,e)*T1(m,b)*(hf.spinERI[m][a+no][e+no][j]);
                            sum8 += T1(j,e)*T1(m,b)*(hf.spinERI[m][a+no][e+no][i]);
                        }
                    }

                    for(int e=0; e<nv; e++) {
                        for(int m=0; m<no; m++) {
                            //P_(ij)
                            sum9 += T1(i,e)*hf.spinERI[a+no][b+no][e+no][j];
                            sum9 -= T1(j,e)*hf.spinERI[a+no][b+no][e+no][i];
                        }
                    }

                    for(int m=0; m<no; m++) {
                        //P(ab)
                        sum10 += T1(m,a)*hf.spinERI[m][b][i][j];
                        sum10 -= T1(m,b)*hf.spinERI[m][a][i][j];
                    }

                    D_ijab[i][j][a][b] = hf.spinFock(i,i) + hf.spinFock(j,j) - hf.spinFock(a,a) - hf.spinFock(b,b);
                    T2[i][j][a][b] = hf.spinERI[i][j][a+no][b+no] + (sum1 - 0.5*sum2) - (sum3 + 0.5*sum4) + 0.5*sum5 + 0.5*sum6 + (sum7 - sum8) + sum9 - sum10;
                    T2[i][j][a][b] /= D_ijab[i][j][a][b];
                }
            }
        }
    }
    return;
}

//CC Energy
void HF::cc_E(HF& hf)
{
    E_cc = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    for(int i=0; i<no; i++) {
        for(int a=0; a<nv; a++) {
            sum1 += hf.spinFock(i,a)*T1(i,a);
            for(int j=0; j<no; j++) {
                for(int b=0; b<nv; b++) {
                    sum2 += hf.spinERI[i][j][a][b]*T2[i][j][a][b];
                    sum3 += hf.spinERI[i][j][a][b]*T1(i,a)*T1(j,b);
                }
            }
            E_cc += sum1 + 0.25*sum2 + 0.5*sum3;
        }
    }

    return;
}

//Delete allocated and used memory
HF::~HF()
{
    for(int i=0; i<2*norb; i++) {
        for(int j=0; j<2*norb; j++) {
            for(int k=0; k<2*norb; k++) {
                delete[] spinERI[i][j][k];
            }
            delete[] spinERI[i][j];
        }
        delete[] spinERI[i];
    }
    delete[] spinERI;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] T2[i][j][a];
            }
            delete[] T2[i][j];
        }
        delete[] T2[i];
    }
    delete[] T2;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] tau[i][j][a];
            }
            delete[] tau[i][j];
        }
        delete[] tau[i];
    }
    delete[] tau;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] tau_t[i][j][a];
            }
            delete[] tau_t[i][j];
        }
        delete[] tau_t[i];
    }
    delete[] tau_t;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] W_mnij[i][j][a];
            }
            delete[] W_mnij[i][j];
        }
        delete[] W_mnij[i];
    }
    delete[] W_mnij;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] W_abef[i][j][a];
            }
            delete[] W_abef[i][j];
        }
        delete[] W_abef[i];
    }
    delete[] W_abef;


    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] W_mbej[i][j][a];
            }
            delete[] W_mbej[i][j];
        }
        delete[] W_mbej[i];
    }
    delete[] W_mbej;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] T2[i][j][a];
            }
            delete[] T2[i][j];
        }
        delete[] T2[i];
    }
    delete[] T2;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] old_T2[i][j][a];
            }
            delete[] old_T2[i][j];
        }
        delete[] old_T2[i];
    }
    delete[] old_T2;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] D_ijab[i][j][a];
            }
            delete[] D_ijab[i][j];
        }
        delete[] D_ijab[i];
    }
    delete[] D_ijab;
}