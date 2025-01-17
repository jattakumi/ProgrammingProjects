//
// Created by Joshua Atta-Kumi on 5/21/23.
//
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "hartreefock.h"
#include "molecule.h"
#include "mp2.h"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

int main()
{
    Molecule mol("geom.dat", 0);
    HF hf("s.dat");

    //Read and print nuclear repulsion energy
    hf.nre = hf.read_nre("enuc.dat");

    //Read and print one electron integrals
    hf.read_overlap(hf, "s.dat");
    hf.read_kei(hf, "t.dat");
    hf.read_nai(hf, "v.dat");
    hf.form_core(hf);

    //Read and print two electron integral
    hf.read_eri(hf, "eri.dat");

    //Build Orthogonalization Matrix
    hf.build_orthog(hf);

    //Build Initial Guess Density with guess Fock Matrix
    hf.F = hf.H;
    hf.build_density(hf, mol.electron_count());
    hf.old_D = hf.D;

    //Compute the Initial SCF energy
    hf.compute_SCF(hf);
    hf.old_SCF = hf.SCF;

    printf("The initial SCF electronic energy is %12.12f Hartrees.\n", hf.SCF);
    printf("The total energy (sum of the SCF electronic energy and nuclear repulsion energy) is %12.12f Hartrees.\n", hf.tot_E);

    //Build the new Density matrix and iterate Step 7-9 until convergence is reached
    double tol = 1e-12;
    hf.iter_max = 100;
    hf.iter = 0;
    double delta_E = 1.0;
    while(abs(delta_E) > tol && hf.iter < hf.iter_max) {
        hf.iter+=1;
        double rms_D = 0.0;

        //Build the new Fock matrix (F) for the SCF procedure
        hf.update_Fock(hf);
        //Build the new Density matrix (D)
        hf.build_density(hf, mol.electron_count());
        //Compute the new SCF energy
        hf.compute_SCF(hf);

        //Test the SCF electronic energy for convergence
        delta_E = hf.SCF - hf.old_SCF;
        hf.old_SCF = hf.SCF;

        //Test the root-mean-squared difference in Densities for convergence
        for(int i=0; i<hf.D.rows(); i++) {
            for(int j=0; j<hf.D.cols(); j++) {
                rms_D += (hf.D(i,j) - hf.old_D(i,j))*(hf.D(i,j) - hf.old_D(i,j));
            }
        }
        rms_D = sqrt(rms_D);
        hf.old_D = hf.D;

        //Print out iterations and corresponding SCF and total electronic energy, along with Delta E and RMS_D
        if(hf.iter == 1) {
            printf("%-12s %-20s %-20s %-20s %-20s \n", "Iter", "E(elec)", "E(tot)", "Delta(E)", "RMS(D)");
        }
        printf("%04d %20.12f %20.12f %20.12f %20.12f \n", hf.iter, hf.SCF, hf.tot_E, delta_E, rms_D);
    }
    hf.print_matrix("The molecular orbital (MO) coefficients in the non-orthogonal atomic orbital (AO) basis: \n", hf.C);
    hf.print_vector("The molecular orbital energies (\u03B5): \n", hf.E_p);

    //Transform the 2e- integrals into the MO basis
    hf.transform_AO_2_MO(hf);
    hf.smart_transform_AO_2_MO(hf);
    hf.print_matrix("The two electron integrals in the MO basis (transformed from the AO basis): \n", hf.ERI_MO);

    //Calculate the MP2 energy using the TEI in the MO basis
    hf.smart_MP2_calc(hf, mol.electron_count());

    return 0;
}