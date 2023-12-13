//
// Created by Joshua Atta-Kumi on 5/23/23.
//

#include "mp2.h"
#include "hartreefock.h"
#include "molecule.h"

#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;  //Eigen::RowMaj    or makes it so that Eigen stores matrices by default in row-major order vs its actual default column-maj    or order
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;


MP2::MP2(const char *filename) {

}

int update_Fock(HF& hf, MP2& mp2){
    return 0;
}

int MP2::transform_AO_2_MO(HF& hf, MP2& mp2)
{
    return 0;
}

MP2::~MP2() {

}