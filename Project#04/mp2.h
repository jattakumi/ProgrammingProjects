//
// Created by Joshua Atta-Kumi on 5/23/23.
//

#ifndef PROJECT_04_MP2_H
#define PROJECT_04_MP2_H

#include <string>
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "hartreefock.h"

//#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;  //Eigen::RowMajor makes it so that Eigen stores matrices by defaut in row-major order vs its actual default column-major order
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

class MP2{
public:
    Vector ERI;
    Vector ERI_MO;

    int update_Fock(HF& hf, MP2& mp2);
    int transform_AO_2_MO(HF& hf, MP2& mp2);

    MP2(const char *filename);
    ~MP2();
};

#endif //PROJECT_04_MP2_H
