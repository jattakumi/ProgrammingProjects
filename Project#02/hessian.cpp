//
// Created by Joshua Atta-Kumi on 5/20/23.
//
#include "hessian.h"
#include <fstream>
#include <cassert>
#include <iostream>

using namespace std;

Hessian::Hessian(const char *filename)
{
    ifstream hessian(filename);

    if (!hessian.is_open()) {
        // Handle error when the file cannot be opened
        cout << "Error: Unable to open file " << filename << endl;
        return;
    }
    hessian >> natom;

    H = new double* [natom * 3];
    for(int i = 0; i < natom*3; i++)
        H[i] = new double[natom*3];

    for(int i = 0; i < natom * 3; i++) {
        for(int j = 0; j <natom; j++){
            //fscanf(hessian, "%lf %lf %lf", &H[i][3*j], &H[i][3*j+1], &H[i][3*j+2]);
            hessian >> H[i][3*j] >> H[i][3*j+1] >> H[i][3*j+2];
        }
    }
    hessian.close();
}
void Hessian::print_H(){
    for(int i=0; i<natom*3; i++) {
        for(int j=0; j<natom*3; j++) {
            printf("%14.6f", H[i][j]);
        }
        printf("\n");
    }
}
Hessian::~Hessian() {
    for (int i = 0; i < natom * 3; i++)
        delete[] H[i];
    delete H;
}