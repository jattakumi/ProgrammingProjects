//
// Created by Joshua Atta-Kumi on 5/18/23.
//
#include "molecule.h"
#include <cstdio>
#include <cassert>
#include <fstream>

Molecule::Molecule(const char *filename, int q)
{
    charge = q;

    // Open the file
    std::ifstream is(filename);
    //assert(is.good());
    is >> natom;        // read the number of atoms from the file

    // allocate space
    zvals = new int[natom];
    geom = new double* [natom];
    for(int i = 0; i < natom; i++)
        geom[i] = new double[3];
    for(unsigned int i = 0; i < natom; i++)
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

    is.close();
}

void Molecule::print_geom() {
    for(int i = 0; i < natom; i++){
        printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
    }
}

double Molecule::bond(int i, int j)
{
    return sqrt(pow((geom[i][0]-geom[j][0]),2)+pow((geom[i][1]-geom[j][1]),2)+pow((geom[i][2]-geom[j][2]),2));
}

void Molecule::translate(double x, double y, double z){
    for(int i = 0; i < natom; i++){
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

double Molecule::unit(int cart, int a, int b){
    return -(geom[a][cart] - geom[b][cart])/bond(a,b);
}

double Molecule::angle(int i, int j, int k){
    return acos(unit(0,j,i)* unit(0, j, k)
                +unit(1,j,i)*unit(1,j,k)
                +unit(2,j,i)*unit(2,j,k));
}

double Molecule::oop(int a, int b, int c, int d){
    double ebcd_x = (unit(1,c,b)*unit(2,c,d)-unit(2,c,b)*unit(1,c,d));
    double ebcd_y = (unit(2,c,b)*unit(0,c,d)-unit(0,c,b)*unit(2,c,d));
    double ebcd_z = (unit(0,c,b)*unit(1,c,d)-unit(1,c,b)*unit(0,c,d));

    double exx = ebcd_x * unit(0,c,a);
    double eyy = ebcd_y * unit(1,c,a);
    double ezz = ebcd_z * unit(2,c,a);

    double theta = (exx + eyy + ezz)/sin(angle(b,c,d));

    if(theta < -1.0) theta = asin(-1.0);
    else if(theta > 1.0) theta = asin(1.0);
    else theta = asin(theta);

    return theta;
}

double Molecule::torsion(int a, int b, int c, int d)
{
    double eabc_x = (unit(1,b,a)*unit(2,b,c)-unit(2,b,a)*unit(1,b,c));
    double eabc_y = (unit(2,b,a)*unit(0,b,c)-unit(0,b,a)*unit(2,b,c));
    double eabc_z = (unit(0,b,a)*unit(1,b,c)-unit(1,b,a)*unit(0,b,c));

    double ebcd_x = (unit(1,c,b)*unit(2,c,d)-unit(2,c,b)*unit(1,c,d));
    double ebcd_y = (unit(2,c,b)*unit(0,c,d)-unit(0,c,b)*unit(2,c,d));
    double ebcd_z = (unit(0,c,b)*unit(1,c,d)-unit(1,c,b)*unit(0,c,d));

    double exx = eabc_x * ebcd_x;
    double eyy = eabc_y * ebcd_y;
    double ezz = eabc_z * ebcd_z;

    double tau = (exx + eyy + ezz)/(sin(angle(a,b,c))*sin(angle(b,c,d)));

    if(tau < -1.0) tau = acos(-1.0);
    else if(tau > 1.0) tau = acos(1.0);
    else tau = acos(tau);

    double cross_x = eabc_y*ebcd_z - eabc_z*ebcd_y;
    double cross_y = eabc_z*ebcd_x - eabc_x*ebcd_z;
    double cross_z = eabc_x*ebcd_y - eabc_y*ebcd_x;
    double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;
    double sign = 1.0;
    double dot = cross_x*unit(0,b,c)+cross_y*unit(1,b,c)+cross_z*unit(2,b,c);
    if(dot < 0.0) sign = -1.0;

    return tau*sign;
}

Molecule::~Molecule() {
    delete[] zvals;
    for(int i = 0; i < natom; i++){
        delete[] geom[i];
    }
    delete[] geom;
}

