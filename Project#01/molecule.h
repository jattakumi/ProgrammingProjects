//
// Created by Joshua Atta-Kumi on 5/18/23.
//

#ifndef PROJECT_01_MOLECULE_H
#define PROJECT_01_MOLECULE_H

#include <string>

using namespace std;

class Molecule
{
public:
    int natom;
    int charge;
    int *zvals;
    double **geom;
    string point_group;

    void print_geom();
    void rotate(double phi);
    double unit(int cart, int a, int b);
    double oop(int a, int b, int c, int d);
    void translate(double x, double y, double z);
    double bond(int atom1, int atom2);
    double angle(int atom1, int atom2, int atom3);
    double torsion(int atom1, int atom2, int atom3, int atom4);

    Molecule(const char *filename, int q);
    ~Molecule();
};
#endif //PROJECT_01_MOLECULE_H
