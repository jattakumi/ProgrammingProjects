//
// Created by Joshua Atta-Kumi on 5/20/23.
//
#include <string>

using namespace std;
class Hessian
{
public:
    int natom;
    double **H;

    void print_H();

    Hessian(const char *filename);
    ~Hessian();
};
