
#include <string>

class Molecule
{
	public:
		int natom;
		int charge;
		int nCoords;
		double *zvals;
		double **geom;
        int electron;

		void print_geometry() const;
		void translate(double x, double y, double z) const;
		double bond(int atom1, int atom2) const;
		double angle(int atom1, int atom2, int atom3) const;
        double torsion(int atom1, int atom2, int atom3, int atom4) const;
		double unit(int cart, int atom1, int atom2) const;
		double oop(int atom1, int atom2, int atom3, int atom4) const;
        int electron_count();

    Molecule(const char *filename, int q);
    ~Molecule();
};
