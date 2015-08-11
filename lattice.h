#ifndef LATT_XYZ_H
#define LATT_XYX_H

#include <vector>
#include <string>
#include "memory.h"

class LattXYZ {
public:
  
  LattXYZ(FILE *fp, const char *);     // read xyz file
  ~LattXYZ();

  Memory *memory;

  char *fname;         // xyz file name
  char *comment;       // comment line of the xyz file
 
  int natom, ntype;    // number of atoms and atomic types in a unit cell
  int    *attyp;       // array to store atomic types
  double alat;         // lattice constant of lattice
  double latvec[3][3]; // lattice vectors, in unit of alat
  double **atpos;      // fractional coordinate for atoms in unit cell
  std::vector<std::string> names;

  int cartesian;       // 1, cartesian; 0, fractional
  int initialized;     // 0, un-initialized; 1, atomic positions initialized

  void car2dir();
  void dir2car();
  void dir2car(const double axis[][3],int ntm, double s[][3]);

private:
  int name2type(const char *);

  void display(const char *);      // method to display lattice info

  void GaussJordan(const int, const double *, double *);

  int count_words(const char *);
  double DotProd(double *, double *);
  void Cross(double *, double *, double *);
};
#endif
