#ifndef DRIVER_H
#define DRIVER_H

#include "rmsd.h"
#include "lattice.h"
#include <vector>

extern "C"{
#include "spglib/spglib.h"
}

using namespace std;
class Driver {
public:
  Driver(int, char **);
  ~Driver();

private:
  void check_config(const int);
  void readxyz(const char *);
  void help();

  void compare();
  void symmetry();

  void assignOne(LattXYZ *, int &, double [3][3], double [][3], int []);
  void writexyz(FILE *, int, double [3][3], double [][3], int [], LattXYZ *);

  int nlat, nmax;
  LattXYZ *one;
  std::vector<LattXYZ *> all;
};

#endif
