#include "lattice.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <map>

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* ----------------------------------------------------------------------
   constructor does nothing
------------------------------------------------------------------------- */
LattXYZ::LattXYZ(FILE *fp, const char *file)
{
  initialized = 0;
  natom = ntype = 0;

  alat  = 1.;
  attyp = NULL;
  atpos = NULL;
  memory = NULL;
  comment = NULL;

  names.clear();

  fname = new char[strlen(file)+1];
  strcpy(fname, file);

  memory = new Memory();

  // read number of atoms
  char str[MAXLINE];
  if (fgets(str, MAXLINE, fp) == NULL) return;
  char *ptr = strtok(str," \n\t\r\f");
  if (ptr == NULL) return;
  natom = atoi(ptr);
  if (natom < 1) return;

  // allocate memory
  attyp = memory->create(attyp,natom,"attyp");
  atpos = memory->create(atpos,natom,3,"atpos");

  // read comment line
  fgets(str, MAXLINE, fp);
  int n = strlen(str)+1;
  comment = memory->create(comment,n,"comment");
  strcpy(comment, str);

  // read atomic info
  int flag[3]; flag[0] = flag[1] = flag[2] = 0;
  for (int i=0; i<natom; i++){
    fgets(str, MAXLINE, fp);
    n = count_words(str);
    char *ptr = strtok(str," \n\t\r\f");
    attyp[i] = name2type(ptr);
    for (int j=0; j<3; j++) atpos[i][j] = atof(strtok(NULL," \n\t\r\f"));
    if (n >= 9){
      strtok(NULL," \n\t\r\f");
      int idim = atoi(strtok(NULL," \n\t\r\f"))-1;
      if (idim >= 0 &&  idim < 3){
        for (int j=0; j<3; j++) latvec[idim][j] = atof(strtok(NULL," \n\t\r\f"));
        flag[idim] = 1;
      }
    }
  }

  // check lattice vectors info
  if (flag[0]+flag[1]+flag[2] < 3){
    printf("\nAtomic configuration read from file: %s\n", fname);
    printf("Lattice vector info is however insufficient, please input them now:\n");
    for (int i=0; i<3; i++){
      while (1){
        printf("Please input vector A%d: ", i+1);
        n = count_words(fgets(str,MAXLINE,stdin));
        if (n < 3) continue;
        latvec[i][0] = atof(strtok(str, " \n\t\r\f"));
        latvec[i][1] = atof(strtok(NULL," \n\t\r\f"));
        latvec[i][2] = atof(strtok(NULL," \n\t\r\f"));

        break;
      }
    }
  }

  cartesian = 1;
  initialized = 1;

  // display lattice info
  //display(file);

return;
}

/* ----------------------------------------------------------------------
   deconstructor destroy any allocated memory
------------------------------------------------------------------------- */
LattXYZ::~LattXYZ()
{
  delete []fname;
  names.clear();
  if (atpos) memory->destroy(atpos);
  if (attyp) memory->destroy(attyp);
  if (comment) memory->destroy(comment);
  
  delete memory;
}

/* ----------------------------------------------------------------------
   Display LattXYZ basic info
------------------------------------------------------------------------- */
void LattXYZ::display(const char * file)
{
  if (!initialized) return;

  printf("\n"); for (int i=0; i<28; i++) printf("=");
  printf(" Lattice Info "); for (int i=0; i<28; i++) printf("="); printf("\n");
  printf("Lattice info read from file.......: %s\n", file);
  printf("Number of atom in system..........: %d\n", natom);
  printf("Number of atom types in system....: %d\n", ntype);
  printf("Elements found in your system.....: ");
  for (int i=0; i<ntype; i++) printf(" %s", names[i].c_str());
  printf("\nLattice vectors:\n");
  for (int i=0; i<3; i++){
    for(int j=0; j<3; j++) printf("%lf ", latvec[i][j]*alat);
    printf("\n");
  }
  for (int i=0; i<70; i++) printf("-");
  printf("\nBasis (Fractional coordinate & type):\n");
  for (int i=0; i<natom; i++){
    printf("%lf %lf %lf %d\n", atpos[i][0], atpos[i][1], atpos[i][2], attyp[i]);
  }
  for (int i=0; i<70; i++) printf("="); printf("\n\n");
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int LattXYZ::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}

/*------------------------------------------------------------------------------
 * Private method to get the cross product of two vectors; C[3] = A[3] X B[3]
 *----------------------------------------------------------------------------*/
void LattXYZ::Cross(double *A, double *B, double *C)
{
   C[0] = A[1] * B[2] - A[2] * B[1];
   C[1] = A[2] * B[0] - A[0] * B[2];
   C[2] = A[0] * B[1] - A[1] * B[0];

return;
}

/*------------------------------------------------------------------------------
 * Private method to get the dot product of two vectors of dim 3
 *----------------------------------------------------------------------------*/
double LattXYZ::DotProd(double *A, double *B)
{
return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

/*------------------------------------------------------------------------------
 * Private method to get the atomic type number
 *----------------------------------------------------------------------------*/
int LattXYZ::name2type(const char *ename)
{
  std::string namestr; namestr.assign(ename);

  int ip = -1;
  for (int i=0; i<ntype; i++){
    if (strcmp(ename, names[i].c_str()) == 0){ ip = i; break; }
  }
  if (ip < 0){
    ip = ntype++;
    names.push_back(namestr);
  }
  namestr.clear();

return ip+1;
}

/*------------------------------------------------------------------------------
 * Private method to do matrix inversion
 *----------------------------------------------------------------------------*/
void LattXYZ::GaussJordan(const int n, const double *MatA, double *Mat)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
  int indxc[n],indxr[n],ipiv[n];
  double big, dum, pivinv;

  for (int i=0; i<n*n; i++) Mat[i] = MatA[i];

  for (i=0; i<n; i++) ipiv[i] = 0;
  for (i=0; i<n; i++){
    big = 0.;
    for (j=0; j<n; j++){
      if (ipiv[j] != 1){
        for (k=0; k<n; k++){
          if (ipiv[k] == 0){
            idr = j*n+k;
            if (fabs(Mat[idr]) >= big){
              big  = fabs(Mat[idr]);
              irow = j;
              icol = k;
            }
          }else if (ipiv[k] >1){
            printf("\nError: Singular matrix in double GaussJordan!\n");
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l=0; l<n; l++){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == 0.) printf("\nError: Singular matrix in double GaussJordan!\n");
    
    pivinv = 1./ Mat[idr];
    Mat[idr] = 1.;
    idr = icol*n;
    for (l=0; l<n; l++) Mat[idr+l] *= pivinv;
    for (ll=0; ll<n; ll++){
      if (ll != icol){
        idc = ll*n+icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l=0; l<n; l++) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }
  for (l=n-1; l>=0; l--){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k=0; k<n; k++){
        idr = k*n+rl;
        idc = k*n+cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }

return;
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinate into fractional
 *----------------------------------------------------------------------------*/
void LattXYZ::car2dir()
{
  double invaxis[3][3];
  GaussJordan(3,&latvec[0][0],&invaxis[0][0]);
  double x[natom][3], **s = atpos;
  for (int i=0; i<natom; i++)
  for (int idim=0; idim<3; idim++) x[i][idim] = atpos[i][idim];

  for (int i=0; i<natom; i++){
    for (int idim=0; idim<3; idim++){
      s[i][idim] = 0.;
      for (int jdim=0; jdim<3; jdim++) s[i][idim] += x[i][jdim]*invaxis[jdim][idim];

      while (s[i][idim] >= 1.) s[i][idim] -= 1.;
      while (s[i][idim] <  0.) s[i][idim] += 1.;
    }
  }

  cartesian = 0;
return;
}

/*------------------------------------------------------------------------------
 * Method to convert fractional coordinate into cartesian
 *----------------------------------------------------------------------------*/
void LattXYZ::dir2car()
{
  double s[natom][3], **x = atpos;
  for (int i=0; i<natom; i++)
  for (int idim=0; idim<3; idim++) s[i][idim] = atpos[i][idim];

  for (int i=0; i<natom; i++){
    for (int idim=0; idim<3; idim++){
      x[i][idim] = 0.;
      for (int jdim=0; jdim<3; jdim++) x[i][idim] += s[i][jdim]*latvec[jdim][idim];
    }
  }
  cartesian = 1;
return;
}

/*------------------------------------------------------------------------------
 * Method to convert fractional coordinate into cartesian
 * Note: axis[][i] carries the i'th base vector
 *----------------------------------------------------------------------------*/
void LattXYZ::dir2car(const double axis[3][3],int ntm, double s[][3])
{
  double x[ntm][3];
  for (int i=0; i<ntm; i++){
    for (int idim=0; idim<3; idim++){
      x[i][idim] = 0.;
      for (int jdim=0; jdim<3; jdim++) x[i][idim] += s[i][jdim]*axis[idim][jdim];
    }
  }
  for (int i=0; i<ntm; i++)
  for (int idim=0; idim<3; idim++) s[i][idim] = x[i][idim];

return;
}
/*----------------------------------------------------------------------------*/
