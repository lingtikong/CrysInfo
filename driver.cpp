#include "driver.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#define PWIDTH  70
#define MAXLINE 256
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))
/* -----------------------------------------------------------------------------
 * driver, main menu
 * -------------------------------------------------------------------------- */
Driver::Driver(int narg, char **arg)
{
  // read all input files, if supplied from command line
  nlat = nmax = 0;
  all.clear();
  for (int iarg=1; iarg<narg; iarg++){
    if (strcmp(arg[iarg], "-h") == 0){
      help();

    } else {
      readxyz(arg[iarg]);
    }
  }

  char str[MAXLINE];
  while (1){
    int job = 0;
    // Main menu
    printf("\n"); for (int i=0; i<PWIDTH; i++) printf("=");
    printf("\nCurrent number of configurations read: %d\n", nlat);
    for (int i=0; i<PWIDTH; i++) printf("-");
    printf("\nPlease select the job you would like to do:\n");
    printf("    1. Compute the RMSD between xyz files;\n");
    printf("    2. Symmetry related info for xyz files;\n");
    printf("    0. exit.\nYour choice [0]: ");
    fgets(str,MAXLINE,stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) job = atoi(ptr);
    printf("Your selection: %d\n", job);
    for (int i=0; i<PWIDTH; i++) printf("="); printf("\n");
  
    // main drivers
    if      (job == 1) compare();
    else if (job == 2) symmetry();
    else break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Deconstructor, free memories
 * -------------------------------------------------------------------------- */
Driver::~Driver()
{
  for (int i=0; i<nlat; i++){
    one = all[i];
    delete one;
  }
  all.clear();
  one = NULL;

return;
}

void Driver::check_config(const int nmin)
{
  char str[MAXLINE];
  while (nlat < nmin){
    printf("The selected job requires at least %d configuration.\n", nmin);
    printf("Please input the file that has the No.%d config: ", nlat+1);
    fgets(str,MAXLINE,stdin);
    char *ptr = strtok(str," \n\t\r\f");
    if (ptr == NULL) continue;

    readxyz(ptr);
  }
return;
}

/* -----------------------------------------------------------------------------
 * To compare the RMSD of xyz configurations with respect to the first one
 * -------------------------------------------------------------------------- */
void Driver::compare()
{
  char str[MAXLINE];
  RMSD *rmsd = new RMSD();

  // ask if not enough configurations
  check_config(2);

  // the first configuraiton is set as reference
  one = all[0];
  if (one->cartesian != 1) one->dir2car();
  int nref = one->natom;
  double pos_ref[nref][3], pos_loc[nref][3];
  double mcom[3], v2ref[3], U[3][3];
  for (int i=0; i<nref; i++)
  for (int j=0; j<3; j++) pos_ref[i][j] = one->atpos[i][j];

  // compute rmsd one by one
  printf("\n"); for (int i=0; i<PWIDTH; i++) printf("=");
  printf("\nRMSD with respect to configurations read from: %s,\n  in unit of A/atom.\n", one->fname);
  printf("Configurations with different # of atom to %s is skipped.\n", one->fname);
  for (int i=0; i<PWIDTH; i++) printf("-");  printf("\n");
  for (int ii=1; ii<nlat; ii++){
    // assign configuration, skip if # of atom mismatch
    one = all[ii];
    if (one->natom != nref) continue;
    if (one->cartesian != 1) one->dir2car();

    for (int i=0; i<nref; i++)
    for (int j=0; j<3; j++) pos_loc[i][j] = one->atpos[i][j];
    
    // to compute the rmsd
    double res;
    rmsd->calculate_rotation_rmsd(pos_ref, pos_loc, nref, mcom, v2ref, U, &res);

    printf("%20.16f   | %s\n", res, one->fname);
  }
  for (int i=0; i<PWIDTH; i++) printf("="); printf("\n\n");

  delete rmsd;

return;
}

/* -----------------------------------------------------------------------------
 * To write out the current configuration in xyz format
 * -------------------------------------------------------------------------- */
void Driver::writexyz(FILE *fp, int ntm, double axis[3][3], double pos[][3], int type[], LattXYZ * latt)
{
  fprintf(fp,"Basis vectors of the lattice:\n");
  for (int i=0; i<3; i++) fprintf(fp,"  %lg %lg %lg  (A%d)\n", axis[0][i], axis[1][i], axis[2][i], i+1);

  for (int i=0; i<ntm; i++)
  for (int j=0; j<3; j++) pos[ntm+i][j] = pos[i][j];
  latt->dir2car(axis, ntm, &pos[ntm]); // fractional to cartesian

  fprintf(fp,"Atomic configuration in xyz format:\n");
  fprintf(fp,"  %d\n  Atomic configuration from: %s\n", ntm, latt->fname);
  for (int i = 0; i < MIN(3, ntm); ++i) fprintf(fp,"  %2s %20.14g %20.14g %20.14g crystal_vector %d %lg %lg %lg\n", latt->names[type[i]-1].c_str(),
    pos[ntm+i][0], pos[ntm+i][1], pos[ntm+i][2], i+1, axis[i][0], axis[i][1], axis[i][2]);

  for (int i = MIN(3, ntm); i < ntm; ++i) fprintf(fp,"  %2s %20.14g %20.14g %20.14g\n", latt->names[type[i]-1].c_str(),
    pos[ntm+i][0], pos[ntm+i][1], pos[ntm+i][2]);

return;
}

void Driver::readxyz(const char *file)
{
  FILE *fp = fopen(file, "r");
  if (fp == NULL){
    printf("\nError: file %s not found!\n", file);
    return;
  }

  while (! feof(fp)){
    one = new LattXYZ(fp, file);

    //printf("\nFrom readxyz: initialized = %d\n", one->initialized);
    if (one->initialized == 1){
      all.push_back(one);
      nmax = MAX(nmax, one->natom);
    } else {
      delete one;

      break;
    }
  }
  fclose(fp);

  nlat = all.size();
return;
}
/* -----------------------------------------------------------------------------
 * To evaluate symmetry related info for xyz configurations
 * -------------------------------------------------------------------------- */
void Driver::symmetry()
{
  // ask if not enough configuration read
  check_config(1);

  char str[MAXLINE];
  // precision in atomic coordinate when evaluating symmetry
  double symprec = 1.e-5;
  printf("\nPlease define your desired precison in the atomic coordinates\n");
  printf("when computing symmetry [%g]: ", symprec);
  fgets(str, MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) symprec = atof(ptr);
  
  char choice[MAXLINE];
  strcpy(choice,"1a");

  while (1){
    // menu
    printf("\n"); for (int i=0; i<PWIDTH; i++) printf("=");
    printf("\nPlease select your job:\n");
    for (int i=0; i<PWIDTH; i++) printf("-"); printf("\n");
    printf("  1. Show                      |    |\n");
    printf("  2. Find symmetry operations  |    |  a. original cell;\n");
    printf("  3. Find group symbol         |    |  b. refined cell;\n");
    printf("  4. Find Schoenflies symbol   | of |  c. primitive cell;\n");
    printf("  5. Find irreducible k-points |    |  d. primitive of refined cell;\n");
    printf("  6. Assign Wyckoff letter     |    |\n");
    for (int i=0; i<PWIDTH; i++) printf("-");
    printf("\n  7. Redefine precision;              0. Return;\n");
    for (int i=0; i<PWIDTH; i++) printf("-");
    printf("\nYour choice [%s]: ", choice);
    fgets(str, MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) strcpy(choice, ptr);
    printf("Your selection: %s\n", choice);
    for (int i=0; i<PWIDTH; i++) printf("="); printf("\n");

    if (strlen(choice) < 1) break;
    int job  = choice[0]-'1';

    if (job < 0 || job > 6) break;
    else if (job == 6){
      printf("\nPlease define your desired precison in the atomic coordinates\n");
      printf("when computing symmetry [%g]: ", symprec);
      fgets(str, MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) symprec = atof(ptr);

      continue;
    }

    int cell = choice[1]-'a';

    // the result could be written to file or stdout
    FILE *fp = stdout;
    char *fout; fout = NULL;
    int flag = 0;
    if (job > 0){
      printf("Write results to file? (y/n)[n]: ");
      fgets(str, MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr && (strcmp(ptr,"y")==0 || strcmp(ptr,"Y")==0 || strcmp(ptr,"yes")!=0)){
        printf("Please input the output file name [result.dat]: ");
        fgets(str, MAXLINE, stdin);
        ptr = strtok(str, " \n\t\r\f");
        if (ptr){
          fout = new char[strlen(ptr)+1];
          strcpy(fout, ptr);
        } else {
          fout = new char [11];
          strcpy(fout, "result.dat");
        }
        fp = fopen(fout,"w");
        flag = 1;
      }
    }

    // local memory space, make it large enough for all
    int natom;
    int attyp[nmax*10];
    double latvec[3][3], atpos[nmax*10][3];

    int flag_mesh = 0, ngrid;
    int mesh[3], shift[3], time_reversal;
    shift[0] = shift[1] = shift[2] = time_reversal = 0;

    fprintf(fp,"\n");for (int i=0; i<PWIDTH; i++) fprintf(fp,"="); fprintf(fp,"\n");
    // loop over all read configurations
    for (int ii=0; ii<nlat; ii++){
      one = all[ii];
      natom = one->natom;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++) latvec[i][j] = one->latvec[j][i];
    
      if (one->cartesian) one->car2dir(); // fractional required by spg
      for (int i=0; i<natom; i++){
        attyp[i] = one->attyp[i];
        for (int j=0; j<3; j++) atpos[i][j] = one->atpos[i][j];
      }

      // get refined if required
      if (cell == 1 || cell == 3){
        int nref = spg_refine_cell(latvec,atpos,attyp,natom,symprec);
        if (nref > 0) natom = nref;
        //else assignOne(one, &natom, latvec, atpos, attyp);
        fprintf(fp,"Refined cell,   originally from file: %s (%d atoms)\n", one->fname, one->natom);
        writexyz(fp, natom, latvec, atpos, attyp, one);
      }
  
      // get primitive if required
      if (cell == 2 || cell == 3){
        int npriv = spg_find_primitive(latvec,atpos,attyp,natom,symprec);
        if (npriv > 0) natom = npriv;
        fprintf(fp,"Primitive cell, originally from file: %s (%d atoms)\n", one->fname, one->natom);
        writexyz(fp, natom, latvec,atpos,attyp,one);
      }

      if (job == 0){ // show the current configurations
        if (cell == 0) writexyz(fp, natom, latvec,atpos,attyp,one);

      } else if (job == 1){ // all symmetry operations

        // max # of possible operations
        int nmax_mul = spg_get_max_multiplicity(latvec,atpos,attyp,natom,symprec);

        // to do the job
        if (nmax_mul > 0){
          int rotation[nmax_mul][3][3];
          double translation[nmax_mul][3];
          int nop = spg_get_symmetry(rotation,translation,nmax_mul,latvec,atpos,attyp,natom,symprec);

          // write the result
          fprintf(fp,"\nSymmetry operations of the current configuration:\n");
          for (int i=0; i<nop; i++){
            fprintf(fp,"Operation %d of %d:\n", i+1, nop);
            fprintf(fp,"  Translation: %lg %lg %lg\n", translation[i][0],translation[i][1],translation[i][2]);
            fprintf(fp,"  Rotation: (");
            for (int j=0; j<3; j++) fprintf(fp," %d %d %d ;", rotation[i][j][0],rotation[i][j][1],rotation[i][j][2]);
            fprintf(fp,")\n");
          }
        }

      } else if (job == 2){ // international symbol
        char symbol[11];

        // get the International group symbol
        int id = spg_get_international(symbol,latvec,atpos,attyp,natom,symprec);

        // output the result
        fprintf(fp,"\nSpace group symbol of current config (%s): ", one->fname);
        fprintf(fp,"%s <--> %d\n", strtok(symbol," \n\t\r\f"),id);

      } else if (job == 3){ // Schoenflies symbol
        char symbol[10];

        // get the result
        int id = spg_get_schoenflies(symbol,latvec,atpos,attyp,natom,symprec);

        // output the result
        fprintf(fp,"\nSchoenflies symbol of current config (%s): ", one->fname);
        fprintf(fp,"%s <--> %d\n", strtok(symbol," \n\t\r\f"),id);

      } else if (job == 4){ //  Find irreducible k-points
        int mesh[3], shift[3], time_reversal;
        shift[0] = shift[1] = shift[2] = time_reversal = 0;

        // ask for k-mesh size
        if (flag_mesh == 0){
          while (1){
            printf("\nPlease input your desired k-mesh size: ");
            fgets(str,MAXLINE,stdin);
            ptr = strtok(str, " \n\t\r\f");
            if (ptr == NULL) continue;
            mesh[0] = atoi(ptr); if (mesh[0] < 1) continue;
    
            ptr = strtok(NULL," \n\t\r\f");
            if (ptr == NULL) continue;
            mesh[1] = atoi(ptr); if (mesh[1] < 1) continue;
    
            ptr = strtok(NULL," \n\t\r\f");
            if (ptr == NULL) continue;
            mesh[2] = atoi(ptr); if (mesh[2] < 1) continue;
    
            break;
          }
          ngrid = mesh[0]*mesh[1]*mesh[2];

          flag_mesh = 1;
        }

        // get the irreducible grids
        int grid_point[ngrid][3], map[ngrid];
        int nir = spg_get_ir_reciprocal_mesh(grid_point, map, mesh, shift, time_reversal, latvec, atpos, attyp, natom, symprec);

        // find the irreducible q-points and their weights
        int *id2ir;
        double **q, *w;
        q = one->memory->create(q,nir,3,"q");
        w = one->memory->create(w,nir,"w");
        id2ir = one->memory->create(id2ir,ngrid,"id2ir");
        for (int i=0; i<nir; i++) w[i] = 0.;

        nir = 0;
        for (int i=0; i<ngrid; i++){
          if (map[i] == i) id2ir[i] = nir++;
        }
        nir = 0;
        for (int i=0; i<ngrid; i++){
          w[id2ir[map[i]]] += 1.;
          if (map[i] == i){
            for (int j=0; j<3; j++) q[nir][j] = double(grid_point[i][j])/double(mesh[j]);
            nir++;
          }
        }
        for (int i=0; i<nir; i++) w[i] /= double(ngrid);

        // output the q-points
        fprintf(fp,"\nBasis vectors of current cell (originally from %s):\n", one->fname);
        for (int i=0; i<3; i++) fprintf(fp,"  %lg %lg %lg\n", latvec[0][i],latvec[1][i],latvec[2][i]);
        fprintf(fp,"Irreducible k-points from mesh %d x %d x %d => %d\n", mesh[0], mesh[1], mesh[2], nir);
        fprintf(fp,"  # qx qy qz (reduced unit) weight\n");
        for (int i=0; i<nir; i++) fprintf(fp,"  %lg %lg %lg %lg\n", q[i][0], q[i][1], q[i][2], w[i]);

        // free memory
        one->memory->destroy(q);
        one->memory->destroy(w);
        one->memory->destroy(id2ir);

      } else if (job == 5){ // Wyckoff letter assignment
        SpglibDataset *dataset;

        // get the result
        dataset = spg_get_dataset(latvec,atpos,attyp,natom,symprec);

        for (int i=0; i<natom; i++)
        for (int j=0; j<3; j++) atpos[i+natom][j] = atpos[i][j];
        one->dir2car(latvec,natom,&atpos[natom]);

        // count wyckoff
        int wyck[natom], wmin = natom, wmax = -natom;
        for (int i=0; i<natom; i++) wyck[i] = 0;
        for (int i=0; i<natom; i++){
          int id = dataset->wyckoffs[i];
          wmin = MIN(wmin,id);
          wmax = MAX(wmax,id);
          wyck[id]++;
        }

        // output the result
        fprintf(fp,"\nSymmetry info of the current cell:\n");
        fprintf(fp,"  Original config file : %s\n", one->fname);
        fprintf(fp,"  Space  group  number : %d\n", dataset->spacegroup_number);
        fprintf(fp,"  International symbol : %s\n", strtok(dataset->international_symbol," \n\r\t\f"));
        fprintf(fp,"  Hall   symbol        : %s\n", strtok(dataset->hall_symbol," \n\r\t\f"));
        fprintf(fp,"  # of symmetry operate: %d\n", dataset->n_operations);
        fprintf(fp,"  # of Wyckoff letters : %d\n", (wmax-wmin)+1);
        fprintf(fp,"  Distinct Wyckoff lett:");
        for (int i=wmin; i<=wmax; i++)  fprintf(fp," %d%c", wyck[i], i+'a');
        fprintf(fp,"\nAtomic configuration (xyz) with Wyckoff info appended:\n");
        fprintf(fp,"  %d\n  # atom x y z wyckoff equiv-atom rx ry rz\n", natom);
        for (int i=0; i<natom; i++) fprintf(fp,"  %4s %12.6f %12.6f %12.6f %c %d -> %d %g %g %g\n", one->names[attyp[i]-1].c_str(),
          atpos[natom+i][0], atpos[natom+i][1], atpos[natom+i][2], dataset->wyckoffs[i]+'a', i+1,dataset->equivalent_atoms[i]+1,
          atpos[i][0], atpos[i][1], atpos[i][2]);

        spg_free_dataset(dataset);
      }
      if (ii<nlat-1){ for (int i=0; i<PWIDTH; i++) fprintf(fp,"-"); fprintf(fp,"\n");}
    }
    for (int i=0; i<PWIDTH; i++) fprintf(fp,"="); fprintf(fp,"\n");

    if (flag){
      fclose(fp);
      printf("\nJob done, view the result? (y/n)[n]: ");
      fgets(str, MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr){
        sprintf(str,"more %s\n", fout);
        system(str);
        delete []fout; fout = NULL;
      }
    }
  }

return;
}

void Driver::help()
{
  exit(0);
}
