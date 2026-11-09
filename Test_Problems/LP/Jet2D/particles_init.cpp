////////////////////////////////////////////////////////////////////////
/*!
 \file
 \brief Initialize particle distribution, set B.C. and injection.

 This file contains routines to initialize particles on the grid,
 assign user-defined boundary conditions and inject particles.
 Particles attributes that can be set are: position, velocity, color.
 In case of evolution with spectra, the initial spectral profile is also
 prescribed here.

 \authors A. Suriano (alessio.suriano@unito.it)\n
          V. Berta (vittoria.berta@unito.it)\n
          A. Mignone (andrea.mignone@unito.it)
         
  \date   March 18, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include <fstream>

#if PARTICLES_INIT_OLD   == NO

long int GetFileSize(FILE *);
char *lastline(char *);
void truncate_file(const char *, int);

////////////////////////////////////////////////////////////////////////
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
*/
////////////////////////////////////////////////////////////////////////
void Particles_Init(Data *data, Grid *grid)
{

  int i, j, k, nv;
  int npCell  = RuntimeGet()->Nparticles_per_block;
  int ncBlock = RuntimeGet()->Ncells_per_block;
  int toBeInsertedGlob = 0;
  int toBeInserted     = 0;
  static int first_call = 1;

  static const double invNpCell = 1.0/static_cast<double>(npCell);
  ParticlesContainer *pc = data->pc;
  RBox box, boxg;

  // 0. Coarse grid for under-sampled distribution 
  Grid pc_grid;
  Runtime &rtime = *(RuntimeGet());  // Copy runtime structure
  SetGrid (&rtime, grid->nproc, &pc_grid, ncBlock);
  #pragma acc enter data copyin(pc_grid)  // Should be present in both cases

  if (first_call == 1){
    pc->initUnits(*(data->run_config));
    pc->insert.create (pc_grid.ntot[IDIR], pc_grid.ntot[JDIR], pc_grid.ntot[KDIR]);
    #pragma acc enter data copyin (pc->insert)
  } 

  // ------------------------------------------------------
  // 1. Initializes particles on a cell by cell basis.
  // ------------------------------------------------------

  int prvSize = pc->Size();

  // 1b. Define the boxes

  RBoxDefine(pc_grid.lbeg, pc_grid.lend, &box);
  RBoxDefine(pc_grid.gbeg, pc_grid.gend, &boxg);
  #pragma acc enter data copyin(box, boxg)

  // ----------------------------------------------------
  // 1a. On first call, we assign initial condition by
  //     setting the pc->insert() array.
  // ----------------------------------------------------

  if (first_call){   // at t = 0
    #pragma acc parallel loop collapse(3) present(data, pc, pc_grid) \
                copy (toBeInserted) reduction (+:toBeInserted)
    BOX_LOOP(box, i, j, k){
      bool doInsert = 0;  

      if(doInsert){
        pc->insert(i,j,k) = true;
        toBeInserted++;
      } else {
        pc->insert(i,j,k) = false;
      }
    }  // BOX_LOOP(box, i, j, k) 
  }

  if (!first_call){  // at t > 0, used for injection
    #pragma acc parallel loop collapse(3) present(box,pc) copy (toBeInserted) \
                reduction (+:toBeInserted)
    BOX_LOOP(box, i, j, k){
      if (pc->insert(i,j,k)) toBeInserted++;
    } 
  }

  toBeInserted *= npCell;

  // ----------------------------------------------------
  // 1b. Count particles from other processors
  // ----------------------------------------------------

  #ifdef PARALLEL
  MPI_Allreduce(&toBeInserted, &toBeInsertedGlob, 1, 
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #else
  toBeInsertedGlob = toBeInserted;
  #endif

  // ----------------------------------------------------
  // 2. Set paticles properties (position, etc..)
  // ----------------------------------------------------

  if (toBeInserted > 0){

    // 2a. Allocate space in the container saving the prevoius size

    pc->Reserve(toBeInserted);
    int index = prvSize;

    // 2b. Insert npCell particles per cell based on the value of
    //     the array insert

    const double time = data->Dts->time;

    #pragma acc parallel loop gang collapse(3) present(data, boxg, pc_grid) \
                         copyin(index)
    BOX_LOOP(box, i, j, k){
      if (pc->insert(i,j,k)){
        double xL[3], xR[3];
        xL[IDIR] = pc_grid.xL[IDIR][i]; xR[IDIR] = pc_grid.xR[IDIR][i];
        xL[JDIR] = pc_grid.xL[JDIR][j]; xR[JDIR] = pc_grid.xR[JDIR][j];
        xL[KDIR] = pc_grid.xL[KDIR][k]; xR[KDIR] = pc_grid.xR[KDIR][k];

        // 2c. Associate a unique index to (i,j,k)

        int ig = pc_grid.beg[IDIR] + i - pc_grid.lbeg[IDIR];
        int jg = pc_grid.beg[JDIR] + j - pc_grid.lbeg[JDIR];
        int kg = pc_grid.beg[KDIR] + k - pc_grid.lbeg[KDIR];
        uint64_t progressive_cell = BOX_INDEX(boxg, ig, jg, kg)*npCell;

        // 2d. Loop over particles in the given cell

        #pragma acc loop seq
        for (int np = 0; np < npCell; np++){  

          int indexLoc;   
          #pragma acc atomic capture
          indexLoc = index++;

          // 2d1. Dispose particles in the cell

          Particles_LoadUniform(np, invNpCell, xL, xR, pc->pos(indexLoc));

          // 2d2. Set particle attributes 

          #if PARTICLES_LP_SPECTRA == YES
          pc->flag(indexLoc)    = 0;
          pc->flagPrv(indexLoc) = 0;
          #if PARTICLES_LP_SHOCK_SPECTRA == YES
          pc->shkGradp(indexLoc) = -1.0;
          NVAR_LOOP(nv) pc->VshkUp(indexLoc)[nv] = -1.0;
          NVAR_LOOP(nv) pc->VshkDn(indexLoc)[nv] = -1.0;
          #endif
          pc->cr(indexLoc) = 0.0;
          pc->lorG(indexLoc) = 1.0;

          Particles_LP_InitSpectra( pc->eng(indexLoc), 
                                    pc->chi(indexLoc), 
                                    pc->nmicro(indexLoc), npCell);
          #endif

          pc->id(indexLoc)   = progressive_cell + np + pc->idCounter;
          pc->tinj(indexLoc) = time;
          #pragma acc loop vector
          for (int nc = 0; nc < PARTICLES_NCOLORS; nc++) {
            double xp = pc->pos(indexLoc)[IDIR];
            double yp = pc->pos(indexLoc)[JDIR];

            pc->color(indexLoc)[nc] = xp/grid->dx[IDIR][0] + 100*(yp+0.00001)/grid->dx[JDIR][0];
          }
        } // for (int np = 0; np < npCell_; np++) 
      } // if(pc->insert(i,j,k)
    } // BOX_LOOP(box, i, j, k) 
  } // if (toBeInserted > 0) 

  // End Cell by cell initialization

  int dir;
  
  // 3. Remove particles accidentally assigned outside the physical domain

  #pragma acc loop seq
  DIM_LOOP(dir) {
    double xbeg = grid->xbeg[dir];
    double xend = grid->xend[dir];
    (pc)->BoundaryRemove(dir, xbeg, -1);
    (pc)->BoundaryRemove(dir, xend, 1);
  }

  Particles_Boundary(data, grid);

  // 4. Store fluid quantities at the newly allocated particle position

  Particles_LP_FixValue(data, grid, true, prvSize);

  #if PARTICLES_LP_SPECTRA == YES
  const uint64_t nptot = pc->Size();
  #pragma acc parallel loop present(data)
  for(uint64_t np = prvSize; np < nptot; np ++){
    #pragma acc loop
    for (int e = 0; e < PARTICLES_LP_NEBINS; e++){
       pc->chi(np)[e] /= pc->density(np); // Chi is a fraction and not number density 
    }
  }
  #endif // PARTICLES_LP_SPECTRA 

  pc->idCounter += toBeInsertedGlob;
  
  #pragma acc update device(pc->idCounter)
  #pragma acc exit data delete(box, boxg, pc_grid)

  first_call = 0;
}


////////////////////////////////////////////////////////////////////////
/*!
   Sets user-defined boundary conditions on particles.
 
   \param [in]  data    Pointer to the PLUTO data structure.
   \param [in]  grid    Pointer to the PLUTO grid structure.
*/
////////////////////////////////////////////////////////////////////////
void Particles_Inject(Data *data, Grid *grid)
{
  int i,j,k;
  ParticlesContainer *pc = data->pc;
  RunConfig &run_config = *(data->run_config);
  RBox box;

  RBoxDefine(0,grid->ntot[IDIR]-1,
             0,grid->ntot[JDIR]-1,
             0,grid->ntot[KDIR]-1, &box);

  double dtinj = run_config.inputParam[INJECTION_DT];
  double t     = data->Dts->time;
  double tnext = data->Dts->time + data->Dts->dt;
  int time_to_inject = int(t/dtinj) != int(tnext/dtinj);
  if (time_to_inject){
    BOX_LOOP(box, i, j, k){
      double x1 = grid->xC[IDIR][i];
      double x2 = grid->xC[JDIR][j];
      double rjc = 0.01;

      pc->insert(i,j,k) = false;
      if ( (sqrt(x2*x2) < rjc) && (x1 >= 0.006) && (x1 <= 0.01)  ){
        pc->insert(i,j,k) = true;
      }
    }
    Particles_Init(data, grid);
  }
}


#if PARTICLES_LP_SPECTRA == YES
////////////////////////////////////////////////////////////////////////
/*!
   Sets user-defined boundary conditions on particles.
 
   \param [in]  data    Pointer to the PLUTO data structure.
   \param [in]  grid    Pointer to the PLUTO grid structure.
*/
////////////////////////////////////////////////////////////////////////
void Particles_LP_InitSpectra(NebinsE &eng, Nebins &chi, double &nmicro, 
                              const int npCell)
{
  constexpr double Emin = 1.0e2;
  constexpr double Emax = 1.0e4;
  const double lnEmin = log10(Emin);
  const double lnEmax = log10(Emax);
  const double dlnE = (lnEmax - lnEmin)/(static_cast<double>(PARTICLES_LP_NEBINS));

  const double alpha1 = 9.0;
  nmicro = 0.000001;

  // ------------------------------------------------------
  // 1. Initialize spectral distribution to power law.
  //    Chi in each bin is assigned as an average over the 
  //    energy interval.
  //    chi_i = \int^Eh_El E^(-p) dE /(Eh - El)
  // ------------------------------------------------------

  const double scrhchi = (nmicro)/(pow(Emax,1.0-alpha1)-pow(Emin,1.0-alpha1));

  eng[0] = Emin;

  // there's dependency on eng[]
  #pragma acc loop seq
  for (int i = 0; i < PARTICLES_LP_NEBINS; i++){

    double scrhln = lnEmin + static_cast<double>(i + 1)*dlnE;

    eng[i+1] = pow(10.0, scrhln);
    chi[i]   = scrhchi*( pow(eng[i+1],1.0-alpha1)
                        -pow(eng[i  ],1.0-alpha1));
    chi[i]  /= eng[i+1] - eng[i];
  }

}

#endif

////////////////////////////////////////////////////////////////////////
/*
 *
*/
////////////////////////////////////////////////////////////////////////
void Particles_UserDefBoundary(Data *data, const int side, const Grid *grid)
{
  int    dir;
  double xbeg[3], xend[3];

  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  if (side == X1_BEG){
    dir = IDIR;
      (data->pc)->BoundaryRemove(dir,xbeg[dir], -1);
  }

  if (side == X1_END){
    dir = IDIR;
      (data->pc)->BoundaryRemove(dir,xend[dir], 1);
  }

  if (side == X2_BEG){
    dir = JDIR;
      (data->pc)->BoundaryRemove(dir,xbeg[dir], -1);
  }

  if (side == X2_END){
    dir = JDIR;
      (data->pc)->BoundaryRemove(dir,xend[dir], 1);
  }

  if (side == X3_BEG){
    dir = KDIR;
      (data->pc)->BoundaryRemove(dir,xbeg[dir], -1);
  }

  if (side == X3_END){
    dir = KDIR;
      (data->pc)->BoundaryRemove(dir,xend[dir], 1);
  }
}


////////////////////////////////////////////////////////////////////////
/*
 *
*/
////////////////////////////////////////////////////////////////////////
void ParticlesAnalysis(Data *data, Grid * g){
#ifdef PARTICLES_LP_STATS_NFIELDS
  char fname[256], metaFname[256];
  sprintf(fname, "%s/update.%d.txt", RuntimeGet()->output_dir, g_rank);
  sprintf(metaFname, "%s/metaUpdate.%d.txt", RuntimeGet()->output_dir, g_rank);

  snprintf(fname, sizeof(fname),"%s/update.%d.txt", 
           RuntimeGet()->output_dir, g_rank);

  snprintf(metaFname, sizeof(metaFname),"%s/metaUpdate.%d.txt",
           RuntimeGet()->output_dir, g_rank);

  printf("> Writing userdef output %s\n", fname);
  static long int nCalls = -1;

  if(g_stepNumber != 0 && nCalls == -1){
    int nDeleteU = 0;
    int nDeleteM = 0;
    //reding nCalls
    FILE *fp;
    fp = fopen(metaFname, "r");
    char line[64];

    while (fgets(line, sizeof(line), fp)) {
      int nL, stepL;
      int n = sscanf(line, "%ld;%d;%d", &nCalls, &stepL, &nL);
      if(stepL >= g_stepNumber) {
        nDeleteU += nL;
        nDeleteM++;
      }
    }


    fclose(fp);

    //delete last nLines
    truncate_file(fname, nDeleteU);
    truncate_file(metaFname, nDeleteM);
    nCalls -= nDeleteM;
  }
  nCalls++;

  uint64_t nptot = data->pc->Size();
  int i = 0;
  int j = 0;
  double *dump_arr;
  constexpr int N = PARTICLES_NCOLORS - PARTICLES_LP_STATS_OFFSET;
  constexpr int PARTICLES_LP_STATS_NFIELDS_LOC = PARTICLES_LP_STATS_NFIELDS + 1;
  const int SIZE = (PARTICLES_LP_STATS_NFIELDS_LOC)*nptot*PARTICLES_LP_STATS_MAX_EVENTS;

  dump_arr = new double[SIZE];

  #pragma acc enter data create(dump_arr[:SIZE])


  int shkTot = 0;
  int nshMax = -1;
  #pragma acc parallel loop present(data, dump_arr[:SIZE]) \
  copy(shkTot, nshMax) reduction(max:nshMax)
  for (uint64_t np = 0; np < nptot; np++){
    int nshk = static_cast<int>(data->pc->color(np)[N]);
    nshMax = std::max(nshMax, nshk);
    if (nshk > 0){
      #pragma acc loop seq
      for(int shk = 0; shk < MIN(nshk,PARTICLES_LP_STATS_MAX_EVENTS); shk++){
        int jj = 0, n;
        #pragma acc atomic capture
        n = shkTot++;
        #pragma acc loop seq
        for(int field = 0; field < PARTICLES_LP_STATS_NFIELDS; field++){
          dump_arr[n*PARTICLES_LP_STATS_NFIELDS_LOC + (jj++)] = data->pc->color(np)[shk*PARTICLES_LP_STATS_NFIELDS + field];
        }
        dump_arr[n*PARTICLES_LP_STATS_NFIELDS_LOC + (jj++)] = nshk;
      }
      // dump_arr[n*PARTICLES_LP_STATS_NFIELDS + (jj++)] = nshk;
    }
  }

  // WRITING data according to particles_lp_spectral_update.cpp line 1400 
  #pragma acc exit data copyout(dump_arr[:SIZE])

  if (nCalls == 0){
    FILE *fp;
    fp = fopen(fname, "w");
    fprintf(fp, "#ID;sigma;beta;Emax;Btot;Brec;p;diffEfl;dEpsilon;normType;DealtaE;finalE;t_rec;nshk;nstep\n");
    fclose(fp);
    printf("%s\n", fname);
    fp = fopen(metaFname, "w");
    fclose(fp);
  }

  int nLines = 0;
  FILE *fp;
  fp = fopen(fname, "a");
  for (int n = 0; n < shkTot; n++){
    int jj = 0;
    for(int field = 0; field < PARTICLES_LP_STATS_NFIELDS_LOC; field++){
      fprintf(fp, "%f;", dump_arr[n*PARTICLES_LP_STATS_NFIELDS_LOC + (jj++)]);
    }
    fprintf(fp, "%d\n", g_stepNumber);
    nLines++;
  }

  fclose(fp);

  fp = fopen(metaFname, "a");
  fprintf(fp, "%ld;%d;%d\n",nCalls, g_stepNumber, nLines);

  fclose(fp);


  #pragma acc parallel loop present(data)
  for (uint64_t np = 0; np < nptot; np++){
    data->pc->color(np)[N + 1] += data->pc->color(np)[N];
    for(int nc = 0; nc < N + 1; nc++){
      data->pc->color(np)[nc] = 0.0;
    }
  }
  delete [] dump_arr;

#endif
}


////////////////////////////////////////////////////////////////////////
/*
 *
*/
////////////////////////////////////////////////////////////////////////
long int GetFileSize(FILE *fp)
{
  long fsize = 0;

  fseek(fp,0,SEEK_END);
  fsize = ftell(fp);
  fseek(fp,0,SEEK_SET);//reset stream position!!

  return fsize;
}

////////////////////////////////////////////////////////////////////////
/*
 *
*/
////////////////////////////////////////////////////////////////////////
char *lastline(char *filepath)
{
  FILE *fp;
  char buff[4096+1];
  int size,i;
  long fsize;
  if(NULL==(fp=fopen(filepath, "r"))){
      perror("file cannot open at lastline");
      return NULL;
  }
  fsize= -1L*GetFileSize(fp);
  if(size=fseek(fp, std::max(fsize, -4096L), SEEK_END)){
      perror("cannot seek");
      exit(1);
  }
  size=fread(buff, sizeof(char), 4096, fp);
  fclose(fp);
  buff[size] = '\0';
  i=size-1;
  if(buff[i]=='\n'){
      buff[i] = '\0';
  }
  while(i >=0 && buff[i] != '\n')
      --i;
  ++i;
  return strdup(&buff[i]);
}


void truncate_file(const char *fname, int n) {

  std::string line;  // Declaring a string variable to store each line of text
  int lineCount = 0;  // Initializing a variable to count lines
  {
    std::ifstream inputStream(fname);
    while (std::getline(inputStream, line)) {  // Loop through each line in the file
      lineCount++;  // Incrementing line count for each line read
    }
    inputStream.close();
  }

  char fnameTmp[128];
  strcpy(fnameTmp, fname);

  char Tmp[128];
//  sprintf(Tmp, "tmp");
  snprintf(Tmp, sizeof(Tmp), "tmp");

  strcat(fnameTmp, Tmp);
  std::ifstream inputStream(fname);
  std::ofstream outputStream(fnameTmp);

  std::getline(inputStream, line);
  int count = 0;
  std::string tmp;
  do{
    std::getline(inputStream,tmp);
    outputStream << line << '\n';
    count++;
    line.swap(tmp);
  } while ( count < lineCount - n);
  outputStream.close();
  inputStream.close();
  std::rename(fnameTmp, fname);
}

#endif // PARTICLES_INIT_NEW == YES


