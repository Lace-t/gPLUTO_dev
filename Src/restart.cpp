////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Restart PLUTO from binary or HDF5 double precision data files.
 
  This file collects the necessary functions for restarting PLUTO 
  from a double precision binary or HDF5 file in the static grid
  version of the code.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   June 1, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if defined(STAGGERED_MHD) && defined(PARALLEL)
static void RestartStaggeredPatch (const Ary3D (&Vs)[6], Grid *);
#endif

#define MAX_BUF_SIZE 512
////////////////////////////////////////////////////////////////////////
/*!
  Read input binary / hdf5 data.
  
  \param [in] data       pointer to Data structure
  \param [in] ini        pointer to Runtime structure
  \param [in] Dts        pointer to timeStep structure
  \param [in] nrestart   number of restart file
  \param [in] type       specifies the output data type (type should be 
                         either DBL_OUTPUT or DBL_H5_OUTPUT).
  \param [in] grid       pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void RestartFromFile (Data *data, Runtime *ini, timeStep *Dts, 
                      int nrestart, int type, Grid *grid)
{ 
  int   nx1_tot = grid->ntot[IDIR];
  int   nx2_tot = grid->ntot[JDIR];
  int   nx3_tot = grid->ntot[KDIR];
  int   SZ_stag[3] = {SZ_stagx, SZ_stagy, SZ_stagz};
  int   nv, single_file, origin, nlines=0;
  int dummy;  // dummy variable to suppress unused variable warning
  int   swap_endian=0;
  char  fname[4*MAX_STRING_SIZE], fout[MAX_STRING_SIZE+16];
  char  str[MAX_STRING_SIZE+16];
  double  dbl;
  void   *Vpt;
  RBox   box, box_tot;
  Output *output;
  FILE   *fbin;

  RBoxDefine (grid->lbeg, grid->lend, &box);
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &box_tot);

  // ------------------------------------------------------
  // 1. Get the pointer to the output format specified by
  //    "type"
  // ------------------------------------------------------

  for (nv = 0; nv < MAX_OUTPUT_TYPES; nv++){
    output = ini->output + nv;
    if (output->type == type) break;
  }

  // ------------------------------------------------------
  // 2. Compare the endianity of the restart file (by reading
  //    the corresponding entry in dbl.out or dbl.h5.out) 
  //    with that of the current architecture.
  //    Turn swap_endian to 1 if they're different.
  // ------------------------------------------------------

  if (g_rank == 0){
    if (type == DBL_OUTPUT) {
      snprintf (fout, sizeof(fout), "%s/dbl.out",ini->output_dir);
      fbin = fopen (fout, "r");
    } else if (type == DBL_H5_OUTPUT) {
      snprintf (fout, sizeof(fout), "%s/dbl.h5.out",ini->output_dir);
      fbin = fopen (fout, "r");
    }
    if (fbin == NULL){
      QuitPLUTO("from RestartFromFile(): cannot find dbl.out or dbl.h5.out",1);
    }
    while (fgets(str, sizeof(str), fbin) != 0) nlines++;  // Count lines in dbl.out

    rewind(fbin);

    if (nrestart > nlines-1){
      string msg =   "from RestartFromFile(): output #" + to_string(nrestart) 
                   + " does not exist in file " + fout;
      QuitPLUTO(msg, 1);
    }
    origin = (nrestart >= 0 ? nrestart:(nlines+nrestart));
    for (nv = origin; nv--;   ) while ( fgetc(fbin) != '\n'){}
    dummy = fscanf(fbin, "%d  %lf  %lf  %d  %s  %s\n",
                          &nv, &dbl, &dbl, &nv, str, str);
    if ( (!strcmp(str,"big")    &&  IsLittleEndian()) ||
         (!strcmp(str,"little") && !IsLittleEndian())) {
      swap_endian = 1;
      print ("> RestartFromFile(): endianity is reversed\n");
    }
    fclose(fbin);

  }
  #ifdef PARALLEL
  MPI_Bcast (&swap_endian, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  // ------------------------------------------------------
  // 3. Read restart.out and get Restart structure
  // ------------------------------------------------------

  RestartGet (data, ini, Dts, nrestart, type, swap_endian);

  if (type == DBL_H5_OUTPUT){
    #ifdef USE_HDF5
    ReadHDF5 (output, grid);
    #endif
  }

  // ------------------------------------------------------
  // 4. For .dbl output, read data from disk
  // ------------------------------------------------------

  if (type == DBL_OUTPUT){
    print ("> Restarting from file #%d (dbl)\n",output->nfile);
    single_file = strcmp(output->mode,"single_file") == 0;
    if (single_file){
      #if USE_AL == NO 
      snprintf (fname, sizeof(fname), "%s/data.%04d.dbl", output->dir, output->nfile);

      // 4a. Check if the file exists

      int file_exists = 0;
      if (g_rank == 0) {
        file_exists = std::filesystem::exists(fname);
        if (!file_exists){
          string msg =   "from RestartFromFile(): file " + string(fname) 
                      + " does not exist";
          QuitPLUTO(msg,1);
        }
      }

      // 4b. Read buffer

      FileOffset offset = 0;

      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;

        int stag_dir = output->stag_var[nv];
        if (stag_dir == -1) {  //  cell-centered data
          Vpt = (void *)output->Vc_buf;
        } else  { // Staggered data 
          Vpt = (void *)output->Vs_buf[stag_dir];
        }

        FileReadData(Vpt, sizeof(double), 
                        stag_dir, swap_endian, grid, fname, offset);

        // Advance offset by this variable's global block size
        int ioff = (stag_dir == 0);
        int joff = (stag_dir == 1);
        int koff = (stag_dir == 2);
        offset += (FileOffset)(grid->nint_glob[IDIR] + ioff)
                            *(grid->nint_glob[JDIR] + joff)
                            *(grid->nint_glob[KDIR] + koff)*sizeof(double);

        ReadOutputBuffer (nv, output, &box);
      }
      #else   // Old version with ArrayLib
      int  sz;
      
      snprintf (fname, sizeof(fname), "%s/data.%04d.dbl", output->dir, output->nfile);
      #ifdef PARALLEL
      long long offset;
      offset = 0;
      #else
      fbin = FileOpen (fname, 0, "r");
      #endif
      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;

        if      (output->stag_var[nv] == -1) {  //  cell-centered data
          sz = SZ;
          Vpt = (void *)output->Vc_buf_AL[0][0];
        } else  { // Staggered data 
          int stag_dir = output->stag_var[nv];
          sz  = SZ_stag[stag_dir];
          Vpt = (void *)output->Vs_buf_AL[stag_dir][0][0];
        }

        #ifdef PARALLEL
        fbin = FileOpen (fname, sz, "r");
        AL_Set_offset(sz, offset);
        #endif
        FileReadData_AL (Vpt, sizeof(double), sz, fbin,
                      output->stag_var[nv], swap_endian, grid);
        #ifdef PARALLEL
        offset = AL_Get_offset(sz);
        FileClose(fbin, sz);
        #endif

        ReadOutputBuffer_AL (nv, output, &box_tot);
      }
      #ifndef PARALLEL
      FileClose(fbin, sz);
      #endif
  #endif
    }else{   // Multiple files
      #if USE_AL == NO
      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;

        snprintf(fname, sizeof(fname), "%s/%s.%04d.%s",
                  output->dir, output->var_name[nv], output->nfile, output->ext);

        int file_exists = 0;
        if (g_rank == 0) {
          file_exists = std::filesystem::exists(fname);
          if (!file_exists){
            string msg =   "from RestartFromFile(): file " + string(fname) 
                        + " does not exist";
            QuitPLUTO(msg,1);
          }
        }

        int stag_dir = output->stag_var[nv];
        if (stag_dir == -1) {  // Cell-centered data
          Vpt = (void *)output->Vc_buf;
        } else { // Staggered data 
          Vpt = (void *)output->Vs_buf[stag_dir];
        }

        FileReadData(Vpt, sizeof(double), 
                          stag_dir, swap_endian, grid, fname, 0);

        ReadOutputBuffer(nv, output, &box);
      }
      #else         // Old version with ArrayLib
      int  sz;

      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        snprintf (fname, sizeof(fname), "%s/%s.%04d.%s", output->dir, 
                  output->var_name[nv], output->nfile, output->ext);

        if      (output->stag_var[nv] == -1) {  // Cell-centered data
          sz = SZ;
          Vpt = (void *)output->Vc_buf_AL[0][0];
        } else { // Staggered data 
          int stag_dir = output->stag_var[nv];
          sz  = SZ_stag[stag_dir];
          Vpt = (void *)output->Vs_buf_AL[stag_dir][0][0];
        }
        fbin = FileOpen (fname, sz, "r");
        FileReadData_AL (Vpt, sizeof(double), sz, fbin,
                        output->stag_var[nv], swap_endian, grid);
        FileClose (fbin, sz);

        ReadOutputBuffer_AL (nv, output, &box_tot);
      }
  #endif
    }
  }

  data->run_config->time = Dts->time;
  // #pragma acc update device(data->run_config)
  // #pragma acc update device(data->Dts)
  // #pragma acc update device(data)
  #ifdef FARGO
  FARGO_Restart(data, output->dir, output->nfile, swap_endian, grid);
  #endif

  // ------------------------------------------------------
  // 5. IMPORTANT (in PARALLEL): owing to the way ArrayLib is 
  //    conceived, the leftmost staggered point of the array will 
  //    not be properly set. 
  //    Example, for a processor with 4 cells:
  //
  //     |---|---|---|---|
  //     A   B   C   D   E
  //
  //    "A" will be undefined at this stage. 
  //    We solve this issue by calling RestartStaggeredPatch()
  //    until AL is redesigned.
  // ------------------------------------------------------

  #if USE_AL == YES
  #if defined(STAGGERED_MHD) && defined(PARALLEL)
  RestartStaggeredPatch (data->Vs, grid);
  #endif
  #endif

  // ------------------------------------------------------
  // 6. Convert to conservative (for standard 2nd order)
  //    or primitive (for high order)
  // ------------------------------------------------------

  #ifndef HIGH_ORDER

  #if CONS_DBL == YES
  ConsToPrim3D (data, 0, grid, &box);
  #else
  PrimToCons3D (data, 0, grid, &box);
  #endif

  #endif
  #ifdef HIGH_ORDER
  #pragma acc enter data copyin (box)  
  ConsToPrim3D (data, 0, grid, &box);
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost)
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc
  #pragma acc exit data delete(box)
  #endif

}

static int counter = -1;
////////////////////////////////////////////////////////////////////////
/*!
  Read file \c restart.out and collect information needed for 
  (potential) later restarts.

  \param [in] data         pointer to Data structure
  \param [in] ini          pointer to Runtime structure
  \param [in] Dts          pointer to timeStep structure
  \param [in] nrestart     number of restart file
  \param [in] out_type     specifies the output data type (type should be 
                           either DBL_OUTPUT or DBL_H5_OUTPUT).
  \param [in] swap_endian  flag to reverse endiantiy
*/
////////////////////////////////////////////////////////////////////////
void RestartGet (Data *data, Runtime *ini, timeStep *Dts, 
                 int nrestart, int out_type, int swap_endian)
{
  int  origin, n, k;
  char fout[512];
  int dummy;
  Restart restart;
  FILE *fr;

  if (nrestart < 0){
    QuitPLUTO ("from RestartGet(): negative restart file temporarily disabled",1);
  }

  // 1. Open "restart.out" and scan line by line until
  //    the output type specified by "out_type" has
  //    nfile = nrestart.
  //    counter will contain the line number where this
  //    occurs. Processor 0 does the actual reading.

  if (g_rank == 0) {

    snprintf (fout, MAX_BUF_SIZE, "%s/restart.out",ini->output_dir);
    fr = fopen (fout, "rb");
    if (fr == NULL){
      QuitPLUTO("from RestartGet(): cannot find restart.out",1);
    }

    origin = (nrestart < 0 ? SEEK_END:SEEK_SET);
    k = 0;
    while (counter == -1){
      if (feof(fr)){
        QuitPLUTO ("from RestartGet(): end of file encountered",1);
      }
      fseek (fr, k*sizeof(Restart), origin);
      dummy = fread (&restart, sizeof (Restart), 1, fr);
      for (n = 0; n < MAX_OUTPUT_TYPES; n++){
        if (swap_endian) SWAP_VAR(restart.nfile[n]);
        if (swap_endian) SWAP_VAR(restart.nfileParticles[n]);
        if (ini->output[n].type == out_type && restart.nfile[n] == nrestart){
          counter = k;
        }
      }
      k++;
    }
    fclose(fr);
    if (swap_endian){
      SWAP_VAR(restart.t);
      SWAP_VAR(restart.dt);
      SWAP_VAR(restart.nstep);
      for (n = 0; n < MAX_RESTART_PARAMS; n++) SWAP_VAR(restart.param[n]);
    }
  }

#ifdef PARALLEL
  MPI_Bcast (&restart, sizeof (Restart), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  Dts->time = restart.t;
  Dts->dt   = restart.dt;
  g_stepNumber = restart.nstep;
  // #pragma acc update device(g_stepNumber)

  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    ini->output[n].nfile = restart.nfile[n];
    ini->output[n].nfileParticles = restart.nfileParticles[n];
  }
  // Check that tstop > time 

  if (RuntimeGet()->tstop < Dts->time){
    string msg =   "from RestartGet(): tstop = " + to_string(RuntimeGet()->tstop) 
                 + " < Dts.t = " + to_string (Dts->time);
    QuitPLUTO(msg,1);
  }

  // Copy restart parameters back on data->>
  for (n = 0; n < MAX_RESTART_PARAMS; n++) {
    data->restartParam[n] = restart.param[n];
//    print ("  restart parameter #%dd = %12.6e\n", n, restart.param[n]);
  }

}
////////////////////////////////////////////////////////////////////////
/*!
  Write restart file \c restart.out

  \param [in] ini          pointer to Runtime structure
  \param [in] restartParam pointer to array of parameters
  \param [in] Dts          pointer to timeStep structure
*/
////////////////////////////////////////////////////////////////////////
void RestartDump (Runtime *ini, double *restartParam, timeStep *Dts)
{
  int n;
  char fout[512];
  Restart restart;
  FILE *fr;

  // 1. Define restart structure elements here

  restart.t     = Dts->time;
  restart.dt    = Dts->dt;
  restart.nstep = g_stepNumber;
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    restart.nfile[n]          = ini->output[n].nfile;
    restart.nfileParticles[n] = ini->output[n].nfileParticles;
  }  

  for (n = 0; n < MAX_RESTART_PARAMS; n++){
    restart.param[n] = restartParam[n];
  }
  
  // 2. Dump structure to disk

  counter++;
  if (g_rank == 0) {   // Only processor 0 does the writing 
    snprintf (fout,MAX_BUF_SIZE, "%s/restart.out",ini->output_dir); // File name 
    if (counter == 0) {
      fr = fopen (fout, "wb");
    }else {
      fr = fopen (fout, "r+b");
      fseek (fr, counter*sizeof(Restart), SEEK_SET); 
    }

    fwrite (&restart, sizeof(Restart), 1, fr);
    fclose(fr);
  }
}

#if defined(STAGGERED_MHD) && defined(PARALLEL)
////////////////////////////////////////////////////////////////////////
/*!
  Transfer staggered normal component between processors, in
  order to fill leftmost interface value:
  
  A1  A2  A3  A4  A5
  |---|---|---|---|
  B1  B2  B3  B4  B5
  |---|---|---|---|
  C1  C2  C3  C4  C5
  |---|---|---|---|
  
  Here A1, B1 and C1 (which are undefined owing to the way AL handles 
  staggered array), must be overwritten with the value coming from 
  the left processor, e.g. 
  
  A5 --> B1
  B5 --> C1
  
  \param [in,out] Vs  pointer to staggered array
  \param [in]  grid   pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
void RestartStaggeredPatch (const Ary3D (&Vs)[6], Grid *grid)
{
  int  i,j,k;
  int  par_dim[3] = {0, 0, 0};
  MPI_Status status;

  static int **nrnks;
  static int ***nrnks_edg;
  static int ***nrnks_crn;

  static double *recv_bufL;
  static double *send_bufR;

  int size_xy = grid->ntot[IDIR]*grid->ntot[JDIR];
  int size_yz = grid->ntot[JDIR]*grid->ntot[KDIR];
  int size_zx = grid->ntot[KDIR]*grid->ntot[IDIR];

  int size = MAX(size_xy, size_yz);
  size = MAX(size, size_zx);

  // 0a. Allocate arrays

  if (nrnks == NULL){
    nrnks    = ARRAY_2D(3, 2, int);
    nrnks_edg = ARRAY_3D(3, 2, 2, int);
    nrnks_crn = ARRAY_3D(2, 2, 2, int);

    recv_bufL = ARRAY_1D(size, double);
    send_bufR = ARRAY_1D(size, double);

    #pragma acc enter data create(send_bufR[:size], recv_bufL[:size])

    GetNeighbourRanks(grid, nrnks, nrnks_edg, nrnks_crn);
  }

  // 0b. Check the number of processors in each direction

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

  // 1. Loop over dimensions and exchange only normal 
  //    component of staggered field(s).

  int staggered_fields = DIMENSIONS;
  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  staggered_fields += DIMENSIONS;
  #endif

  for (int ksf = 0; ksf < staggered_fields; ksf++){

    int dir = ksf%DIMENSIONS;

    if (!par_dim[dir]) continue;  // Direction is not parallel, skip

    int l;
    int beg[3], end[3];
    long ind, countR;
    int rankL = nrnks[dir][0];
    int rankR = nrnks[dir][1];
    RBox send_boxR;
    RBox recv_boxL; 
    RBox *box;

    // 1a. Define boxes 

    for (l = 0; l < 3; l++) {  // Initialize to whole domain
      beg[l] = 0;
      end[l] = grid->ntot[l] - 1;
    }
    beg[dir] = grid->lend[dir]+1;
    end[dir] = beg[dir];
    RBoxDefine (beg, end, &send_boxR);

    beg[dir] = grid->lbeg[dir];
    end[dir] = beg[dir];
    RBoxDefine (beg, end, &recv_boxL);

    // 1b. Pack cell-centered data 

    box = &send_boxR;
    BOX_LOOP(*box,i,j,k) {
      ind = BOX_INDEX(*box,i,j,k);
      send_bufR[ind] = Vs[ksf](i,j,k);
    }
    countR = box->nx1*box->nx2*box->nx3; 

    // 1c. Send / recv data

    MPI_Sendrecv (send_bufR, countR, MPI_DOUBLE, rankR, 2,
                  recv_bufL, countR, MPI_DOUBLE, rankL, 2,  MPI_COMM_WORLD, &status);

    // 1d. Unpack staggered data

    if (rankL != MPI_PROC_NULL){
      box = &recv_boxL;
      BOX_LOOP(*box,i,j,k){
        ind = BOX_INDEX(*box,i,j,k);
        Vs[ksf](i,j,k) = recv_bufL[ind];
      }
    }
  } //  End loop on dir
}
#endif
