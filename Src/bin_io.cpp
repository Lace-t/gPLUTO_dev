////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Functions for handling binary I/O.

  This file provides a number of handy functions for opening, reading
  and writing binary files using single or double precision in serial
  or parallel mode.
  It is employed by the following output formats: .dbl, .flt and .vtk.
  
  In parallel mode these functions work as wrappers to the actual
  parallel implementations contained in AL_io.c.
 
  Pointer to data array must be cast into (void *) and are assumed to
  start with index 0 for both cell-centered and staggered data arrays.
  This means that if V3D is a 3D array then:
  - for a cell-centered array, <tt>V = (void *) V3D[0][0]</tt>
  - for a x-staggered array,   <tt>V = (void *)(V3D[0][0]-1)</tt>
  - for a y-staggered array,   <tt>V = (void *)(V3D[0][-1])</tt>
  - for a z-staggered array,   <tt>V = (void *)(V3D[-1][0])</tt>

  \authors A. Mignone (andrea.mignone@unito.it)\n

  \date   May 28, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*
  Convert the a double-precision 3D array into single precision.
  Swap endianity if swap_endian == 1.
  
  \param [in] Vdbl         pointer to a 3D double precision aray
  \param [in] unit         a multiplicative constant typically used 
  to write in c.g.s units.
  \param [in] swap_endian  when set to 1, swap endianity during the 
  conversion.
  \return a pointer to a 3D array in single precision.
*/
////////////////////////////////////////////////////////////////////////
float *Convert_dbl2flt (double *Vdbl, double unit, int swap_endian, Grid *grid)
{
  int i, j, k, ind;
  float  flt;
  static float *Vflt;
  
  int  ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int  jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int  kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  RBox box;
  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);

  if (Vflt == NULL) {
    Vflt = ARRAY_1D(grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR], float);
  }

  if (!swap_endian){
    BOX_LOOP(box, i,j,k) {
      ind = BOX_INDEX(box, i, j, k);
      Vflt[ind] = (float)(Vdbl[ind]*unit);
    }
  }else{
    BOX_LOOP(box, i,j,k) {
      ind = BOX_INDEX(box, i, j, k);
      flt = (float)(Vdbl[ind]*unit);
      SWAP_VAR(flt);
      Vflt[ind] = flt;
    }
  }

  return Vflt;
}

////////////////////////////////////////////////////////////////////////
/*
  Convert the a double-precision 3D array into single precision.
  Swap endianity if swap_endian == 1.
  
  \param [in] Vdbl         pointer to a 3D double precision aray
  \param [in] unit         a multiplicative constant typically used 
  to write in c.g.s units.
  \param [in] swap_endian  when set to 1, swap endianity during the 
  conversion.
  \return a pointer to a 3D array in single precision.
*/
////////////////////////////////////////////////////////////////////////
float ***Convert_dbl2flt_AL (double ***Vdbl, double unit, int swap_endian, Grid *grid)
{
  int i, j, k;
  float  flt;
  static float ***Vflt;
  
  int  ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int  jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int  kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  if (Vflt == NULL) {
    Vflt = ARRAY_3D(grid->ntot[KDIR], grid->ntot[JDIR], grid->ntot[IDIR], float);
  }

  if (!swap_endian){
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      Vflt[k][j][i] = (float)(Vdbl[k][j][i]*unit);
    }}}
  }else{
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      flt = (float)(Vdbl[k][j][i]*unit);
      SWAP_VAR(flt);
      Vflt[k][j][i] = flt;
    }}}
  }

  return Vflt;
}

////////////////////////////////////////////////////////////////////////
/*!
  Close file.
  
  \param [in] fbin     pointer to the FILE that needs to be closed (serial
  mode only)
  \param [in] sz       the distributed array descriptor for parallel mode
  
  \return Returns 0 on success
*/
////////////////////////////////////////////////////////////////////////
int FileClose (FILE *fbin, int sz) 
{
  #ifdef PARALLEL 
#if USE_AL == YES
  AL_File_close(sz); 
#endif
  #else
  fclose (fbin);
  #endif
  return(0);
}

////////////////////////////////////////////////////////////////////////
/*!
  Open a file for write/reading in binary mode.
  
  \param [in] filename a valid file name
  \param [in] sz       the distributed array descriptor. This parameter
  replaces dsize in parallel mode
  \param [in] mode     a string giving the opening mode (only "w" or 
  "r" are allowed)
  
  \return The pointer to the file.
*/
////////////////////////////////////////////////////////////////////////
FILE *FileOpen (char *filename, int sz, const char *mode)  
{
  FILE *fp;

  // ----------------------------------------------
  // when reading, check if file exists
  // ----------------------------------------------
   
  if (strcmp(mode,"r") == 0 && g_rank == 0){
    fp = fopen(filename, "rb");
    if (fp == NULL){
      QuitPLUTO ("from FileOpen(): file " + string(filename) + "does not exists", 1);
    }
    fclose(fp);
  }
  
  // ------------------------------------------------ 
  // file exists, keep going
  // ------------------------------------------------
  
#ifdef PARALLEL
#if USE_AL == YES
/*  CANNOT USE IT HERE SINCE FILE IS OPENED AND CLOSED FOR APPENDING SEVERAL TIMES
  if (strcmp(mode,"w") == 0){
print ("Deleting file...\n");
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                  MPI_INFO_NULL, &fh);
    MPI_File_close(&fh);
  }  
*/      
  AL_File_open(filename, sz); 
#endif
  return NULL;
#else
  if      (strcmp(mode,"w") == 0) fp = fopen(filename, "wb");
  else if (strcmp(mode,"r") == 0) fp = fopen(filename, "rb");
  else if (strcmp(mode,"a") == 0) fp = fopen(filename, "ab");
  if (fp == NULL){
    printLog ("! FileOpen(): cannot find file %s\n",filename);
    QuitPLUTO("",1);
  }
  return (fp);
#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Read a 3D grid data array V[k][j][i] from a binary file. 
  
  \param [in] V      pointer to a 3D array,
  V[k][j][i] -->  V[i + NX1*j + NX1*NX2*k]. Must
  start with index 0
  \param [in] dsize  the size of the each buffer element   
  (sizeof(double) or sizeof(float)). This
  parameter is used only in serial mode.
  \param [in] sz     the distributed array descriptor. This parameter
  replaces dsize in parallel mode
  \param [in] fl     a valid FILE pointer
  \param [in] istag  a flag to identify cell-centered (istag = -1) or
  staggered field data (istag = 0,1 or 2 for staggering 
  in the x1, x2 or x3 directions) 
  \param [in] swap_endian  a flag for swapping endianity
  
  \return This function has no return value.
  
  \remark The data array V is assumed to start \a always with index 0
  for  both cell-centered and staggered arrays.
*/
////////////////////////////////////////////////////////////////////////
void FileReadData_AL (void *V, size_t dsize, int sz, FILE *fl, int istag,
                   int swap_endian, Grid *grid) 
{
#if USE_AL == YES
  int i, j, k;
  int ioff, joff, koff;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  size_t dummy;  // dummy variable to suppress unused variable warning */
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];

  // ---------------------------------------
  // parallel reading handled by ArrayLib
  // ---------------------------------------

#ifdef PARALLEL 
  AL_Read_array (V, sz, istag);
#else

  // ---------------------------------------
  // serial reading
  // ---------------------------------------

  char * Vc = (char *) V;
  ioff = (istag == 0); 
  joff = (istag == 1); 
  koff = (istag == 2);

  for (k = kbeg; k <= kend + koff; k++) {
  for (j = jbeg; j <= jend + joff; j++) {
    i = ibeg + (nx1_tot + ioff)*(j + (nx2_tot + joff)*k);
    dummy = fread (Vc + i*dsize, dsize, grid->nint[IDIR] + ioff, fl);
  }}
#endif

  // -----------------------------------
  // now swap endian if necessary
  // -----------------------------------

  if (swap_endian){
    double *Vd;
    Vd = (double *) V;
    ioff = (istag == 0); 
    joff = (istag == 1); 
    koff = (istag == 2);

    for (k = kbeg-koff; k <= kend + koff; k++) {
    for (j = jbeg-joff; j <= jend + joff; j++) {
    for (i = ibeg-ioff; i <= iend + ioff; i++) {
      SWAP_VAR(Vd[i + (nx1_tot + ioff)*(j + (nx2_tot + joff)*k)]);
    }}}     
  }
#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Write a 3D grid data array V[k][j][i] to a binary file. 
  
  \param [in] V      single pointer to a 3D array,
  V[k][j][i] -->  V[i + NX1*j + NX1*NX2*k]. Must
  start with index 0
  \param [in] dsize  the size of the each buffer element   
  (sizeof(double) or sizeof(float)). This
  parameter is used only in serial mode.
  \param [in] sz     the distributed array descriptor. This parameter
  replaces dsize in parallel mode
  \param [in] fl     a valid FILE pointer
  \param [in] istag  a flag to identify cell-centered (istag = -1) or
  staggered field data (istag = 0,1 or 2 for staggering 
  in the x1, x2 or x3 directions) 
  
  \return This function has no return value.
  
  \remark The data array is assumed to start \a always with index 0 for
  both cell-centered and staggered arrays.
*/
////////////////////////////////////////////////////////////////////////
void FileWriteData(void *V, size_t dsize, FileDataType  mpi_type,
                   int istag, Grid *grid, char *fname, FileOffset  offset)

/*
void FileWriteData(void *V, size_t dsize, char *fname,
                   int istag, long int offset, Grid *grid)
void FileReadData (void *V, size_t dsize, char *fname, 
                   int istag, int swap_endian, long int offset, Grid *grid)
*/
{
  int ioff = (istag == 0);
  int joff = (istag == 1);
  int koff = (istag == 2);
  int nx_loc  = grid->nint[IDIR] + ioff;
  int ny_loc  = grid->nint[JDIR] + joff;
  int nz_loc  = grid->nint[KDIR] + koff;
  
  #ifdef PARALLEL

  int nx_glob = grid->nint_glob[IDIR] + ioff;
  int ny_glob = grid->nint_glob[JDIR] + joff;
  int nz_glob = grid->nint_glob[KDIR] + koff;

  int gi0 = grid->beg[IDIR] - grid->nghost[IDIR];
  int gj0 = grid->beg[JDIR] - grid->nghost[JDIR];
  int gk0 = grid->beg[KDIR] - grid->nghost[KDIR];

  int gsizes[3] = {nz_glob, ny_glob, nx_glob};
  int lsizes[3] = {nz_loc,  ny_loc,  nx_loc };
  int starts[3] = {gk0,     gj0,     gi0    };

  MPI_Datatype filetype;
  MPI_File     fh;
  MPI_Status   status;

  MPI_Type_create_subarray(3, gsizes, lsizes, starts,
                           MPI_ORDER_C, mpi_type, &filetype);
  MPI_Type_commit(&filetype);

  MPI_File_open(MPI_COMM_WORLD, fname,
                MPI_MODE_WRONLY | MPI_MODE_APPEND,
                MPI_INFO_NULL, &fh);

  // offset=0: use APPEND position (dbl case)   
  // offset>0: use broadcast ftell value (vtk case) 
  MPI_Offset start_pos;
  if (offset == 0) {
    MPI_File_get_position(fh, &start_pos);
  } else {
    start_pos = offset;
  }

  MPI_File_set_view(fh, start_pos, mpi_type, filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, V, nx_loc*ny_loc*nz_loc, mpi_type, &status);

  MPI_Type_free(&filetype);
  MPI_File_close(&fh);

  #else

  FILE *fbin = fopen(fname, "ab");
  fwrite(V, dsize, nx_loc*ny_loc*nz_loc, fbin);
  fclose(fbin);

  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Read a 3D grid data array from a binary file into a 1D flat buffer
  containing only active computational cells.
  No AL distributed array descriptor is used here.
  In parallel mode, it uses MPI-IO subarray views to read from the file.
*/
////////////////////////////////////////////////////////////////////////
void FileReadData(void *V, size_t dsize, int istag, int swap_endian, 
                  Grid *grid, char *fname, FileOffset offset)
{
  int ioff    = (istag == 0);
  int joff    = (istag == 1);
  int koff    = (istag == 2);

  int nx_loc  = grid->nint[IDIR] + ioff;
  int ny_loc  = grid->nint[JDIR] + joff;
  int nz_loc  = grid->nint[KDIR] + koff;

#ifdef PARALLEL
  int nx_glob = grid->nint_glob[IDIR] + ioff;
  int ny_glob = grid->nint_glob[JDIR] + joff;
  int nz_glob = grid->nint_glob[KDIR] + koff;

  int gi0 = grid->beg[IDIR] - grid->nghost[IDIR];
  int gj0 = grid->beg[JDIR] - grid->nghost[JDIR];
  int gk0 = grid->beg[KDIR] - grid->nghost[KDIR];

  int gsizes[3] = {nz_glob, ny_glob, nx_glob};
  int lsizes[3] = {nz_loc,  ny_loc,  nx_loc };
  int starts[3] = {gk0,     gj0,     gi0    };

  MPI_Datatype filetype;
  MPI_File     fh;
  MPI_Status   status;
  MPI_Offset   start_pos;

  MPI_Type_create_subarray(3, gsizes, lsizes, starts,
                           MPI_ORDER_C, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  MPI_File_set_view(fh, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, V, nx_loc*ny_loc*nz_loc, MPI_DOUBLE, &status);
  
  MPI_Type_free(&filetype);
  MPI_File_close(&fh);

#else  // Serial read
  FILE *fbin = fopen(fname, "rb");

  fseek(fbin, offset, SEEK_SET);
  fread(V, dsize, nx_loc*ny_loc*nz_loc, fbin);
  fclose(fbin);
#endif

  // -----------------------------------------------------------
  // Swap endianity if needed - operates on the ghosted array.
  // ----------------------------------------------------------- 

  if (swap_endian) {
    double *Vd  = (double *) V;
    for (long int i = 0; i < (nx_loc*ny_loc*nz_loc); i++){
      SWAP_VAR(Vd[i]);
    }
  }
}


////////////////////////////////////////////////////////////////////////
/*!
  Write a 3D grid data array V[k][j][i] to a binary file. 
  
  \param [in] V      single pointer to a 3D array,
  V[k][j][i] -->  V[i + NX1*j + NX1*NX2*k]. Must
  start with index 0
  \param [in] dsize  the size of the each buffer element   
  (sizeof(double) or sizeof(float)). This
  parameter is used only in serial mode.
  \param [in] sz     the distributed array descriptor. This parameter
  replaces dsize in parallel mode
  \param [in] fl     a valid FILE pointer
  \param [in] istag  a flag to identify cell-centered (istag = -1) or
  staggered field data (istag = 0,1 or 2 for staggering 
  in the x1, x2 or x3 directions) 
  
  \return This function has no return value.
  
  \remark The data array is assumed to start \a always with index 0 for
  both cell-centered and staggered arrays.
*/
////////////////////////////////////////////////////////////////////////
void FileWriteData_AL (void *V, size_t dsize, int sz, FILE *fl, int istag, Grid *grid)
{
#if USE_AL == YES
#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  AL_Write_array (V, sz, istag);
  return;
#else
  int i,j,k;
  int ioff, joff, koff;
  int ibeg = grid->lbeg[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  char *Vc = (char *) V;
  ioff = (istag == 0); 
  joff = (istag == 1); 
  koff = (istag == 2);

  for (k = kbeg; k <= kend + koff; k++) {
  for (j = jbeg; j <= jend + joff; j++) {
    i = ibeg + (nx1_tot + ioff)*(j + (nx2_tot + joff)*k);
    fwrite (Vc + i*dsize, dsize, grid->nint[IDIR] + ioff, fl);
  }}
#endif
#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Delete a file (parallel mode only).
  Useful on some architecture where MPI_File_open will just overwrite
  the file without erasing it.
  
  
  \param [in] filename a valid file name
  
  \return The pointer to the file.
*/
////////////////////////////////////////////////////////////////////////
int FileDelete(char *fname)
{
  int status = 0;

#ifdef PARALLEL
  MPI_File fh;

  status = MPI_File_open(MPI_COMM_WORLD, fname,
                         MPI_MODE_CREATE | MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fh);

  if (status == MPI_SUCCESS) {
    status = MPI_File_set_size(fh, 0);   // collective
    MPI_File_close(&fh);
  }

  MPI_Barrier(MPI_COMM_WORLD);

#else
  status = std::remove(fname);
#endif

  return status;
}

////////////////////////////////////////////////////////////////////////
/*!
  Write an array to file using raw binary format.
  No AL distributed array descriptor is used here.
  
  \param [in]   buf      address of buffer to be written
  \param [in]   offset   file offset in bytes [Only for parallel writing]
  \param [in]   cnt      the number of element to be written
  \param [in]   size     the sizeof() datatype
  \param [in]   fname    the file name
*/
////////////////////////////////////////////////////////////////////////
void FileWriteArray(void *buf, long int offset, long int cnt, size_t size,
                    char *fname)
{
  char *cbuf = (char *) buf;
  cnt *= size;  // Convert number of elements to byte 

#ifdef PARALLEL
  MPI_File fhw;
  MPI_Status status;
  MPI_Offset start;   // Store end of file position 
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_open(MPI_COMM_WORLD, fname, 
                MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fhw);
  MPI_File_get_position(fhw, &start);
   
  MPI_File_set_view (fhw, start+offset, MPI_BYTE, MPI_CHAR, "native", MPI_INFO_NULL);
  MPI_File_write_all(fhw, cbuf, cnt, MPI_CHAR, &status);
  MPI_File_close(&fhw);
#else
  FILE *fp = fopen (fname,"ab"); // Open for appending in binary mode
  fwrite (cbuf, sizeof(char), cnt, fp);
  fclose (fp);    
#endif
  
}

////////////////////////////////////////////////////////////////////////
/*!
  Write or Append an ASCII header to a file.
  
  \param [in]   buffer  the input buffer (a string). 
  \param [in]   fname   the file name
  \param [in]   mode    an integer specifying if the file has to be
  opened for the first time (mode = -1) or
  appended (mode >= 0)
*/
////////////////////////////////////////////////////////////////////////
void FileWriteHeader(char *buffer, char fname[], int mode)
{
  int nelem = strlen(buffer);
#ifdef PARALLEL
  MPI_File fhw;
  MPI_Status status;
  MPI_Offset start; // Used to store end of file position 

  if (mode < 0){
  
  // -- Open file and delete it -- 
  
    MPI_File_open(MPI_COMM_WORLD, fname,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                  MPI_INFO_NULL, &fhw);
    MPI_File_close(&fhw);
    
  // -- Open file for writing --     

    MPI_File_open(MPI_COMM_WORLD, fname, 
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
    start = 0;
  }else{
    MPI_File_open(MPI_COMM_WORLD, fname, 
                  MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fhw);
    MPI_File_get_position(fhw, &start);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_set_view(fhw, start, MPI_BYTE, MPI_CHAR, "native", MPI_INFO_NULL);
  if( g_rank == 0 ){
    MPI_File_write(fhw, buffer, nelem, MPI_CHAR, &status); 
  }
  MPI_File_close(&fhw);
  MPI_Barrier(MPI_COMM_WORLD);
#else
  FILE *fp;
  if (mode < 0){
    fp = fopen (fname,"wb"); // -- Open for the first time -- 
  }else{
    fp = fopen(fname, "ab"); // -- for appending in binary mode -- 
  }
  // print(fname);
  fwrite (buffer, sizeof(char), nelem, fp);
  fclose (fp);    
#endif
}
