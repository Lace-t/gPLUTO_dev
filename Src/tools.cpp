////////////////////////////////////////////////////////////////////////
/* Collect various tools.

   This files collect various tools such as memory allocation, wrapper
   to prim/cons variables and so forth.

   This code is distributed under the BSD 3-Clause License (see below).
   Last modified: 14 May 2020 by A. Mignone                           */
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
unsigned int BitPosition( uint32_t x )
/*!
  Count push off bits to right until 0
  101 => 10 => 1 => 0
  which means hibit was 3rd bit, its value is 2^3
  
  Bit counter starts at 0                                            */
////////////////////////////////////////////////////////////////////////
{
  unsigned int log2Val = 0 ;

  while( x>>=1 ) log2Val++;  // div by 2 until find log2. 
  return log2Val;
}


////////////////////////////////////////////////////////////////////////
void BinaryPrint(uint32_t n)
/*! 
  Print integer in binary format, eg..
  
  BinaryPrint(5)  -->  "0000000101"
*/
////////////////////////////////////////////////////////////////////////
{
  int l = sizeof(n) * 8;
  for (int i = l - 1 ; i >= 0; i--) {
    print ("%x", (n & (1 << i)) >> i);
  }
  print ("  [Read rgt->lft]\n");
}

#ifdef PARALLEL
////////////////////////////////////////////////////////////////////////
void GetNeighbourRanks (const Grid *grid, int **nranks, int ***nranks_edge, 
                        int ***nranks_corner)
/*! 
  Find the ranks of neighbour processors direction by direction.
  Store them in the array \c nrank[dir][s] where \c dir is the
  direction and <tt> s = 0,1 </tt> stands for the left (0) or right
  (1) processor.
  If a processor touches the physical boundary that is not periodic,
  the corresponding neighbour rank will be set to -1.
  
  \param [in]   grid    pointer to an array of grid structures.
  \param [out]  nranks  an array of integers containing the ranks of
  the neighbouring processors
*/
////////////////////////////////////////////////////////////////////////
{
  int dir, dir1, dir2, coords[3], rnk;
  MPI_Comm cartcomm;
  
#if USE_AL == NO
cartcomm = new_SZ.cart_comm;
#else
  AL_Get_cart_comm(SZ, &cartcomm);
#endif

  // ------------------------------------------------------------
  // Neighbour exists when there's no physical boundary, or
  // when PERIODIC or SHEARING boundary conditions are imposed
  // along that direction.
  // ------------------------------------------------------------

  for (dir = 0; dir < 3; dir++){

  // -- Obtain local processor coordinates -- 

    coords[IDIR] = grid->rank_coord[IDIR];
    coords[JDIR] = grid->rank_coord[JDIR];
    coords[KDIR] = grid->rank_coord[KDIR];

    nranks[dir][0] = nranks[dir][1] = MPI_PROC_NULL;

    if (PROC_EXIST(grid->lbound[dir], grid->nproc[dir])){
      coords[dir] = grid->rank_coord[dir] - 1;
      MPI_Cart_rank(cartcomm, coords, &rnk);
      nranks[dir][0] = rnk;
    }

    if (PROC_EXIST(grid->rbound[dir], grid->nproc[dir])){
      coords[dir] = grid->rank_coord[dir] + 1;
      MPI_Cart_rank(cartcomm, coords, &rnk);
      nranks[dir][1] = rnk;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // ------------------------------------------------------------
  // EDGES
  // ------------------------------------------------------------

  for (dir = 0; dir < 3; dir++){
    dir1 = (dir+1)%3;

  // -- Obtain local processor coordinates -- 

    DIM_EXPAND(coords[IDIR] = grid->rank_coord[IDIR];  ,
               coords[JDIR] = grid->rank_coord[JDIR];  ,
               coords[KDIR] = grid->rank_coord[KDIR];)

    nranks_edge[dir][0][0] = nranks_edge[dir][0][1] = nranks_edge[dir][1][0] = nranks_edge[dir][1][1] = MPI_PROC_NULL;

    if (PROC_EXIST(grid->lbound[dir], grid->nproc[dir])){
      coords[dir] = grid->rank_coord[dir] - 1;
      if (PROC_EXIST(grid->lbound[dir1], grid->nproc[dir1])){
        coords[dir1] = grid->rank_coord[dir1] - 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_edge[dir][0][0] = rnk;
      }
      if (PROC_EXIST(grid->rbound[dir1], grid->nproc[dir1])){
        coords[dir1] = grid->rank_coord[dir1] + 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_edge[dir][0][1] = rnk;
      }
    }

    if (PROC_EXIST(grid->rbound[dir], grid->nproc[dir])){
      coords[dir] = grid->rank_coord[dir] + 1;
      dir1 = (dir+1)%3;
      if (PROC_EXIST(grid->lbound[dir1], grid->nproc[dir1])){
        coords[dir1] = grid->rank_coord[dir1] - 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_edge[dir][1][0] = rnk;
      }
      if (PROC_EXIST(grid->rbound[dir1], grid->nproc[dir1])){
        coords[dir1] = grid->rank_coord[dir1] + 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_edge[dir][1][1] = rnk;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  // ------------------------------------------------------------
  // CORNERS
  // ------------------------------------------------------------

  dir  = 0;
  dir1 = 1;
  dir2 = 2;

  // -- Obtain local processor coordinates -- 

  DIM_EXPAND(coords[IDIR] = grid->rank_coord[IDIR];  ,
             coords[JDIR] = grid->rank_coord[JDIR];  ,
             coords[KDIR] = grid->rank_coord[KDIR];)

  nranks_corner[0][0][0] = nranks_corner[0][0][1] = MPI_PROC_NULL;
  nranks_corner[0][1][0] = nranks_corner[0][1][1] = MPI_PROC_NULL;
  nranks_corner[1][0][0] = nranks_corner[1][0][1] = MPI_PROC_NULL;
  nranks_corner[1][1][0] = nranks_corner[1][1][1] = MPI_PROC_NULL;

  if (PROC_EXIST(grid->lbound[dir], grid->nproc[dir])){ // 0 
    coords[dir] = grid->rank_coord[dir] - 1;
    if (PROC_EXIST(grid->lbound[dir1], grid->nproc[dir1])){ // 00 
      coords[dir1] = grid->rank_coord[dir1] - 1;
      if (PROC_EXIST(grid->lbound[dir2], grid->nproc[dir2])){ // 000 
        coords[dir2] = grid->rank_coord[dir2] - 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[0][0][0] = rnk;
      }
      if (PROC_EXIST(grid->rbound[dir2], grid->nproc[dir2])){ // 001 
        coords[dir2] = grid->rank_coord[dir2] + 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[0][0][1] = rnk;
      }
    }
    if (PROC_EXIST(grid->rbound[dir1], grid->nproc[dir1])){ // 01 
      coords[dir1] = grid->rank_coord[dir1] + 1;
      if (PROC_EXIST(grid->lbound[dir2], grid->nproc[dir2])){ // 010 
        coords[dir2] = grid->rank_coord[dir2] - 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[0][1][0] = rnk;
      }
      if (PROC_EXIST(grid->rbound[dir2], grid->nproc[dir2])){ // 011 
        coords[dir2] = grid->rank_coord[dir2] + 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[0][1][1] = rnk;
      }
    }
  }

  if (PROC_EXIST(grid->rbound[dir], grid->nproc[dir])){ // 1 
    coords[dir] = grid->rank_coord[dir] + 1;
    if (PROC_EXIST(grid->lbound[dir1], grid->nproc[dir1])){ // 10 
      coords[dir1] = grid->rank_coord[dir1] - 1;
      if (PROC_EXIST(grid->lbound[dir2], grid->nproc[dir2])){ // 100 
        coords[dir2] = grid->rank_coord[dir2] - 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[1][0][0] = rnk;
      }
      if (PROC_EXIST(grid->rbound[dir2], grid->nproc[dir2])){ // 101 
        coords[dir2] = grid->rank_coord[dir2] + 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[1][0][1] = rnk;
      }
    }
    if (PROC_EXIST(grid->rbound[dir1], grid->nproc[dir1])){ // 11 
      coords[dir1] = grid->rank_coord[dir1] + 1;
      if (PROC_EXIST(grid->lbound[dir2], grid->nproc[dir2])){ // 110 
        coords[dir2] = grid->rank_coord[dir2] - 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[1][1][0] = rnk;
      }
      if (PROC_EXIST(grid->rbound[dir2], grid->nproc[dir2])){ // 111 
        coords[dir2] = grid->rank_coord[dir2] + 1;
        MPI_Cart_rank(cartcomm, coords, &rnk);
        nranks_corner[1][1][1] = rnk;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  return;
}
#endif  // PARALLEL 

#ifdef GNUPLOT_HEADER
////////////////////////////////////////////////////////////////////////
void GnuplotSetting(Runtime *runtime, Grid *grid)
/*!
 * Set-up a gnuplot script containing grid info, variable names, etc...
*/
////////////////////////////////////////////////////////////////////////
#define GP_VAR_NAME_SIZE 32
{
  int  k, nv;
  int  ibeg = grid->lbeg[IDIR];
  int  jbeg = grid->lbeg[JDIR];
  int  kbeg = grid->lbeg[KDIR];
  int  dbl_out = 0;
  int  flt_out = 0;
  int  tab_out = 0;
  // char *var_name[NVAR];
  Output *output;
  time_t  now;
  FILE *fp;

  if (g_rank != 0) return;

  // 0. Allocate memory for an array of strings containing
  //    the names of the variables.
  //    Get the variable names being written to disk.

  // for (nv = 0; nv < NVAR; nv++) var_name[nv] = ARRAY_1D(GP_VAR_NAME_SIZE,char);

  // 1. Determine which output datatypes have been enabled

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){
    output = runtime->output + k;
    if (output->type == DBL_OUTPUT && (output->dt > 0 || output->dn > 0)){
      dbl_out = 1;
    }
    if (output->type == FLT_OUTPUT && (output->dt > 0 || output->dn > 0)){
      flt_out = 1;
    }
    if (output->type == TAB_OUTPUT && (output->dt > 0 || output->dn > 0)){
      tab_out = 1;
    }
  }

  // 2. Open file and write relevant info.
  
  time(&now);
  fp = fopen("pluto.gp","w");

  fprintf (fp, "# ***********************************************************\n");
  fprintf (fp, "#\n");
  fprintf (fp, "# Initialize simulation quantities for reading data   \n");
  fprintf (fp, "# with Gnuplot. \n");
  fprintf (fp, "# Generated by PLUTO on  %s",asctime(localtime(&now)));
  fprintf (fp, "#\n");
  fprintf (fp, "# ***********************************************************\n");

  // 2a. Select default data type and index name

  if (dbl_out){
    fprintf (fp, "if (!exists('dtype')) {\n");
    fprintf (fp,"  dtype = 'dbl' # Default data type\n");
    fprintf (fp,"}\n");
  }else if (flt_out){
    fprintf (fp, "if (!exists('dtype')) {\n");
    fprintf (fp,"  dtype = 'flt' # Default data type\n");
    fprintf (fp,"}\n");
  }
  fprintf (fp,"print '> dtype set to ',dtype\n");
 
  // 2b. Set grid size

  fprintf (fp, "\n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Grid setting                            \n");
  fprintf (fp, "# --------------------------------------- \n");

  fprintf (fp, "nx1 = %d\n", grid->nint_glob[IDIR]);
  fprintf (fp, "nx2 = %d\n", grid->nint_glob[JDIR]);
  fprintf (fp, "nx3 = %d\n", grid->nint_glob[KDIR]);
  
  fprintf (fp, "x1beg = %f; x1end = %f\n", grid->xbeg_glob[IDIR],
                                           grid->xend_glob[IDIR]);

  fprintf (fp, "x2beg = %f; x2end = %f\n", grid->xbeg_glob[JDIR],
                                           grid->xend_glob[JDIR]);

  fprintf (fp, "x3beg = %f; x3end = %f\n", grid->xbeg_glob[KDIR],
                                           grid->xend_glob[KDIR]);
  fprintf (fp, "dx1   = %f\n", grid->dx[IDIR][ibeg]);
  fprintf (fp, "dx2   = %f\n", grid->dx[JDIR][jbeg]);
  fprintf (fp, "dx3   = %f\n", grid->dx[KDIR][kbeg]);

  fprintf (fp, "\n");

  // 2c. Variable names are written in different ways depending on the 
  //     selected output.

  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Variable indices                        \n");
  fprintf (fp, "# --------------------------------------- \n");

  // 2c1. float output

  if (flt_out){
    char all_var[512] = {"> Var Names: "};

    fprintf (fp, "if (dtype eq 'flt') {\n");

    output = &(runtime->output[FLT_OUTPUT]);
    for (nv = 0; nv < output->nvar; nv++){
      if (output->dump_var[nv]){
        fprintf (fp, "  %s = %d\n",output->var_name[nv], nv);
        strcat (all_var, output->var_name[nv]);
        strcat (all_var, "  ");
      }
    }
    fprintf (fp, "  print '%s'\n",all_var);
    fprintf (fp, "}\n");
  }

  // 2c2. double output

  if (dbl_out){
    char all_var[512] = {"> Var Names: "};

    fprintf (fp, "if (dtype eq 'dbl') {\n");

    output = &(runtime->output[DBL_OUTPUT]);
    for (nv = 0; nv < output->nvar; nv++){
      if (output->dump_var[nv]){
        fprintf (fp, "  %s = %d\n",output->var_name[nv], nv);
        strcat (all_var, output->var_name[nv]);
        strcat (all_var, "  ");
      }
    }
    fprintf (fp, "  print '%s'\n",all_var);
    fprintf (fp, "}\n");
  }

  // 2c3. tab output 

  if (tab_out){
    char all_var[512] = {"> Var Names: "};

    fprintf (fp, "if (dtype eq 'tab') {\n");

    output = &(runtime->output[FLT_OUTPUT]);
    for (nv = 0; nv < output->nvar; nv++){
      if (output->dump_var[nv]){
        // Variables start from the 3rd column
        fprintf (fp, "  %s = %d\n",output->var_name[nv], nv+3); 
        strcat (all_var, output->var_name[nv]);
        strcat (all_var, "  ");
      }
    }
    fprintf (fp, "  print '%s'\n",all_var);
    fprintf (fp, "}\n");
  }

  fprintf (fp, "if (!exists('nvar')) {\n");
  fprintf (fp, "  nvar = 0\n");
  fprintf (fp, "  print '> nvar = ',nvar\n");
  fprintf (fp, "}\n");

  // 3. Make gnuplot count the number of lines in "flt.out"
  //    and "dbl.out"

  fprintf (fp, "\n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Output datatype                         \n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp,"#set term push       # remember current terminal\n");
  fprintf (fp,"#set term unknown\n");

  if (dbl_out){
    fprintf (fp,"#plot 'dbl.out' u 1\n");
    fprintf (fp,"#ndbl = int(GPVAL_DATA_Y_MAX)\n");
  }
  if (flt_out){
    fprintf (fp,"#plot 'flt.out' u 1\n");
    fprintf (fp,"#nflt = int(GPVAL_DATA_Y_MAX)\n");
  }
  fprintf (fp,"#set term pop         # restore terminal\n");

  // 4. Load some other gnuplot scripts

  fprintf (fp, "\n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Load gnuplot scripts                    \n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "load 'macros.gp'\n");
  fprintf (fp, "if (nx2 > 1) load 'pm3d_setting.gp'\n");

  fclose(fp);
}
#endif  // GNUPLOT_HEADER

////////////////////////////////////////////////////////////////////////
int IsLittleEndian (void) 
/*!
  Return 1 if the current architecture has little endian order
*/
////////////////////////////////////////////////////////////////////////
{
  int TestEndian = 1;
  return *(char*)&TestEndian;
}

////////////////////////////////////////////////////////////////////////
void PrintColumnLegend(const char *legend[], int nfields, FILE *fp)
/*!
  Print nicely formatted legend to an ASCII datfile header.
  
  \param [in] legend    an array of strings describing the fields
  \param [in] nfields   the number of fields (also the dimensions of
  legend)
  \param [in] fp        a valid file descriptor.
*/
////////////////////////////////////////////////////////////////////////
{
  int n;
  char num_string[16];

  for (n = 0; n < nfields; n++){
    snprintf(num_string, sizeof(num_string), "[%d]", n);

    if (n == 0) fprintf (fp,"# %-13s ",num_string);
    else        fprintf (fp,"%-13s  ",num_string);
  }
  fprintf (fp,"\n");

  for (n = 0; n < nfields; n++){
    if (n == 0) fprintf (fp,"# %-13s ",legend[n]);
    else        fprintf (fp,"%-13s  ",legend[n]);
  }
  fprintf (fp,"\n");
}

// helper function that converts an integer to a binary string
static string binary_string_int(int value) {
    std::string out;
    unsigned int mask = 1U << (sizeof(int) * 8 - 1); // start with MSB
    bool started = false;

    for (; mask; mask >>= 1) {
        if (value & mask) {
            out.push_back('1');
            started = true;
        } else if (started) {
            out.push_back('0');
        }
    }
    if (!started) out = "0"; // handle value = 0
    return out;
}

////////////////////////////////////////////////////////////////////////
void QuitPLUTO(string err_string, int e_code)
/*!
  Quit the code by printing an error message and exit code in binary
*/
////////////////////////////////////////////////////////////////////////
{
  string bin = binary_string_int(e_code);
  printLog("! QuitPLUTO(): %s [exit code = %s]\n", err_string.c_str(), bin.c_str());

#ifdef PARALLEL
  LogFileFlush();
  MPI_Abort(MPI_COMM_WORLD, e_code);
  MPI_Finalize();        
  exit(e_code);
#else
  exit(e_code);
#endif
}


////////////////////////////////////////////////////////////////////////
void SwapEndian (void *x, const int nbytes) 
/*!
  Swap the byte order of x.
  
  \param [in] x      pointer to the variable being swapped
  \param [in] nbytes data type size
  \return This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
{
  int k;
  static char Swapped[16];
  char *c;

  c = (char *) x;

  for (k = nbytes; k--; ) Swapped[k] = *(c + nbytes - 1 - k);
  for (k = nbytes; k--; ) c[k] = Swapped[k];
}
////////////////////////////////////////////////////////////////////////
int StringArrayIndex (char *str_arr[], int size, char *str)
/*!
  Find the index of an array of strings whose value matches the
  input string \c *str.
  
  \param [in]  str_arr  the array of strings
  \param [in]  size     the size of the array
  \param [in]  str      the string to search for
*/
////////////////////////////////////////////////////////////////////////
{
  int i;

  for (i = 0; i < size; i++){
    if (!strcmp(str_arr[i],str)) return i;
  }
  return -1;
}