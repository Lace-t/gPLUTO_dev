////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Write data in VTK format.

  Collection of basic functions to write legacy VTK binary files.
  Files can be written in serial or parallel mode and consists of the
  basic five parts :

  1. File version and identifier
  2. Header, consisting of a string
  3. File format
  4. Dataset structure: describes the geometry and topology of the
     dataset
  5. Dataset attributes. This section is used to write the actual binary
     data as a vector or scalar data

  The mesh topology and the variable datasets are written using single
  precision (4 bytes) binary format.
  VTK file are written following big endian order.
  Therefore, we swap endianity only if local architecture has little
  endian ordering.

  The WriteVTK_Header() function provides the basic functionality for
  steps 1, 2, 3 and 4. Only processor 0 does the actual writing.
  For cartesian/cylindrical geometries the default grid topology is
  <em>rectilinear grid</em> (\c VTK_RECTILINEAR_GRIDS) whereas for
  polar/spherical we employ <em>structured grid</em> (\c VTK_STRUCTURED_GRID)
  to provide a convenient mapping to a cartesian mesh.

  <b> Note for 2D datasets</b>: in order to produce a strictly 2D
  dataset we always set the third coordinate (x3) to zero.
  For this reason, in 2D spherical coordinates we swap the role of the
  "y" and "z" coordinates.

  The WriteVTK_Scalar() is fully parallel and is used to write data with
  the scalar attribute (by default these include density, pressure, tracer
  and user-defined variables).

  \b Reference

  http://www.vtk.org/VTK/img/file-formats.pdf

  \author A. Mignone (andrea.mignone@unito.it)
  \date   May 28, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define VTK_RECTILINEAR_GRID    14
#define VTK_STRUCTURED_GRID     35

#ifndef VTK_FORMAT
  #if GEOMETRY == CARTESIAN || (GEOMETRY == CYLINDRICAL && INCLUDE_JDIR == NO)
    #define VTK_FORMAT  VTK_RECTILINEAR_GRID
  #else
    #define VTK_FORMAT  VTK_STRUCTURED_GRID
  #endif
#endif

#ifndef VTK_TIME_INFO       // Enable this keyword only if you're using 
 #define VTK_TIME_INFO  NO  // VisIt to display data results            
#endif

////////////////////////////////////////////////////////////////////////
void WriteVTK_Header (char *filename, timeStep *Dts, Grid *grid)
/*!
  Write VTK header in parallel or serial mode.
  In parallel mode only processor 0 does the actual writing
  (see al_io.c/AL_Write_header).
  
  \param [in]  filename   vtk filename
  \param [in]  Dts        pointer to timeStep structure (used for time info)
  \param [in]  grid       pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
  long int nx1 = grid->nint_glob[IDIR];
  long int nx2 = grid->nint_glob[JDIR];
  long int nx3 = grid->nint_glob[KDIR];
  float x1, x2, x3;

  int ibeg = grid->lbeg[IDIR];
  int jbeg = grid->lbeg[JDIR];
  int kbeg = grid->lbeg[KDIR];

  FILE *fvtk;

  // 0. Open file for writing (proc #0 only) 

  fvtk = fopen(filename, "wb");
  if (fvtk == NULL) {
    QuitPLUTO ("WriteVTK_Header: invalid pointer",1);
  }

  // 1. Header section starts with 
  //
  //    - file version and identifier
  //    - Header
  //    - file format
  //    - structure (rect./struct.) 

  fprintf (fvtk, "# vtk DataFile Version 2.0\n");
  fprintf (fvtk, "gPLUTO %s VTK Data\n", gPLUTO_VERSION);
  fprintf (fvtk, "BINARY\n");

  #if VTK_FORMAT == VTK_RECTILINEAR_GRID
  fprintf(fvtk, "DATASET %s\n","RECTILINEAR_GRID");
  #elif VTK_FORMAT == VTK_STRUCTURED_GRID
  fprintf(fvtk, "DATASET %s\n","STRUCTURED_GRID");
  #endif

  // 2. Generate time info (VisIt reader only) 

  #if VTK_TIME_INFO == YES
  fprintf (fvtk, "FIELD FieldData 1\n");
  fprintf (fvtk, "TIME 1 1 double\n");
  double tt = Dts->time;
  if (IsLittleEndian()) SWAP_VAR(tt);

  fwrite(&tt, sizeof(double), 1, fvtk);
  fprintf (fvtk, " \n");
  #else
  #endif // VTK_TIME_INFO 

  // --------------------------------------------
  // 3. Dimension & coordinates
  // --------------------------------------------

  fprintf (fvtk, "DIMENSIONS %ld %ld %ld\n",
            nx1 + INCLUDE_IDIR, nx2 + INCLUDE_JDIR, nx3 + INCLUDE_KDIR);

  #if VTK_FORMAT == VTK_RECTILINEAR_GRID

  // 3a1. Allocate memory and define node coordinates only once 

  static float *xnode, *ynode, *znode;
  if (xnode == NULL){
    xnode = ARRAY_1D(nx1 + INCLUDE_IDIR, float);
    ynode = ARRAY_1D(nx2 + INCLUDE_JDIR, float);
    znode = ARRAY_1D(nx3 + INCLUDE_KDIR, float);

    for (int i = 0; i < nx1 + INCLUDE_IDIR; i++){
      x1 = (float)(grid->xL_glob[IDIR][i+ibeg]);
      if (IsLittleEndian()) SWAP_VAR(x1);
      xnode[i] = x1;
    }
    for (int j = 0; j < nx2 + INCLUDE_JDIR; j++){
      x2 = (float)(grid->xL_glob[JDIR][j+jbeg]);
      if (IsLittleEndian()) SWAP_VAR(x2);
      ynode[j] = x2;
    }
    for (int k = 0; k < nx3 + INCLUDE_KDIR; k++){
      x3 = (float)(grid->xL_glob[KDIR][k+kbeg]);
      if (IsLittleEndian()) SWAP_VAR(x3);
      #if !INCLUDE_KDIR
      znode[k] = 0.0;
      #else
      znode[k] = x3;
      #endif
    }
  }

  // 3a2. Write rectilinear grid information

  fprintf (fvtk, "X_COORDINATES %ld float\n", nx1 + INCLUDE_IDIR);
  fwrite(xnode, sizeof(float), nx1+INCLUDE_IDIR, fvtk);

  fprintf (fvtk, "\nY_COORDINATES %ld float\n", nx2 + INCLUDE_JDIR);
  fwrite(ynode, sizeof(float), nx2+INCLUDE_JDIR, fvtk);

  fprintf (fvtk, "\nZ_COORDINATES %ld float\n", nx3 + INCLUDE_KDIR);
  fwrite(znode, sizeof(float), nx3+INCLUDE_KDIR, fvtk);

  #elif VTK_FORMAT == VTK_STRUCTURED_GRID

  // 3b1. Allocate memory and define node_coord

  static float  **node_coord;
 
  if (node_coord == NULL) node_coord = ARRAY_2D(nx1 + INCLUDE_IDIR, 3, float);

  fprintf (fvtk, "POINTS %ld float\n", 
          (nx1+INCLUDE_IDIR)*(nx2+INCLUDE_JDIR)*(nx3+INCLUDE_KDIR));

  // 3b2. Write structured grid information

  x1 = x2 = x3 = 0.0;
  for (int k = 0; k < nx3 + INCLUDE_KDIR; k++){
  for (int j = 0; j < nx2 + INCLUDE_JDIR; j++){
    for (int i = 0; i < nx1 + INCLUDE_IDIR; i++){
      DIM_EXPAND(x1 = grid->xL_glob[IDIR][ibeg + i];  ,
                 x2 = grid->xL_glob[JDIR][jbeg + j];  ,
                 x3 = grid->xL_glob[KDIR][kbeg + k];)

      #if GEOMETRY == CARTESIAN
      node_coord[i][0] = x1;
      node_coord[i][1] = x2;
      node_coord[i][2] = x3;
      #elif GEOMETRY == CYLINDRICAL
      node_coord[i][0] = x1*cos(x2);
      node_coord[i][1] = x1*sin(x2);
      node_coord[i][2] = x3;
      #elif GEOMETRY == SPHERICAL
      #if DIMENSIONS == 2
      node_coord[i][0] = x1*sin(x2);
      node_coord[i][1] = x1*cos(x2);
      node_coord[i][2] = 0.0;
      #elif DIMENSIONS == 3
      node_coord[i][0] = x1*sin(x2)*cos(x3);
      node_coord[i][1] = x1*sin(x2)*sin(x3);
      node_coord[i][2] = x1*cos(x2);
      #endif
      #endif

      if (IsLittleEndian()){
        SWAP_VAR(node_coord[i][0]);
        SWAP_VAR(node_coord[i][1]);
        SWAP_VAR(node_coord[i][2]);
      }
    }
    fwrite (node_coord[0], sizeof(float), 3*(nx1+INCLUDE_IDIR), fvtk);
  }}
  #endif // VTK_FORMAT == VTK_RECTILINEAR_GRID / VTK_STRUCTURED_GRID

  // 4. Dataset attributes [will continue later 
  //    in WriteVTK_Scalar()]

  fprintf (fvtk, "\nCELL_DATA %ld\n", nx1*nx2*nx3);
  fclose(fvtk);
  fflush(fvtk);
}
#undef VTK_STRUCTURED_GRID
#undef VTK_RECTILINEAR_GRID


////////////////////////////////////////////////////////////////////////
/*!
  Write VTK scalar field.
  
  \param [in]   filename   vtk filename
  \param [in]   V          pointer to 3D data array
  \param [in]   unit       the corresponding cgs unit (if specified, 1 otherwise)
  \param [in]   var_name   the variable name appearing in the VTK file
  \param [in]   grid       pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void WriteVTK_Scalar (char *filename, double *V, double unit,
                      char *var_name, Grid *grid)
{
  long int offset = 0;
  float *Vflt;
  FILE *fvtk;

  // 0. Synchronize all ranks

  #ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // 1. Header section done by proc #0 only

  if (g_rank == 0){
    fvtk = fopen(filename,"a");
    fprintf (fvtk, "\nSCALARS %s float\n", var_name);
    fprintf (fvtk, "LOOKUP_TABLE default\n");
    offset = ftell(fvtk);   // byte position right after header 
    fclose(fvtk);
  }

  // 2. Broadcast offset to all ranks

  #ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&offset, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  #endif

  // 3. Convert double -> float

  Vflt = Convert_dbl2flt(V, unit, IsLittleEndian(), grid);

  // 4. Write binary data

  FileWriteData((void *)Vflt, sizeof(float), FILE_FLT, 
                 -1, grid, filename, (FileOffset)offset);

}

////////////////////////////////////////////////////////////////////////
/*!
  Write VTK scalar field.
  
  \param [in]   filename   vtk filename
  \param [in]   V          pointer to 3D data array
  \param [in]   unit       the corresponding cgs unit (if specified, 1 otherwise)
  \param [in]   var_name   the variable name appearing in the VTK file
  \param [in]   grid       pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void WriteVTK_Scalar_AL (char *filename, double ***V, double unit,
                      char *var_name, Grid *grid)
{
#if USE_AL == YES
  long int offset = 0;
  float ***Vflt;
  FILE *fvtk;

  // 0. Synchronize all ranks

  #ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // 1. Header section done by proc #0 only

  if (g_rank == 0){
    fvtk = fopen(filename,"a");
    fprintf (fvtk, "\nSCALARS %s float\n", var_name);
    fprintf (fvtk, "LOOKUP_TABLE default\n");
    offset = ftell(fvtk);  // get size of header
    fclose(fvtk);
  } // end if (g_rank == 0)

  // 2. Broadcast offset to all processors

  #ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
  #endif

  // 3. Re-open file in serial or parallel mode

  fvtk  = FileOpen(filename, SZ_float, "a");
  #ifdef PARALLEL
  AL_Set_offset(SZ_float, offset);
  #endif
  Vflt = Convert_dbl2flt_AL(V, unit, IsLittleEndian(), grid);
  FileWriteData_AL (Vflt[0][0], sizeof(float), SZ_float, fvtk, -1, grid);
  FileClose(fvtk, SZ_float);
#endif
}
