////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Write fluid variables to disk.

  WriteData() is the main driver for writing data arrays in any of
  the available formats (binary, VTK, HDF5, etc...).

  - For .dbl, .flt or .vtk file formats, access to binary files is
    provided by the functions in bin_io.cpp
  - HDF5 files are handled by hdf5_io.cpp
  - image files are handled by write_img.cpp
  - tabulated ascii files are handled by write_tab.cpp

  This function also updates the corresponding <tt>.out</tt> file associated
  with the output data format.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    May 28, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Write data to disk using any of the available formats.
  
  \param [in] d      pointer to PLUTO Data structre
  \param [in] Dts    pointer to timeStep structure
  \param [in] output the output structure corresponding to a given
  format
  \param [in] grid   pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void WriteData (Data *data, timeStep *Dts, Output *output, Grid *grid)
{
  int    nv;
  char   filename[512], sline[512];
  int    SZ_stag[3] = {SZ_stagx, SZ_stagy, SZ_stagz};
  int    *ntot = grid->ntot;
  int    single_file;
  char *dummy;  // dummy variable to suppress unused variable warning
  size_t dsize;
  double units[MAX_OUTPUT_VARS];
  float ***Vpt3;
  void *Vpt;
  // #ifdef FARGO
  // static double ***vphi_bck;
  // double **wA = FARGO_Velocity();
  // #endif
  FILE *fout, *fbin;
  #ifdef PARALLEL
  time_t tbeg, tend;
  #endif
  RunConfig &run_config = *(data->run_config);

  RBox tot_box, box;

  RBoxDefine (0, ntot[IDIR]-1,
              0, ntot[JDIR]-1,
              0, ntot[KDIR]-1, &tot_box);

  RBoxDefine (grid->lbeg, grid->lend, &box);
  #pragma acc enter data copyin (box, tot_box)
  
  // 0. Provides high-order punctual-value output variables
   
#ifdef HIGH_ORDER
  // BoundaryConservative (data, grid);
  // g_intStage = -1;

  // // On 1st step, ho limiter has not been computed 
  // // yet. Alternative: move it in Startup() ?
 
  // #if (HO_LAP_LIMITER != NO) && (GEOMETRY == CARTESIAN)  
  // if (g_stepNumber == 0) DiscontinuityDetect(data, grid, &box);   
  // #endif

  // PointValue (data, grid, &box);
  // Boundary (data, 1, grid);

#endif  // #ifdef HIGH_ORDER 

  // 1. Increment the file number and initialize units

  output->nfile++;
  print ("> Writing file #%d (%s) to disk...\n", output->nfile, output->ext);

  #ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (g_rank == 0) time(&tbeg);
  #endif

  for (nv = 0; nv < MAX_OUTPUT_VARS; nv++) units[nv] = 1.0;
  if (output->cgs) GetCGSUnits(units, run_config);

  // 2. User-def output call

  UserDefOutput (data, grid);

  // 3. With FARGO, we can output total or residual velocity

  #if (defined FARGO) &&  (FARGO_OUTPUT_VTOT == YES)
    // if (vphi_bck == NULL) vphi_bck = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    // BOX_LOOP(tot_box, i,j,k){
    //   vphi_bck[k][j][i] = data->Vc(i,j,k,VX1+SDIR); /* Backup residual */
    // }
    FARGO_ComputeTotalVelocity (data, VX1+SDIR, grid);
  #endif

  // 4. Select the output type

  if (output->type == DBL_OUTPUT) {

    // 4a. DBL Output
    //     - \b DBL_OUTPUT:
    //       Double-precision data files can be written using single or
    //       multiple file mode. 
    //     - single file, serial: we open the file just once before
    //       the main variable loop, dump variables and then close.
    //     - single file, parallel: the distributed array descriptor sz is
    //       different for cell-centered or staggered data type and we
    //       thus have to open and close the file after each variable
    //       has been dumped.
    //     - multiple files: we open, write to and close the
    //       file one each loop cycle.
    //
    //     \note In all cases, the pointer to the data array that has to be 
    //           written must be cast into (void *) and the starting index of 
    //           the array must be zero.
  
    int sz;
    single_file = strcmp(output->mode,"single_file") == 0;
    dsize = sizeof(double);

    if (single_file){  // Single output file
    
      // Build the filename 
      snprintf(filename, sizeof(filename), "%s/data.%04d.%s",
              output->dir, output->nfile, output->ext);
      FileDelete (filename);  // Avoid partial fill of pre-existing files

#if USE_AL == NO
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        // Building the output buffer with only active points
        FillOutputBuffer(nv, output, &box);

        int stag_dir = output->stag_var[nv];
        if (stag_dir == -1) {        // cell-centered 
          Vpt = (void *)(output->Vc_buf);
        } else {                     // staggered     
          Vpt = (void *)(output->Vs_buf[stag_dir]);
        }

        FileWriteData(Vpt, dsize, FILE_DBL, stag_dir, grid, filename, 0);
      }
#else
      #ifdef PARALLEL
      long long offset = 0;
      #else
      fbin = FileOpen (filename, 0, "w");
      #endif

      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        FillOutputBuffer_AL (nv, output, &tot_box);
        if  (output->stag_var[nv] == -1) {  // Cell-centered data
          sz  = SZ;
          Vpt = (void *)(output->Vc_buf_AL[0][0]);
        }else{
          int stag_dir = output->stag_var[nv];
          sz  = SZ_stag[stag_dir];
          Vpt = (void *)(output->Vs_buf_AL[stag_dir][0][0]);
        }

        #ifdef PARALLEL
        fbin = FileOpen (filename, sz, "w");
        AL_Set_offset(sz, offset);
        #endif
        FileWriteData_AL (Vpt, dsize, sz, fbin, output->stag_var[nv], grid);
        #ifdef PARALLEL
        offset = AL_Get_offset(sz);
        FileClose(fbin, sz);
        #endif
      }
      #ifndef PARALLEL
      FileClose(fbin, sz);
      #endif
#endif // USE_AL

    }else{          //  Multiple files

#if USE_AL == NO
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        FillOutputBuffer(nv, output, &box);

        snprintf(filename, sizeof(filename), "%s/%s.%04d.%s",
                output->dir, output->var_name[nv], output->nfile, output->ext);

        FileDelete (filename); // Erase file before writing

        int stag_dir = output->stag_var[nv];
        if (stag_dir == -1) {        // cell-centered 
          Vpt = (void *)(output->Vc_buf);
        } else {                     // staggered     
          Vpt = (void *)(output->Vs_buf[stag_dir]);
        }

        FileWriteData(Vpt, dsize, FILE_DBL, 
                          stag_dir, grid, filename, 0);

      }
#else
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        FillOutputBuffer_AL (nv, output, &tot_box);

        snprintf (filename, sizeof(filename), "%s/%s.%04d.%s", 
                  output->dir, output->var_name[nv], output->nfile, output->ext);

        FileDelete (filename);  // Avoid partial fill of pre-existing files

        if  (output->stag_var[nv] == -1) {  // Cell-centered data
          sz = SZ;
          Vpt = (void *)(output->Vc_buf_AL[0][0]);
        }else{
          int stag_dir = output->stag_var[nv];
          sz  = SZ_stag[stag_dir];
          Vpt = (void *)(output->Vs_buf_AL[stag_dir][0][0]);
        }

        fbin = FileOpen (filename, sz, "w");
        FileWriteData_AL (Vpt, dsize, sz, fbin, output->stag_var[nv], grid);
        FileClose (fbin, sz);
      }
#endif // USE_AL
    }

    // With FARGO, write2D array for orbital vel. to separate file.   
    // Needed also for restarting purposes.

    #ifdef FARGO
    FARGO_Write(data, output->dir, output->nfile, grid);
    #endif

  } else if (output->type == FLT_OUTPUT) {

    // 4b. FLT output for cell-centered data

    single_file = strcmp(output->mode,"single_file") == 0;

    if (single_file){  // Single output file 

      snprintf (filename, sizeof(filename), "%s/data.%04d.%s", 
                output->dir, output->nfile,output->ext);
      FileDelete (filename);  // Avoid partial fill of pre-existing files 

#if USE_AL == NO
      for (nv = 0; nv < output->nvar; nv++) {
        int stag_dir = output->stag_var[nv];

        if (output->dump_var[nv] == NO) continue;

        FillOutputBuffer (nv, output, &box);

        float *Vpt_flt = Convert_dbl2flt(output->Vc_buf, units[nv],0, grid);
        FileWriteData((void *)Vpt_flt, sizeof(float), FILE_FLT, 
                      -1, grid, filename, 0);
      }
#else
      fbin = FileOpen (filename, SZ_float, "w");
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        FillOutputBuffer_AL (nv, output, &tot_box);
        Vpt3 = Convert_dbl2flt_AL(output->Vc_buf_AL, units[nv],0, grid);
        Vpt = (void *)Vpt3[0][0];
        FileWriteData_AL (Vpt, sizeof(float), SZ_float, fbin,
                          output->stag_var[nv], grid);
      }
      FileClose(fbin, SZ_float);
#endif

    }else{            // Multiple files 

#if USE_AL == NO
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        snprintf (filename, sizeof(filename), "%s/%s.%04d.%s", 
                  output->dir, output->var_name[nv],output->nfile, output->ext);

        FileDelete (filename);  // Avoid partial fill of pre-existing files 

        FillOutputBuffer (nv, output, &box);
        float *Vpt_flt = Convert_dbl2flt(output->Vc_buf, units[nv],0, grid);
        FileWriteData ((void *)Vpt_flt, sizeof(float), FILE_FLT, 
                        -1, grid, filename, 0);
      }
#else
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        snprintf (filename, sizeof(filename), "%s/%s.%04d.%s", 
                  output->dir, output->var_name[nv],output->nfile, output->ext);

        FileDelete (filename);  // Avoid partial fill of pre-existing files 
        fbin = FileOpen (filename, SZ_float, "w");

        FillOutputBuffer_AL (nv, output, &tot_box);
        Vpt3 = Convert_dbl2flt_AL(output->Vc_buf_AL, units[nv],0, grid);
        Vpt = (void *)Vpt3[0][0];
        FileWriteData_AL (Vpt, sizeof(float), SZ_float, fbin,
                          output->stag_var[nv], grid);
        FileClose (fbin, SZ_float);
      }
#endif
    }

  }else if (output->type == DBL_H5_OUTPUT || output->type == FLT_H5_OUTPUT){

    // 4c.  HDF5 (static grid) output (single/double precision)

    #ifdef USE_HDF5
    single_file = YES;
    WriteHDF5 (output, Dts, grid);
    #else
    print ("! WriteData: HDF5 library not available\n");
    // #pragma acc exit data delete (box, tot_box)
    // return;
    #endif

  }else if (output->type == VTK_OUTPUT || output->type == VTK_ERR_OUTPUT) {

    // 4d. VTK Output
    //     in order to enable parallel writing, files must be closed and
    //     opened again for scalars, since the distributed array descriptors 
    //     used by ArrayLib (Float_Vect) and (float) are different. This is
    //     done using the AL_Get_offset() and AL_Set_offset() functions.    
    //
    //     Note for VTK_ERR_OUTPUT: the filename is established at the time 
    //     of the call using the step number and the stage numbers, e.g., 
    //     data.<step>.<stage>.err.vtk

    single_file = strcmp(output->mode,"single_file") == 0;
    snprintf (filename, sizeof(filename), "%s/data.%04d.%s", 
              output->dir, output->nfile, output->ext);

    if (single_file){  // Single output file

#if USE_AL == NO

      FileDelete (filename); // Erase file before writing

      if (g_rank == 0) WriteVTK_Header(filename, Dts, grid);

      for (nv = 0; nv < output->nvar; nv++) { // Write scalars
        if (output->dump_var[nv] == NO) continue;
        FillOutputBuffer (nv, output, &box);
        WriteVTK_Scalar (filename, output->Vc_buf, units[nv], 
                         output->var_name[nv], grid);
      }

#else

      if (g_rank == 0) WriteVTK_Header(filename, Dts, grid);
      
      for (nv = 0; nv < output->nvar; nv++) { // Write scalars
        if (output->dump_var[nv] == NO) continue;
        FillOutputBuffer_AL (nv, output, &tot_box);
        WriteVTK_Scalar_AL (filename, output->Vc_buf_AL, units[nv], 
                         output->var_name[nv], grid);
      }

#endif // USE_AL

    }else{          //    Multiple output files

#if USE_AL == NO
      for (nv = 0; nv < output->nvar; nv++) {

        if (output->dump_var[nv] == NO) continue;

        snprintf(filename, sizeof(filename), "%s/%s.%04d.%s",
                output->dir, output->var_name[nv], output->nfile, output->ext);

        FileDelete (filename); // Erase file before writing
        if (g_rank == 0) WriteVTK_Header(filename, Dts, grid);

        // Fill buffer and write scalar */
        FillOutputBuffer(nv, output, &box);
        WriteVTK_Scalar(filename, output->Vc_buf, units[nv],
                            output->var_name[nv], grid);
      }
#else
      for (nv = 0; nv < output->nvar; nv++) {  // Write scalars
        if (output->dump_var[nv] == NO) continue;

        snprintf (filename, sizeof(filename),  "%s/%s.%04d.%s", 
                  output->dir, output->var_name[nv],
                  output->nfile,  output->ext);
        FileDelete (filename);  // Avoid partial fill of pre-existing files 
        if (g_rank == 0) WriteVTK_Header(filename, Dts, grid);

        FillOutputBuffer_AL (nv, output, &tot_box);
        WriteVTK_Scalar_AL (filename, output->Vc_buf_AL, 
                            units[nv], output->var_name[nv], grid);
      }
#endif
    }

  }else if (output->type == TAB_OUTPUT) {

    // 4e. Tabulated (ASCII) output

    single_file = YES;
    snprintf(filename, sizeof(filename), "%s/data.%04d.%s", 
             output->dir, output->nfile, output->ext);

    WriteTabArray (output, filename, grid);

  }else if (output->type == PPM_OUTPUT) {

#if 00
    // ------------------------------------------------------
    // PPM output
    // ------------------------------------------------------

    single_file = NO;
    for (nv = 0; nv < output->nvar; nv++) {
      if (!output->dump_var[nv]) continue;
      sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv],
                                          output->nfile, output->ext);
      WritePPM (output->V[nv], output->var_name[nv], filename, grid);
    }
#endif
  }else if (output->type == PNG_OUTPUT) {
#if 00

    // ------------------------------------------------------
    // PNG  output
    // ------------------------------------------------------

    #ifdef USE_PNG
    single_file = NO;
    for (nv = 0; nv < output->nvar; nv++) {
      if (!output->dump_var[nv]) continue;
      sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv],
                                          output->nfile, output->ext);
      WritePNG (output->V[nv], output->var_name[nv], filename, grid);
    }
    #else
    print ("! PNG library not available\n");
    // #pragma acc exit data delete (box, tot_box)
    // return;
    #endif
#endif
  }

  // 5. Update corresponding ".out" file

  snprintf (filename,sizeof(filename), "%s/%s.out",output->dir, output->ext);

  if (g_rank == 0) {
    if (output->nfile == 0) {
      fout = fopen (filename, "w");
    }else {
      fout = fopen (filename, "r+");
      for (nv = 0; nv < output->nfile; nv++) dummy = fgets (sline, 512, fout);
      fseek (fout, ftell(fout), SEEK_SET);
    }

    // Write a multi-column file 

    fprintf (fout, "%d %12.6e %12.6e %d ",
             output->nfile, Dts->time, Dts->dt, g_stepNumber);

    if (single_file) fprintf (fout,"single_file ");
    else             fprintf (fout,"multiple_files ");

    if (IsLittleEndian()) fprintf (fout, "little ");
    else                  fprintf (fout, "big ");

    for (nv = 0; nv < output->nvar; nv++) {
      if (output->dump_var[nv]) {
        fprintf (fout, "%s ", output->var_name[nv]);
      }
    }
    fprintf (fout,"\n");
    fclose (fout);
  }

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (g_rank == 0){
    time(&tend);
    print ("  [%5.2f sec ]\n",difftime(tend,tbeg));
  }
#endif

  #pragma acc exit data delete (box, tot_box)

}

////////////////////////////////////////////////////////////////////////
/*!
  Compute an array of c.g.s units
  
  \param [in]  u   an array containing the c.g.s units of the primitive
  variables, e.g., <tt> u[RHO] </tt> will be in gr/cm^3,
  <tt> u[VX1] </tt> will be in cm/s, etc...
*/
////////////////////////////////////////////////////////////////////////
void GetCGSUnits (double *u, RunConfig &run_config)
{
  double unit_density  = run_config.unit_density;
  double unit_velocity = run_config.unit_velocity;

  #if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
  unit_velocity = CONST_c;
  #endif

  u[RHO] = unit_density;
  u[VX1] = unit_velocity;
  u[VX2] = unit_velocity;
  u[VX3] = unit_velocity;
  #if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
  double unit_mag = unit_velocity*sqrt(4.0*CONST_PI*unit_density);
  u[BX1] = unit_mag;
  u[BX2] = unit_mag;
  u[BX3] = unit_mag;
  #endif
  #if HAVE_ENERGY
  double unit_pressure = run_config.unit_density*unit_velocity*unit_velocity;
  u[PRS] = unit_pressure;
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Copy the "nv" variable array from the correpsonding Ary4D / Ary3D 
  to a standard array buffer before writing to disk, that is,
  
  V[i][j][k] <-- Ary(i,j,k)
  
  \param [in]     nv      the output variable to be copied
  \param [in,out] output  pointer to output structure
  \param [in]     box     the copy region
*/
////////////////////////////////////////////////////////////////////////
void FillOutputBuffer (int nv, Output *output, RBox *box)
{
  int ind,i,j,k;
  double *Vc_buf, *Vs_buf[3];

  Vc_buf = output->Vc_buf;
  Vs_buf[0] = output->Vs_buf[0];
  Vs_buf[1] = output->Vs_buf[1];
  Vs_buf[2] = output->Vs_buf[2];

  RBox boxs;

  if (nv < NVAR) {  // These variables are all cell-centered by default
    Ary4D *Vary4D = output->Vary4D;
    BOX_LOOP(*box, i,j,k) {
      ind = BOX_INDEX(*box,i,j,k);
      Vc_buf[ind] = (*Vary4D)(i,j,k,nv);
    }
  }else{
    if (output->stag_var[nv] == -1){   // Cell centered, but user-def
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      BOX_LOOP(*box, i,j,k) {
        ind = BOX_INDEX(*box,i,j,k);
        Vc_buf[ind] = (*Vary3D)(i,j,k);
      }
    }else{      
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      int is = (nv-NVAR)%3;
      is = output->stag_var[nv];

      RBoxDefine(box->ibeg, box->iend+(is==IDIR), 
                 box->jbeg, box->jend+(is==JDIR), 
                 box->kbeg, box->kend+(is==KDIR), &boxs);

      BOX_LOOP(boxs, i,j,k) {
        ind = BOX_INDEX(boxs,i,j,k);
        Vs_buf[is][ind] = (*Vary3D)(i,j,k);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
/*!
  Copy the "nv" variable array from the correpsonding Ary4D / Ary3D 
  to a standard array buffer before writing to disk, that is,
  
  V[i][j][k] <-- Ary(i,j,k)
  
  \param [in]     nv      the output variable to be copied
  \param [in,out] output  pointer to output structure
  \param [in]     box     the copy region
*/
////////////////////////////////////////////////////////////////////////
void FillOutputBuffer_AL (int nv, Output *output, RBox *box)
{
  int i,j,k;
  double ***Vc_buf, ***Vs_buf[3];

  Vc_buf = output->Vc_buf_AL;
  Vs_buf[0] = output->Vs_buf_AL[0];
  Vs_buf[1] = output->Vs_buf_AL[1];
  Vs_buf[2] = output->Vs_buf_AL[2];

  if (nv < NVAR) {  // These variables are all cell-centered by default
    Ary4D *Vary4D = output->Vary4D;
    BOX_LOOP(*box, i,j,k) Vc_buf[k][j][i] = (*Vary4D)(i,j,k,nv);
  }else{
    if (output->stag_var[nv] == -1){   // Cell centered, but user-def
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      BOX_LOOP(*box, i,j,k) Vc_buf[k][j][i] = (*Vary3D)(i,j,k);
    }else{                             // Staggered fields      
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      int is = (nv-NVAR)%3;
      is = output->stag_var[nv];
      BOX_LOOP(*box, i,j,k) Vs_buf[is][k][j][i] = (*Vary3D)(i,j,k);
    }
  }
}


////////////////////////////////////////////////////////////////////////
/*!
  This function does the opposite job of FillOutputBuffer()
  
  Ary(i,j,k) <-- V[i][j][k]
*/
////////////////////////////////////////////////////////////////////////
void ReadOutputBuffer_AL (int nv, Output *output, RBox *box)
{
  int i,j,k;
  double ***Vc_buf, ***Vs_buf[3];

  Vc_buf = output->Vc_buf_AL;
  Vs_buf[0] = output->Vs_buf_AL[0];
  Vs_buf[1] = output->Vs_buf_AL[1];
  Vs_buf[2] = output->Vs_buf_AL[2];

  if (nv < NVAR) {      // These variables are all cell-centered by default
    Ary4D *Vary4D = output->Vary4D;
    BOX_LOOP(*box, i,j,k) (*Vary4D)(i,j,k,nv) = Vc_buf[k][j][i];
  }else{
    if (output->stag_var[nv] == -1){
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      BOX_LOOP(*box, i,j,k) (*Vary3D)(i,j,k) = Vc_buf[k][j][i];
    }else{
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      int is = (nv-NVAR)%3;
      is = output->stag_var[nv];
      BOX_LOOP(*box, i,j,k) (*Vary3D)(i,j,k) = Vs_buf[is][k][j][i];
    }
  }
}

////////////////////////////////////////////////////////////////////////
/*!
  This function does the opposite job of FillOutputBuffer()
  
  Ary(i,j,k) <-- V[ind]
*/
////////////////////////////////////////////////////////////////////////
void ReadOutputBuffer (int nv, Output *output, RBox *box)
{
  int ind,i,j,k;
  double *Vc_buf, *Vs_buf[3];

  Vc_buf = output->Vc_buf;
  Vs_buf[0] = output->Vs_buf[0];
  Vs_buf[1] = output->Vs_buf[1];
  Vs_buf[2] = output->Vs_buf[2];

  RBox boxs;

  if (nv < NVAR) {  // These variables are all cell-centered by default
    Ary4D *Vary4D = output->Vary4D;
    BOX_LOOP(*box, i,j,k) {
      ind = BOX_INDEX(*box,i,j,k);
      (*Vary4D)(i,j,k,nv) = Vc_buf[ind];
    }
  }else{
    if (output->stag_var[nv] == -1){   // Cell centered, but user-def
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      BOX_LOOP(*box, i,j,k) {
        ind = BOX_INDEX(*box,i,j,k);
        (*Vary3D)(i,j,k) = Vc_buf[ind];
      }
    }else{      
      Ary3D *Vary3D = output->Vary3D[nv-NVAR];
      int is = (nv-NVAR)%3;
      is = output->stag_var[nv];

      RBoxDefine(box->ibeg, box->iend+(is==IDIR), 
                 box->jbeg, box->jend+(is==JDIR), 
                 box->kbeg, box->kend+(is==KDIR), &boxs);

      BOX_LOOP(boxs, i,j,k) {
        ind = BOX_INDEX(boxs,i,j,k);
        (*Vary3D)(i,j,k) = Vs_buf[is][ind];
      }
    }
  }
}

