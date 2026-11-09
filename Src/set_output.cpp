////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Set/retrieve output data attributes.

  The function SetOutput() sets, for each output data type (DBL, FLT,
  VTK etc..), the default attributes of the corresponding ::Output structure.
  These include the variable name, a pointer to the actual
  \c Ary4D or \c Ary3D class, the centering of the variable (center/staggered), 
  a conditional inclusion flag (telling if the corresponding variable has
  to be written in the specified format), and so on.

  Variables names are stored in <tt>output->var_name[nv]</tt> in 
  sequential order with the following convention:

      nv = {0 ... NVAR-1}  {NVAR ... NVAR + ns*DIM-1} {NVAR + ns*DIM + Vuser...}
           Zone-centered       Staggered Fields           User-defined
                ^                      ^                        ^               
  Vary4D           *Vary3D[nv-NVAR]           *Vary3D[nv-NVAR]         

  \c Vary4D points to different \c Ary4D depending on the integration method:
 
      output->Vary4D = *(data->Uc);  (High-Order)
      output->Vary4D = *(data->Vpc); (2nd order)

  The function SetOutputVar() can be used to include or exclude a given
  variable to be written using a particular output format.

  The function LinkUserVar() should be used from the ComputeUserVar() 
  function to link a user-defined \c Ary3D into a specified output.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    May 28, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifndef VTK_VECTOR_DUMP
 #define VTK_VECTOR_DUMP NO
#endif

static Output *g_output; // Make it visibile to whole file

////////////////////////////////////////////////////////////////////////
/*!
  Set default attributes (variable names, pointers to data structures,
  filename extensions, etc...) of the output structures.
  
  \param [in] d        pointer to Data structure
  \param [in] runtime  pointer to Runtime structure
  \param [in] runtime  pointer to Grid structure
*/
////////////////////////////////////////////////////////////////////////
void SetOutput (Data *data, Runtime *runtime, Grid *grid)
{
  int nv, type;
  int *ntot = grid->ntot;
  int *nint = grid->nint;
  static double ***Vc_buf_AL, ***Vs_buf_AL[3];
  static double *Vc_buf, *Vs_buf[3];
  Output *output;

  g_output = output = runtime->output;

  // 0. Allocate buffer memory for I/O operations 
  //    (either serial or parallel)

  if (Vc_buf_AL == NULL){

    Vc_buf = ARRAY_1D(nint[KDIR]*nint[JDIR]*nint[IDIR], double); 
    Vc_buf_AL = ARRAY_3D(ntot[KDIR], ntot[JDIR], ntot[IDIR], double); 
  
    DIM_EXPAND(
      Vs_buf_AL[0] = ARRAY_3D(ntot[KDIR],   ntot[JDIR],   ntot[IDIR]+1, double); ,
      Vs_buf_AL[1] = ARRAY_3D(ntot[KDIR],   ntot[JDIR]+1, ntot[IDIR],   double); ,
      Vs_buf_AL[2] = ARRAY_3D(ntot[KDIR]+1, ntot[JDIR],   ntot[IDIR],   double); 
    )

    DIM_EXPAND(
      Vs_buf[0] = ARRAY_1D(nint[KDIR]*nint[JDIR]*(nint[IDIR]+1), double); ,
      Vs_buf[1] = ARRAY_1D(nint[KDIR]*(nint[JDIR]+1)*nint[IDIR], double); ,
      Vs_buf[2] = ARRAY_1D((nint[KDIR]+1)*nint[JDIR]*nint[IDIR], double); 
    )
  }

  // 1. Set output file extensions

  strcpy(output[DBL_OUTPUT].ext, "dbl");
  strcpy(output[FLT_OUTPUT].ext, "flt");
  strcpy(output[DBL_H5_OUTPUT].ext, "dbl.h5");
  strcpy(output[FLT_H5_OUTPUT].ext, "flt.h5");
  strcpy(output[VTK_OUTPUT].ext, "vtk");
  strcpy(output[VTK_ERR_OUTPUT].ext, "vtk_err");
  strcpy(output[TAB_OUTPUT].ext, "tab");
  strcpy(output[PPM_OUTPUT].ext, "ppm");
  strcpy(output[PNG_OUTPUT].ext, "png");

  // 2. Loop on output types.

  for (type = 0; type < MAX_OUTPUT_TYPES; type++){

    output = runtime->output + type;
    output->var_name = ARRAY_2D(MAX_OUTPUT_VARS, 32, char);
    strcpy(output->dir, runtime->output_dir); // output directory is the same     
                                              // for all outputs (easy to change) 
    output->nfile    = -1;

    // 2a. Exclude particle outputs

    #if PARTICLES != NO
    if (output->type == PARTICLES_DBL_OUTPUT ||
        output->type == PARTICLES_FLT_OUTPUT ||
        output->type == PARTICLES_VTK_OUTPUT ||
        output->type == PARTICLES_TAB_OUTPUT ||
        output->type == PARTICLES_HDF5_OUTPUT) continue;
    #endif

    // 2b. All outputs share the same memory buffer

    output->Vc_buf_AL = Vc_buf_AL;  
    output->Vs_buf_AL[0] = Vs_buf_AL[0];  
    output->Vs_buf_AL[1] = Vs_buf_AL[1];  
    output->Vs_buf_AL[2] = Vs_buf_AL[2];  

    output->Vc_buf    = Vc_buf;  
    output->Vs_buf[0] = Vs_buf[0];  
    output->Vs_buf[1] = Vs_buf[1];  
    output->Vs_buf[2] = Vs_buf[2];  

    // 2c. Set variable names for cell center data 
    //     with 0 <= nv < NVAR. 
    //     Also, set output Ary pointers.

    SetDefaultVarNames(output);

    NVAR_LOOP(nv){
      #ifdef HIGH_ORDER
      if (type == DBL_OUTPUT) output->Vary4D = &(data->Uc);
      else                    output->Vary4D = &(data->Vpc);
      #else
      output->Vary4D = &(data->Vc);
      #if CONS_DBL == YES
      output->Vary4D = &(data->Uc);
      #endif
      #endif

      output->dump_var[nv] = YES;
      output->stag_var[nv] = -1; // -- means cell centered -- 
    }

    // 2d. Exclude some varibales from being dumped to disk (like entropy)
    //     For GLM: never dump scalars unless output is used for restart.

    #if ENTROPY_SWITCH
    output->dump_var[ENTR] = NO; // Exclude entropy from being written
    #endif

    #ifdef GLM_MHD
     if (output->type == DBL_OUTPUT || output->type == DBL_H5_OUTPUT) {
       output->dump_var[PSI_GLM] = YES;
       #ifdef PHI_GLM
       output->dump_var[PHI_GLM] = YES;
       #endif
     } else {
       output->dump_var[PSI_GLM] = NO;
       #ifdef PHI_GLM
       output->dump_var[PHI_GLM] = NO;
       #endif
     }
    #endif

    // 2e. Add staggered fields.
    //      Note: exclude them if the output type is
    //      not DBL_OUTPUT or DBL_H5_OUTPUT.

    nv = NVAR;
    #ifdef STAGGERED_MHD
    int dir;
    size_t svar_name = sizeof(output->var_name[nv]);
    DIM_LOOP(dir){
      snprintf (output->var_name[nv], svar_name, "Bx%01ds",1+dir);

      #ifdef HIGH_ORDER
      if (type == DBL_OUTPUT) output->Vary3D[nv-NVAR] = &(data->Vs[BX1s+dir]);
      else                    output->Vary3D[nv-NVAR] = &(data->Vps[BX1s+dir]); 
      #else
      output->Vary3D[nv-NVAR] = &(data->Vs[BX1s+dir]);
      #endif

      bool flag = (type == DBL_OUTPUT || type == DBL_H5_OUTPUT);
      output->dump_var[nv] = flag;
      output->stag_var[nv] = dir;         
      nv++;
    }

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_LOOP(dir){
      snprintf (output->var_name[nv], svar_name, "Ex%01ds",1+dir);

      #ifdef HIGH_ORDER
      if (type == DBL_OUTPUT) output->Vary3D[nv-NVAR] = &(data->Vs[EX1s+dir]);
      else                    output->Vary3D[nv-NVAR] = &(data->Vps[EX1s+dir]); 
      #else
      output->Vary3D[nv-NVAR] = &(data->Vs[EX1s+dir]);
      #endif

      bool flag = (type == DBL_OUTPUT || type == DBL_H5_OUTPUT);
      output->dump_var[nv] = flag;
      output->stag_var[nv] = dir;         
      nv++;
    }
    #endif
    #endif  // STAGGERED_MHD

    // 2f. Add vector potential [TODO]

    #if UPDATE_VECTOR_POTENTIAL == YES
     #if DIMENSIONS == 3
      strcpy(output->var_name[nv], "Ax1");
      output->stag_var[nv++] = -1;  // -- vector potential is
                                          // computed at cell center -- 

      strcpy(output->var_name[nv], "Ax2");
      output->stag_var[nv++] = -1;
     #endif
     strcpy(output->var_name[nv], "Ax3");
     output->stag_var[nv++] = -1;
    #endif
    output->nvar = nv;

    // 2g. Add user-defined variables, for nv > NVAR, to each output type.
    //     User-def variables will be written only if the corresponding 
    //     output->dump_var[nv] is set to YES from UserDefOutput().  
    //     User-def variables are Assumed to be cell centered

    for (int k = 0; k < runtime->user_var; k++){
      strcpy (output->var_name[nv + k], runtime->user_var_name[k]);
      output->stag_var[nv + k] = -1; // -- assume cell-centered -- 
      output->dump_var[nv + k] = NO;
    }

    // 2h. Update total number of variables

    output->user_var = runtime->user_var; 
    output->nvar    += runtime->user_var;

  } // end loop on type

  // 3. Defaults: dump density only in ppm and png formats 

  SetOutputVar ("rho", PPM_OUTPUT, YES);
  SetOutputVar ("rho", PNG_OUTPUT, YES);
}

////////////////////////////////////////////////////////////////////////
/*!
  Include ('flag == YES') or exclude ('flag == NO') the variable
  corresponding to 'var_name' in or from the output type 'type'.
  If 'type' corresponds to an image (ppm or png), create a
  correspdonding Image structure.
  
  \param [in] var_name   the name of the variable (e.g. "rho", "vx1",...)
  \param [in] type       select the output type (e.g., DBL_OUTPUT,
  VTK_OUTPUT, and so forth)
  \param [in] flag       an integer values (YES/NO).
*/
////////////////////////////////////////////////////////////////////////
int SetOutputVar (const char *var_name, int type, int flag)
{
  int nv;
  Output *output;

  output = &(g_output[type]);

  for (nv = output->nvar; nv--; ) {
    if (strcmp(output->var_name[nv], var_name) == 0) {
      output->dump_var[nv] = flag;
      return 0;
    }
  }

  print ("! SetOutputVar(): var_name '%s' cannot be set/unset for writing\n",
             var_name);
  return(1);
}

////////////////////////////////////////////////////////////////////////
/*!
  Find the name index of the user-supplied Ary \c A, and create 
  a pointer to an \c Ary3D in the specified output type.
  
  \param[in] A         pointer to user-defined Ary3D
  \param[in] var_name  the name of the variable corresponding to \c A
  \param[in] type      the output type
*/
////////////////////////////////////////////////////////////////////////
void LinkUserVar(Ary3D *A, const char *var_name, int type)
{
  int nv = 0;

  if (A->isAllocated() == NO){
    QuitPLUTO("! LinkUserVar(): uservar '%s' is not allocated\n", 1);
  }

  // Search variable name index 

  while (strcmp(g_output->var_name[nv], var_name)){
    nv++;
    if (nv == MAX_OUTPUT_VARS) return; 
  }

  g_output[type].Vary3D[nv-NVAR] = A;
  g_output[type].dump_var[nv]    = YES;

}
