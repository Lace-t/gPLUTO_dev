////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief 3D wrapper for conservative/primitive conversion.

  Provide 3D wrappers to the standard 1D conversion functions
  ConsToPrim() and PrimToCons().

  \authors A. Mignone  (andrea.mignone@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it)
  \date    Feb 24, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Convert V (primitive) to U (conservative).
  
  \param [in]  data        pointer to data structure
  \param [in]  iarr        integer specifying which array is to be converted:
  0 means Vc --> Uc; 1 means Vpc --> Upc
  \param [in]  grid        pointer to Grid structure
  \param [in]  box         pointer to RBox structure
*/
////////////////////////////////////////////////////////////////////////
void PrimToCons3D(const Data *data, int iarr, Grid *grid, RBox *box)
{
  int i,j,k, nv;
  double v[NVAR], u[NVAR];

  const Ary4D& Vc = (iarr == 0 ? data->Vc: data->Vpc);  
  const Ary4D& Uc = (iarr == 0 ? data->Uc: data->Upc);  
  RunConfig &run_config = *(data->run_config);
  #pragma omp target enter data map(to: box[:1], Uc, Vc,run_config)

  mynvtxstart_("PrimToCons3D",RED);

  #pragma omp target teams distribute parallel for collapse(3) private(v, u) //num_teams(1240) thread_limit(128) 
  #pragma acc parallel loop collapse(3) present(Uc,Vc) \
              private(v[:NVAR], u[:NVAR]) copyin(box)
  BOX_LOOP(*box, i,j,k){

    #pragma acc loop
    for (nv = 0; nv < NVAR; nv++)  v[nv] = Vc(i,j,k,nv);

    PrimToCons (v, u, run_config);

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    OrthoToCons(u, metric);
    #endif
   
    #pragma acc loop
    for (nv = 0; nv < NVAR; nv++)  Uc(i,j,k,nv) = u[nv];
  }

  #if GEOMETRY == CARTESIAN
  #endif

  mynvtxstop_();
  #pragma omp target exit data map(release: box[:1], Uc, Vc, run_config)
}

////////////////////////////////////////////////////////////////////////
/*!
 * Convert conservative array U to primitive V on the given *box.
 *
  \param [in]  data        pointer to data structure
  \param [in]  iarr        integer selecting the array to be converted:
 *                          0 means Uc; 1 means Upc
  \param [in]  grid        pointer to Grid structure
  \param [in]  box         pointer to RBox structure
*/ ///////////////////////////////////////////////////////////////////// */
void ConsToPrim3D(Data *data, int iarr, Grid *grid, RBox *box)
{
  int i,j,k, nv;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag  
  double v[NVAR];
  double u[NVAR];
  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];
  #ifdef HIGH_ORDER
  const Ary4D& Vc      = (iarr == 0 ? data->Vc: data->Vpc);  
  const Ary4D& Uc      = (iarr == 0 ? data->Uc: data->Upc);  
  #else
  const Ary4D& Vc      = data->Vc;  
  const Ary4D& Uc      = data->Uc;  
  #endif

  RunConfig &run_config = *(data->run_config);

  mynvtxstart_("ConsToPrim3D",RED);

  // ------------------------------------------------------
  // 1. Loop over box and convert to primitive 
  // ------------------------------------------------------

  //#pragma omp target data map(tofrom: err_kernel)
  {

  #pragma acc parallel loop collapse(3) present (data, grid, run_config, box) \
          private(v[:NVAR], u[:NVAR]) copy(err_kernel) 
  #pragma omp target teams distribute parallel for collapse(3) private(v,u)
  BOX_LOOP(*box, i,j,k){
    State state;
    STATE_SET(state, v, u, NULL, x1[i], x2[j], x3[k]);
    state.flag = &(data->flag(i,j,k));

    // 1a. Make a local copy, U() -> u,  
    //     check for NaNs

    data->err(i,j,k) = 0;
    #pragma acc loop seq
    for (nv = 0; nv < NVAR; nv++) {
      u[nv] = Uc(i,j,k,nv);
      if (isnan(u[nv])) {  
        ATOMIC_ERR_UPDATE (err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE (data->err(i,j,k), ERR_NAN, ERR_QUIT);
      }
    }

    // 1b. Convert conservative to primitive

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    ConsToOrtho (u, metric); 
    #endif

    uint32_t local_err = ConsToPrim (state, run_config);

    #if GEOMETRY != CARTESIAN
    OrthoToCons (u, metric);  
    #endif

    // 1c. Fallback convex blending

    #if FALLBACK_MODE == FALLBACK_BLENDING
    if (data->integration_mode == FALLBACK_BLENDING){
      if (local_err){
        printLog ("Error found in FALLBACK mode, zone [%d, %d, %d]\n",i,j,k);
        ErrorShow(local_err); 
        QuitPLUTO("",1);
      }
    }
    if (data->integration_mode == DEFAULT){
      if (local_err){
        double th = 0.0; 
        if (g_intStage == 1) th = 1.0;
        if (g_intStage == 2) th = 0.5;
        if (g_intStage == 3) th = 1.0;
        NVAR_LOOP(nv) {
          u[nv] = data->Uc0(i,j,k,nv)*(1.0 - th) + th*data->Usafe(i,j,k,nv);
        }
        printLog("Mappers3D(): Replacing zone [%d, %d, %d] in ConsToPrim\n",
                  i,j,k);
        local_err = ConsToPrim (state, run_config);
        if (local_err) {
          printLog ("Error: conversion did not succeed using failsafe values\n");    
          ErrorShow(local_err); 
          QuitPLUTO("",1);
        }
      }
    }
    #endif

    // 1d. Update error flag after conversion

    ATOMIC_ERR_UPDATE(err_kernel, local_err);
    KERNEL_ERR_UPDATE(data->err(i,j,k), local_err);
    if (local_err & ERR_CONS2PRIM_FAIL) {
      ATOMIC_ERR_UPDATE(err_kernel,  ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_QUIT);
    }

    // 1e. Copy back v -> V(). Also, make sure 
    //     any change done on u is transferred to U()

    if (local_err){ 
      #pragma acc loop seq
      for (nv = 0; nv < NVAR; nv++)  Uc(i,j,k,nv) = u[nv];
    }else{
      #if ENTROPY_SWITCH != NO
      Uc(i,j,k,ENG)  = u[ENG];
      Uc(i,j,k,ENTR) = u[ENTR];
      #endif
    }

    #pragma acc loop seq
    for (nv = 0; nv < NVAR; nv++)  Vc(i,j,k,nv) = v[nv];
  }
}

  // ------------------------------------------------------
  // 2. Error Check condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("ConsToPrim3D()", data, err_kernel, box, grid);
  #endif

  // 3. Quit check: by default, quit if ERR_QUIT is raised
  //    unless FALLBACK is activated and we're currently 
  //    running using the DEFAULT method.

#if WARNING_MESSAGES >= 0 
  int quit_condition = (err_kernel & ERR_QUIT);

  #if FALLBACK_MODE != NO
  quit_condition *= (data->integration_mode != DEFAULT);
  Fallback_Trigger(data->err, data->flag, box);
  #endif

  if (quit_condition) QuitPLUTO("from ConsToPrim3D()",err_kernel);

    // ----------------------------------------------------------------
    // IDEALLY ONE WOULD LIKE TO CALL WriteData from here, producing one file
    // for diagnostic:
    // 
    // data.<step>.<stage>.err.vtk
    // 
    // This data contains both err and flag within a "single file"
    // This output does not use output->nfile.
    // WriteData (output[ERR_OUTPUT], )
    // 
    // To this end:
    // 
    // - Add one more output format called "err.vtk"  (VTK_ERR_OUTPUT)
    // - make sure err(i,j,k) and flag(i,j,k) reside in data->()
    // - Add a pointer to Runtime inside Data struct
    // - Adapt WriteVTK() (to use unsigned int instead of float)
    // - Modify WriteData accordingly
    // - Call 
    // WriteData (data, data->Dts, data->runtime->output + VTK_ERR_OUTPUT, grid);
    // --------------------------------------------------------------------------
  
#endif

  mynvtxstop_();
}
