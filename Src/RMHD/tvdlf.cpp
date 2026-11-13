////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Rusanov Lax-Friedrichs flux for RMHD Eqns

  Solve riemann problem for the relativistic MHD equations
  using th Lax-Friedrichs Rusanov solver.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \authors A. Mignone  (andrea.mignone@unito.it)  \n
           M. Rossazza (marco.rossazza@unito.it)
           M. Bugli    (matteo.bugli@cea.fr)
  \date    Jun 09, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir, bool fofc>
static void LF_SolverImpl(Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void LF_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) LF_SolverImpl<IDIR, false>(data, Dts, grid, box);
  else if (dir == JDIR) LF_SolverImpl<JDIR, false>(data, Dts, grid, box);
  else if (dir == KDIR) LF_SolverImpl<KDIR, false>(data, Dts, grid, box);
}

#if PCP_FOFC == YES
// Reuse the ordinary LF kernel with marked faces and cell-centered states.
// Dissipation and CFL use the same local characteristic-speed bound.
void PCP_LF_Solver(Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) LF_SolverImpl<IDIR, true>(data, Dts, grid, box);
  else if (dir == JDIR) LF_SolverImpl<JDIR, true>(data, Dts, grid, box);
  else if (dir == KDIR) LF_SolverImpl<KDIR, true>(data, Dts, grid, box);
}
#endif

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir, bool fofc>
static void LF_SolverImpl(Data *data, timeStep *Dts, Grid *grid, RBox *box)
 {
  int      i, j, k, nv;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int      ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  RunConfig &run_config = *(data->run_config);
  #ifdef GLM_MHD
  double glm_ch = run_config.glm_ch;
  double ch2    = glm_ch*glm_ch;
  #endif
  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  // ------------------------------------------------------
  // 0. Compatibility checks
  // ------------------------------------------------------

  #if (defined STAGGERED_MHD) 
  #if (CT_EMF_AVERAGE == UCT_GFORCE) || (CT_EMF_AVERAGE == UCT_HLLD)
  QuitPLUTO ("from LF_Solver(): not compatible with UCT_GFORCE or UCT_HLLD \n", 1);
  #endif
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value


  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts)\
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    // 1a. Compute offsets and assign local pointers

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D PSIn(&data->sweep->PSIn[offset],ni);
    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary1D cs2(&data->sweep->cs2[offset],ni);

    Ary1D press(&data->sweep->press[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    // 2. Loop over normal direction

    int ti, tj, tk;
    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl, a2L;
    double vr[NVAR], ur[NVAR], fr[NVAR], Sr, pr, a2R;
    State stateL, stateR;

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    if constexpr (fofc) {
      Ary1Duint32_t mask(&data->sweep->flag[offset],ni);
      if (!((mask(i-1) | mask(i)) & FLAG_PCP_FOFC)) continue;
      Ary2D vc(&data->sweep->vC[offset1],ni,NVAR);
      for (nv = 0; nv < NVAR; nv++) {
        vl[nv] = vc(i-1,nv);
        vr[nv] = vc(i,nv);
      }
#ifdef STAGGERED_MHD
      vl[BXn] = vr[BXn] = Bn(i);
#endif
    } else {
      for (nv = 0; nv < NVAR; nv++) {
        vl[nv] = vL(i,nv);
        vr[nv] = vR(i,nv);
      }
    }

    data->err(ti,tj,tk) = 0;   // reset error 

    #if PHYSICS == GRMHD
    Metric G = DEFAULT_METRIC;
    Metric Gr = DEFAULT_METRIC;
    if (!SingularMetric(x1[ti],x2[tj],x3[tk])){
      ComputeLocalMetric(x1[ti], x2[tj], x3[tk], &G, run_config);
      ComputeLocalMetric(x1[ti+1], x2[tj], x3[tk], &Gr, run_config);
      /* ComputeSqrtDet(&G); */
    }
    #endif

    // 2a. GLM pre-conditioner

    #ifdef GLM_MHD
    GLM_Solve (vl, vr, run_config, dir);
    #endif
    #if PHYSICS != GRMHD
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);
    #else
    PrimToCons (vl, ul, run_config, &G);
    PrimToCons (vr, ur, run_config, &G);
    #endif

    // 2b. Compute sound speed & signal speed

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    a2L = stateL.a2;
    a2R = stateR.a2;
    // Use the same local characteristic bound for LLF dissipation and CFL.
    // FOFC uses cell-centered states; the ordinary solver uses face states.
    double sl_min, sl_max, sr_min, sr_max;

    #if PHYSICS != GRMHD
    MaxSignalSpeed<dir> (stateL, sl_min, sl_max);
    MaxSignalSpeed<dir> (stateR, sr_min, sr_max);
    #else
    MaxSignalSpeed<dir> (stateL, sl_min, sl_max, &G);
    MaxSignalSpeed<dir> (stateR, sr_min, sr_max, &G);
    #endif

    double cL = MAX(fabs(sl_max), fabs(sl_min));
    double cR = MAX(fabs(sr_max), fabs(sr_min));

    cmax(i) = MAX(cL, cR);
    cs2(i)  = 0.5*(a2L + a2R);
    Sl = -cmax(i);
    Sr =  cmax(i);

    // 2c. Compute Flux & Powell's interface state

    #if PHYSICS != GRMHD
    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);
    #else
    Flux<dir> (stateL, fl, pl, run_config, &G);
    Flux<dir> (stateR, fr, pr, run_config, &G);
    #endif

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) {
      flux(i,nv) = 0.5*(fl[nv] + fr[nv]) - 0.5*cmax(i)*(ur[nv] - ul[nv]);
    }
    press(i) = 0.5*(pl + pr);

    #if (DIVB_CONTROL == EIGHT_WAVES)
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    // 2d. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2e. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2f. Store emf-flux components for CT 

    #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
    #if (CT_EMF_AVERAGE == ARITHMETIC) || (CT_EMF_AVERAGE == CT_CONTACT)
    CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == CT_FLUX
    CT_Store_Flux<dir>(data, fl, fr, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_HLL
    #if PHYSICS != GRMHD
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk);
    #else
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk, &G);
    #endif
    #endif
    #endif

  }}}  // End for (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("LF_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from LF_Solver()",err_kernel);
  #endif

}
