////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLL Riemann solver for MHD.

  Solve the Riemann problem for the MHD equations using the
  single-state HLL solver by Toro.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b References:
   - "Riemann Solver and Numerical Methods for Fluid Dynamics"
      by E.F. Toro (Chapter 10)
   - "On Godunov-Type Method near Low Densities"
     by B. Einfeldt, C.D. Munz, P.L. Roe, JCP 92, 273-295 (1991)

  \authors A. Mignone  (andrea.mignone@unito.it)  \n
           M. Rossazza (marco.rossazza@unito.it)  \n
           A. Suriano  (alessio.suriano@unito.it) \n
           C. Zanni    (zanni@oato.inaf.it)       \n 
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void HLL_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void HLL_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) HLL_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) HLL_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) HLL_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLL_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
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
  QuitPLUTO ("! HLL_Solver(): not compatible with UCT_GFORCE or UCT_HLLD \n", 1);
  #endif
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value


  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data, run_config, Dts, grid)\
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  
    // Compute offsets and assign local pointers

    long int offset  = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D PSIn(&data->sweep->PSIn[offset],ni);
    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary1D cs2(&data->sweep->cs2[offset],ni);
    Ary1D press(&data->sweep->press[offset],ni);

    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);
    Ary2D vL (&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR (&data->sweep->vR[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl;
    double ur[NVAR], vr[NVAR], fr[NVAR], Sr, pr;
    double bgf[3];

    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR)	{ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR)	{ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR)	{ti = j; tj = k; tk = i;}

    #pragma acc loop seq
    NFLUID_LOOP(nv) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    STATE_SET(stateL, vl, ul, bgf, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, bgf, x1[ti], x2[tj], x3[tk]);

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. GLM pre-conditioner

    #ifdef GLM_MHD
    GLM_Solve (vl, vr, run_config, dir);
    #endif
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    // 2b. Compute Background field, sound velocity, HLL speeds
    //     and max propagation speed for time step computation.

    #if BACKGROUND_FIELD == YES
    BackgroundField (stateL, run_config); // stateL->bgf = stateR->bgf
    #endif

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);

    cmax(i) = MAX(fabs(Sl), fabs(Sr));    
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    // 2c. Compute HLL flux

    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    double scrh = 1.0 / (bR - bL);
    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) {
      flux(i,nv)  = bL*bR*(ur[nv] - ul[nv]) + bR*fl[nv] - bL*fr[nv];
      flux(i,nv) *= scrh;
    }
    press(i) = (bR*pl - bL*pr)*scrh;

    // 2d. Compute Powell's / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    // 2e. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2f. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2g. Store emf-flux components for CT

    #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
    #if (CT_EMF_AVERAGE == ARITHMETIC) || (CT_EMF_AVERAGE == CT_CONTACT)
    CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == CT_FLUX
    CT_Store_Flux<dir>(data, fl, fr, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_HLL
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk);
    #endif
    #endif

    // 2h. Conservative Energy correction step [Experimental]

    #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
    double EgL, EgR;
    double fEgL = HD_EnergyFlux(vl, &EgL, dir, run_config);
    double fEgR = HD_EnergyFlux(vr, &EgR, dir, run_config);
    flux(i,ENG) = (bR*fEgL - bL*fEgR + bL*bR*(EgR - EgL))*scrh; 
    #elif CT_ENERGY_CORRECTION == 5
    flux(i,ENG) = (bR*fEgL - bL*fEgR)*scrh; 
    #endif

  }}} // End Domain loop

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("HLL_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from HLL_Solver()",err_kernel);
  #endif
}
