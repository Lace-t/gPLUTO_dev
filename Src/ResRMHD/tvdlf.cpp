////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Lax-Friedrechs (Rusanov) Riemann solver for resistive RMHD.

  Solve the Riemann problem for resistive relativistic MHD
  equation using the Lax-Friedrichs Rusanov Riemann solver with local
  maximum characteristic speed equal to the speed of light:
  \f[
      \hat{F}_{i+\HALF} = \frac{F^L_{i+\HALF} + F^R_{i+\HALF}}{2}
   - c^{\rm max}_{i+\HALF}\frac{U^R_{i+\HALF} - U^L_{i-\HALF}}{2}
  \f]

  \authors A. Mignone  (andrea.mignone@unito.it)
           M. Rossazza (marco.rossazza@unito.it)
           A. Suriano  (alessio.suriano@unito.it)
  \date    Jan 06, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void LF_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
void LF_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  if      (dir == IDIR) LF_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) LF_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) LF_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
template<int dir>
void LF_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
/*!
  Solve the Riemann problem between L/R states using a
  Rusanov-Lax Friedrichs flux.
*/
////////////////////////////////////////////////////////////////////////
 {
  int      i, j, k, nv;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int      ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  RunConfig &run_config = *(data->run_config);
  #if DIVB_CONTROL == DIV_CLEANING
  double glm_ch = run_config.glm_ch;
  double ch2    = glm_ch*glm_ch;
  #endif

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  // ------------------------------------------------------
  // 0. Compatibility checks
  // ------------------------------------------------------

  #if   (defined STAGGERED_MHD) \
     && (CT_EMF_AVERAGE != CT_MAXWELL) && (CT_EMF_AVERAGE != ARITHMETIC)
  QuitPLUTO("! LF Solver(): EMF_AVERAGE not supported", 1);
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, EX1, EXn, EXt, EXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------
  
  # pragma omp target teams distribute parallel for collapse(3)
  # pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts)\
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
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl, a2L;
    double ur[NVAR], vr[NVAR], fr[NVAR], Sr, pr, a2R;
    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    for (nv = 0; nv < NVAR; nv++) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. GLM pre-conditioner

    #ifdef GLM_MHD
    GLM_Solve (vl, vr, run_config, dir);
    #endif
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    // 2b. Compute sound speed, cons variables, fluxes,

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);

    a2L = stateL.a2;
    a2R = stateR.a2;

    stateL.eta = stateR.eta = Resistive_eta (stateL, run_config);
    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    cmax(i) = 1.0;
    Sl = -1.0;
    Sr = +1.0;
    {
      double bL = MIN(Sl, 0.0);
      double bR = MAX(0.0, Sr);
      cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);
    }

    // 2c. compute Rusanov flux
       
    for (nv = 0; nv < NFLX; nv++) {
      flux(i,nv) = 0.5*(fl[nv] + fr[nv]) - 0.5*cmax(i)*(ur[nv] - ul[nv]);
    }

    press(i) = 0.5*(pl + pr);

    #if (defined STAGGERED_MHD) && (defined CRG)
    flux(i,CRG) = 0.0;
    #endif

    // 2d. Check for NaNs in flux
 
    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2e. Compute Powell's / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    //  2f. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2g. Store emf-flux components for CT

    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    #if CT_EMF_AVERAGE == ARITHMETIC 
    CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == CT_MAXWELL
    CT_Store_Maxwell<dir>(data, fl, fr, flux, i, ti, tj, tk);
    #endif
    #endif

    // 2h. Conservative Energy correction step
    //       [Experimental]

    #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
    double EgL, EgR;
    double fEgL = HD_EnergyFlux(vl, &EgL, dir, run_config);
    double fEgR = HD_EnergyFlux(vr, &EgR, dir, run_config);
    flux(i,ENG) = 0.5*(fEgL + fEgR) - 0.5*cmax(i)*(EgR - EgL); 
    #elif CT_ENERGY_CORRECTION == 4
    flux(i,ENG) = 0.5*(fl[ENG] + fr[ENG]);
    #endif

  }}} // * End for (i,j,k)

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
