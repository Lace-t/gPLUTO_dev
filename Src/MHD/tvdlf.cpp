////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Lax-Friedrechs (Rusanov) Riemann solver for MHD.

  Solve the Riemann problem for the adiabatic and isothermal MHD
  equations using the Lax-Friedrichs Rusanov Riemann solver with local
  maximum characteristic speed:
  \f[
      \hat{F}_{i-\HALF} = \frac{F^L_{i-\HALF} + F^R_{i-\HALF}}{2}
   - c^{\rm max}_{i+\HALF}\frac{U^R_{i-\HALF} - U^L_{i-\HALF}}{2}
     \qquad{\rm where}\quad c^{\rm max} = \lambda_{\rm max}\left(
          \frac{V^R_{i+\HALF} + V^L_{i+\HALF}}{2}\right)
  \f]
  where \f$\lambda_{\rm \max}\f$ is a function of the arithmetic
  average between the left and the right states.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b References:
    - "Comparison of Some Flux Correced Transport and Total Variation
       Diminishing Numerical Schemes for Hydrodynamics and Magnetohydrodynamics
       Problems", Toth and Odstrcil, JCP (1996), 128,82

  \authors A. Mignone  (andrea.mignone@unito.it) \n
           M. Rossazza (marco.rossazza@unito.it) \n
           C. Zanni    (zanni@oato.inaf.it)      \n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void LF_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void LF_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) LF_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) LF_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) LF_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!

*/ 
////////////////////////////////////////////////////////////////////////
template<int dir>
void LF_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
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
  QuitPLUTO ("! LF_Solver(): not compatible with UCT_GFORCE or UCT_HLLD", 1);
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
  #pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts)\
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
 
    // Compute offsets and assign local pointers

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
  
    // Private variable declaration

    int ti, tj, tk;

    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl;
    double ur[NVAR], vr[NVAR], fr[NVAR], Sr, pr;
    double vLR[NVAR];
    double bgf[3];

    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

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
    
    // 2b. Background field

    #if BACKGROUND_FIELD == YES
    BackgroundField (stateL, run_config); // stateL->bgf = stateR->bgf
    #endif

    // 2c. Get max eigenvalue and sound speed

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) vLR[nv] = 0.5*(vl[nv] + vr[nv]);

    {
      State stateLR;
      double cminLR, cmaxLR;

      STATE_SET(stateLR, vLR, NULL, bgf, x1[ti], x2[tj], x3[tk]);

      SoundSpeed2 (stateLR, run_config);
      stateL.a2 = stateR.a2 = stateLR.a2;
      MaxSignalSpeed<dir> (stateLR, cminLR, cmaxLR);

      cmax(i) = MAX(fabs(cminLR), fabs(cmaxLR));
      cs2(i)  = stateLR.a2;
      Sl = -cmax(i);
      Sr =  cmax(i);
    }

    // 2d. Compute Lax-Friedrichs flux

    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) {
      flux(i,nv) = 0.5*(fl[nv] + fr[nv]) - 0.5*cmax(i)*(ur[nv] - ul[nv]);
    }
    press(i) = 0.5*(pl + pr);

    // 2e. Compute Powell's / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    // 2f. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2g. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2h. Store emf-flux components for CT

    #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
    #if (CT_EMF_AVERAGE == ARITHMETIC) || (CT_EMF_AVERAGE == CT_CONTACT)
    CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == CT_FLUX
    CT_Store_Flux<dir>(data, fl, fr, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_HLL
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk);
    #endif
    #endif

    // 2i. Conservative Energy correction step
    //      [Experimental]

    #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
    double EgL, EgR;
    double fEgL = HD_EnergyFlux(vl, &EgL, dir, run_config);
    double fEgR = HD_EnergyFlux(vr, &EgR, dir, run_config);
    flux(i,ENG) = 0.5*(fEgL + fEgR) - 0.5*cmax(i)*(EgR - EgL); 
    #elif CT_ENERGY_CORRECTION == 4
    flux(i,ENG) = 0.5*(fl[ENG] + fr[ENG]);
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
