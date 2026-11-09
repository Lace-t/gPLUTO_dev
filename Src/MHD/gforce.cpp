////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief GFORCE Riemann solver for the classical MHD equations.

  Solve the Riemann problem for MHD equations using the
  generalized FORCE flux described in Section 3 of
  Toro \& Titarev, JCP (2006) 216, 403.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b References:
    - "MUSTA Fluxes for systems of conservation laws"
       Toro \& Titarev, JCP (2006) 216, 403

  \authors A. Mignone  (andrea.mignone@unito.it)  \n
           M. Rossazza (marco.rossazza@unito.it)  \n
           S. Truzzi   (stefano.truzzi@unito.it)  \n
           V. Berta    (vittoria.berta@unito.it)  \n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void GFORCE_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void GFORCE_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if     (dir == IDIR) GFORCE_Solver<IDIR>(data, Dts, grid, box);
  else if(dir == JDIR) GFORCE_Solver<JDIR>(data, Dts, grid, box);
  else if(dir == KDIR) GFORCE_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void GFORCE_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
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

  double cfl = RuntimeGet()->cfl;

  // ------------------------------------------------------
  // 0. Compatibility checks
  // ------------------------------------------------------

  #if BACKGROUND_FIELD == YES
  QuitPLUTO ("! GFORCE_Solver(): background field splitting not allowed\n",1);
  #endif

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #if (CT_EMF_AVERAGE == UCT_HLLD)
  QuitPLUTO ("! GFORCE Solver(): not compatible with UCT_HLLD \n", 1);
  #endif
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //      in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data,run_config, grid, Dts)\
                                        copy (err_kernel)
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
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;

    double flux_lf[NVAR];
    double v_lw[NVAR];
    double u_lw[NVAR];
    double f_lw[NVAR];

    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl;
    double ur[NVAR], vr[NVAR], fr[NVAR], Sr, pr;
    double bgf[3];
    double p_lw;

    State stateL, stateR, stateLW;

    // Set state and "true" indices

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    NFLUID_LOOP(nv) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    STATE_SET(stateL,  vl, ul, bgf, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR,  vr, ur, bgf, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateLW, v_lw, u_lw, bgf, x1[ti], x2[tj], x3[tk]);

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. GLM pre-conditioner

    #ifdef GLM_MHD
    GLM_Solve (vl, vr, run_config, dir);
    #endif
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    // 2b. Compute sound & HLL speeds

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    // 2c. Get max propagation speed for dt comp.

    double scrh = MAX(fabs(Sl), fabs(Sr));
    cmax(i) = scrh;
    cs2(i)  = 0.5*(a2L + a2R);

    #ifndef GFORCE_OMEGA
    double om = 1.0/(1.0 + cfl);
    #else
    double om = GFORCE_OMEGA;
    #endif

    // 2d. Compute:
    //      - LF flux
    //      - LW conservative state

    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    double dtdx = 1.0/scrh;
    NFLX_LOOP(nv){
      flux_lf[nv] = 0.5*(fl[nv] + fr[nv]  -    scrh*(ur[nv] - ul[nv]));
      u_lw[nv]    = 0.5*(ul[nv] + ur[nv]) - 0.5*dtdx*(fr[nv] - fl[nv]);
    }

    uint32_t local_flag = 0;  // Use energy during the next ConsToPrim
    stateLW.flag = &local_flag;
    uint32_t local_err = ConsToPrim (stateLW, run_config);

    // 2e. Fallback approach: upon error, set 
    //     om = 0 --> use flux_lf[])

    if (local_err) {
      ATOMIC_ERR_UPDATE(err_kernel, local_err);
      KERNEL_ERR_UPDATE(data->err(ti,tj,tk), local_err);
      om = 0.0;
    } else {
      State state_lw;
      STATE_SET(state_lw, v_lw, u_lw, bgf, x1[ti], x2[tj], x3[tk]);
      Flux<dir> (state_lw, f_lw, p_lw, run_config);
    }

    #if HLL_SWITCH == YES
    if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
      om = 0.0;
    }
    #endif

    // 2f. Compute GFORCE as convex average  
    //      between LF and LW fluxes.

    NFLX_LOOP(nv) flux(i,nv) = om*f_lw[nv] + (1.0 - om)*flux_lf[nv];
    press(i) = om*p_lw + (1.0 - om)*0.5*(pl + pr);

    // 2g. Compute Powell's / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
//    Bn(i) = (flux(i,RHO) > 0.0 ? vl[BXn]:vr[BXn]);
    Bn(i) = 0.5*(vl[BXn]+vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    // 2h. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2i. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2j. Store emf-flux components for CT

    #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
    #if (CT_EMF_AVERAGE == ARITHMETIC) || (CT_EMF_AVERAGE == CT_CONTACT)
    CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == CT_FLUX
    CT_Store_Flux<dir>(data, fl, fr, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_HLL
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_GFORCE
    CT_Store_GForce<dir> (data, vl, vr, u_lw, v_lw, om, cmax(i), ti, tj, tk);
    #endif
    #endif

    // 2k. Conservative Energy correction step
    //      [Experimental]

    #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
    double EgL, EgR, Eg_lw;
    double fEgL   = HD_EnergyFlux(vl,   &EgL,   dir, run_config);
    double fEgR   = HD_EnergyFlux(vr,   &EgR,   dir, run_config);
    double fEg_lw = HD_EnergyFlux(v_lw, &Eg_lw, dir, run_config);
    double fEg_lf = 0.5*(fEgL + fEgR) - 0.5*cmax(i)*(EgR - EgL); 

    flux(i,ENG) = om*fEg_lw + (1.0 - om)*fEg_lf;
    #elif CT_ENERGY_CORRECTION == 4
    flux(i,ENG) = (bR*fEgL - bL*fEgR)*scrh; 
    #endif

  }}} // End for (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("GFORCE_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from GFORCE_Solver()",err_kernel);
  #endif

}
#undef CHECK_ROE_MATRIX
