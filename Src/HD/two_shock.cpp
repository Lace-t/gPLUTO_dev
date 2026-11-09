////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Two-shock Riemann solver for the HD equations.

  Solve the Riemann problem for the Euler equations of gasdynamics
  using the two-shock approximation.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also, compute maximum wave propagation speed (cmax)
  for  explicit time step computation

  \b Reference:
   -  "FLASH: an Adaptive Mesh Hydrodynamics Code for Modeling
      Astrophysical Thermonuclear Flashes"
      Fryxell et al, ApJS(2000) 131:273 (pages 292-294)

  \authors A. Mignone (andrea.mignone@unito.it)\n
           S. Truzzi (stefano.truzzi@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define MAX_ITER      5
#define small_p       1.e-9
#define small_rho     1.e-9
//#define DEBUG_TRIAL   1

template<int dir>
void TwoShock_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
void TwoShock_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  if      (dir == IDIR) TwoShock_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) TwoShock_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) TwoShock_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
  Solve the Riemann problem for the Euler equations of gasdynamics
  using the two-shock approximation.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void TwoShock_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  // ------------------------------------------------------
  // 0a. Compatibility checks
  // ------------------------------------------------------

#if EOS != IDEAL
  QuitPLUTO("! Two_Shock Solver: available only for IDEAL EoS", 1);
#else

  // ------------------------------------------------------
  // 0b. Allocate local variables and pointers
  // ------------------------------------------------------

  int i, j, k, nv;
  int iter, g_maxRiemannIter = MAX_ITER;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int      ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 
  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 
  RunConfig &run_config = *(data->run_config);

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  // Calculate the normal, tangent and bitanget direction

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  // Cycle on directional indexes to preserve sweep-> ordering

  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  double  gamma  = run_config.gamma;
  double  g1_g;

  g1_g = 0.5*(gamma + 1.0)/gamma;
  
  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config, Dts, grid) \
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  
    // Compute offsets and assign local pointers

    long int offset  = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary1D cs2(&data->sweep->cs2[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D vC(&data->sweep->vC[offset1],ni,NVAR);

    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary1D press(&data->sweep->press[offset],ni);

    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    // Private variable declaration

    int ti, tj, tk;

    double scrh1, scrh2, scrh3, scrh4;
    double rho_star, cstar;
    double lambda_s, lambda_star;
    double dp, zeta;
    double pstar;
    double taul, taur;

    double vl[NVAR], vr[NVAR], vc[NVAR]; // primitive
    double f_twoshock[NVAR]; // fluxes
    double ul[NVAR], ur[NVAR], uc[NVAR]; // conservative
    double ql[NVAR], qr[NVAR], qs[NVAR]; //

    double sigma, taus, cs, zs, ustar;
    double cl, cr, cmax_twoshock, a2;
    double press_twoshock;
    double vxl, vxr;

    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR)  {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR)  {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR)  {ti = j; tj = k; tk = i;}

    NFLUID_LOOP(nv) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. Compute conservative variables

    PrimToCons(vl, ul, run_config);
    PrimToCons(vr, ur, run_config);

    // 2b. Compute HLL Flux in case of strong shock

    #if HLL_SWITCH == YES
    double flux_hll[NVAR];
    double fl[NVAR], fr[NVAR];
    double pl, pr;
    double press_hll;
    double cmax_hll;
    if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){

      stateL.a2 = gamma*vl[PRS]/vl[RHO];
      stateR.a2 = gamma*vr[PRS]/vr[RHO];
      HLL_Speed<dir> (stateL, stateR, cl, cr);

      scrh1 = MAX(fabs(cl), fabs(cr));
      cmax_hll = scrh1;

      cl = MIN(0.0, cl);
      cr = MAX(0.0, cr);

      scrh1 = 1.0 / (cr - cl);

      Flux<dir> (stateL, fl, pl);
      Flux<dir> (stateR, fr, pr);

      for (nv = 0; nv < NFLX; nv++){
        flux_hll[nv]  = cl*cr*(ur[nv] - ul[nv]) + cr*fl[nv] - cl*fr[nv];
        flux_hll[nv] *= scrh1;
      }

      press_hll = (cr*pl - cl*pr)*scrh1;

    }
    #endif // HLL_SWITCH

    // ------------------------------------------
    // 2c. Two-shock solver
    // ------------------------------------------

    NFLUID_LOOP(nv) {
      qr[nv]=vr[nv];
      ql[nv]=vl[nv];
    }

    cl = sqrt(gamma*ql[PRS]*ql[RHO]);
    cr = sqrt(gamma*qr[PRS]*qr[RHO]);

    taul = 1.0/ql[RHO];
    taur = 1.0/qr[RHO];

    // 2c1. Provide 1st guess

    pstar = qr[PRS] - ql[PRS] - cr*(qr[VXn] - ql[VXn]);
    pstar = ql[PRS] + pstar*cl/(cl + cr);
    pstar = MAX(small_p , pstar);

    // 2c2. Newton iterative cycle

    for(iter = 1; iter <= MAX_ITER; iter++){
      vxl = cl*sqrt(1.0 + g1_g*(pstar - ql[PRS])/ql[PRS]);
      vxr = cr*sqrt(1.0 + g1_g*(pstar - qr[PRS])/qr[PRS]);

      scrh1 = vxl*vxl;
      scrh1 = 2.0*scrh1*vxl/(scrh1 + cl*cl);

      scrh2 = vxr*vxr;
      scrh2 = 2.0*scrh2*vxr/(scrh2 + cr*cr);

      scrh3 = ql[VXn] - (pstar - ql[PRS])/vxl;
      scrh4 = qr[VXn] + (pstar - qr[PRS])/vxr;

      dp = scrh1*scrh2/(scrh1 + scrh2)*(scrh4 - scrh3);

      pstar -= dp;

      pstar = MAX (small_p, pstar);
      g_maxRiemannIter = MAX(g_maxRiemannIter,iter);
      if (fabs(dp/pstar) < 1.e-6) break;
    }

    scrh3 = ql[VXn] - (pstar - ql[PRS])/vxl;
    scrh4 = qr[VXn] + (pstar - qr[PRS])/vxr;
    ustar = 0.5*(scrh3 + scrh4);

    // 2c3. Select correct * region

    if (ustar > 0.0) {
      sigma = 1.0;
      taus  = taul;
      cs    = cl*taul;
      zs    = vxl;
      NFLUID_LOOP(nv) qs[nv] = ql[nv];
    }else{
      sigma = -1.0;
      taus  = taur;
      cs    = cr*taur;
      zs    = vxr;
      NFLUID_LOOP(nv) qs[nv] = qr[nv];
    }

    rho_star = taus - (pstar - qs[PRS])/(zs*zs);
    rho_star = MAX(small_rho,1.0/rho_star);
    cstar    = sqrt(gamma*pstar/rho_star);
    if (pstar < qs[PRS]){
      lambda_s    = cs    - sigma*qs[VXn];
      lambda_star = cstar - sigma*ustar;
    }else{
      lambda_s    = lambda_star = zs*taus  - sigma*qs[VXn];
    }

    if(lambda_star > 0.0){
      vc[RHO] = rho_star;
      vc[VXn] = ustar;
      vc[PRS] = pstar;
    } else if (lambda_s < 0.0){
      vc[RHO] = qs[RHO];
      vc[VXn] = qs[VXn];
      vc[PRS] = qs[PRS];

    } else{  //  Linearly interpolate rarefaction fan 

      scrh1 = MAX(lambda_s - lambda_star, lambda_s + lambda_star);
      scrh1 = MAX(1.e-12,scrh1); // AM (26/01/2017): avoid division by 0.0 at
                                    // sonic rarefactions when
                                    // lambda_s = lambda_star = 0.0 
      zeta  = 0.5*(1.0 + (lambda_s + lambda_star)/scrh1);

      vc[RHO] = zeta*rho_star + (1.0 - zeta)*qs[RHO];
      vc[VXn] = zeta*ustar    + (1.0 - zeta)*qs[VXn];
      vc[PRS] = zeta*pstar    + (1.0 - zeta)*qs[PRS];
    }

    vc[VXt] = qs[VXt];
    vc[VXb] = qs[VXb];

    PrimToCons (vc, uc, run_config);

    a2 = gamma*vc[PRS]/vc[RHO];
    cstar = sqrt(a2);
    scrh1 = fabs(vc[VXn])/cstar;
    scrh1 = fabs(vc[VXn]) + cstar;
    cmax_twoshock = scrh1;

    {
      State state2s;
      STATE_SET(state2s, vc, uc, NULL, x1[ti], x2[tj], x3[tk]);
      Flux<dir>(state2s, f_twoshock, press_twoshock);
    }

    #if HLL_SWITCH == YES
    if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = flux_hll[nv];
      press(i) = press_hll;
      cmax(i) = cmax_hll;
    } else {
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = f_twoshock[nv];
      press(i) = press_twoshock;
      cmax(i) = cmax_twoshock;
    }

    #else
    for (nv = 0; nv < NFLX; nv++) flux(i,nv) = f_twoshock[nv];
    press(i) = press_twoshock;
    cmax(i) = cmax_twoshock;
    #endif

    cs2(i)  = a2;

    // 2d. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

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

  }}} // End Domain loop (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("TwoShock_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from TwoShock_Solver()",err_kernel);
  #endif

#endif // EOS check
}
