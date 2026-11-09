////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief gPLUTO porting of the Implementation of the two-shock Riemann solver 
         for the relativistic hydro equations.
  
  Solve the Riemann problem for the relativistic hydrodynamics 
  equations using the two-shock approach described by Mignone, Plewa 
  and Bodo (2005).
  The formulation works with IDEAL or TAUB equation of state.

  On input, this function takes left and right primitive state vectors 
  \c stateL->v and \c stateR->v at zone edge \c i+1/2.
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "The Piecewise Parabolic Method for  Multidimensional Relativistic 
       Fluid Dynamics", Mignone, Plewa and Bodo, ApJS (2005) 160,199.
       
  \authors A. Mignone (andrea.mignone@unito.it), 
           S. Truzzi (stefano.truzzi@unito.it)
  \date    Jan 15, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define TWO_SHOCK_MAX_ITER   20
#define TWO_SHOCK_SMALL_PRS    1.e-12
#define TWO_SHOCK_SMALL_RHO    1.e-12
#define INTERPOLATE_RAREFACTION   YES

#define TWO_SHOCK_TOL   1.e-6

void   TwoShock_Shock (double, double, double, double, double, double,
                       double, double &, double &, double &,  
                       double, int);

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
template<int dir>
void TwoShock_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
/*!
  Solve the Riemann problem for the Euler equations of gasdynamics
  using the two-shock approximation.
  
  \param[in, out] d   Datastructure contains all the prim cons and fluxes
  \param[in,out] Dts  Datastructure contains the speeds
  \param[in] grid     Datastructure where the grid is defined
*/
////////////////////////////////////////////////////////////////////////
{ 
  int g_maxRiemannIter = 0;
  int      iter, i, j, k, nv;
  int      zone_fail;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  double gamma;

  RunConfig &run_config = *(data->run_config);

  #if EOS == IDEAL
  gamma = run_config.gamma;
  double gmmr  = gamma/(gamma - 1.0);
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  //    (set ni and nj (row and column indexes multiplier)
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  // Calculate the normal, tangent and bitanget direction

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);

  // Cycle on directional indexes to preserve sweep-> ordering

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
  
    // calculates offset based on cycle pos

    long int offset  = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D vC(&data->sweep->vC[offset1],ni,NVAR);

    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary1D press(&data->sweep->press[offset],ni);
    Ary1D cmax(&Dts->cmax[offset],ni);
    
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    // Private variable declaration

    int    ti, tj, tk;    
    int    iside;
    double ul[NVAR], ur[NVAR], uc[NVAR]; // conservative
    double ql[NVAR], qr[NVAR], qc[NVAR]; //
    double fl[NVAR], fr[NVAR], f_twoshock[NVAR]; // fluxes
    double Vstar[NFLX];

    double cminC, cmaxC;

    double a2L, duL, gL, gL1, hL, j2L, pL, scrhL, tauL, wL, wlL, uxL, VL, VL1;
    double a2R, duR, gR, gR1, hR, j2R, pR, scrhR, tauR, wR, wlR, uxR, VR, VR1;
    double am,ap;
    double press_twoshock;
    double p1, wl2;
    
    double vx, vt2, vel2, sroot, delta2;
    double a2local;
    double a0, beta_fix=0.9999;

    State stateL, stateR, state_star;

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    for (nv = 0; nv < NVAR; nv++) {
      ql[nv] = vL(i,nv);
      qr[nv] = vR(i,nv);
    }

    PrimToCons(ql, ul, run_config);
    PrimToCons(qr, ur, run_config);

    STATE_SET(stateL, ql, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, qr, ur, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(state_star, Vstar, NULL, NULL, x1[ti], x2[tj], x3[tk]);

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. Compute sound speed and Enthalpy at zone interf.

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);

    a2L = stateL.a2; hL = stateL.h;
    a2R = stateR.a2; hR = stateR.h;

    // 2b. Compute relevant quantities

    scrhL = ql[VX1]*ql[VX1] + ql[VX2]*ql[VX2] + ql[VX3]*ql[VX3];
    scrhR = qr[VX1]*qr[VX1] + qr[VX2]*qr[VX2] + qr[VX3]*qr[VX3];

    if (scrhL >= 1.0){
      //printf ("! TwoShock_Lorentz(): u2 > 1 (%12.6e) in TwoShock_Lorentz\n", scrhL);
      scrhL = beta_fix/sqrt(scrhL);
      ql[VX1] *= scrhL;
      ql[VX2] *= scrhL;
      ql[VX3] *= scrhL;
      scrhL = beta_fix*beta_fix;
    }

    if (scrhR >= 1.0){
      //printf ("! TwoShock_Lorentz(): u2 > 1 (%12.6e) in TwoShock_Lorentz\n", scrhR);
      scrhR = beta_fix/sqrt(scrhR);
      qr[VX1] *= scrhR;
      qr[VX2] *= scrhR;
      qr[VX3] *= scrhR;
      scrhR = beta_fix*beta_fix;
    }

    wlL = 1.0/(1.0 - scrhL);
    gL = sqrt(wlL);

    wlR = 1.0/(1.0 - scrhR);
    gR = sqrt(wlR);

    tauR = 1.0/qr[RHO];
    tauL = 1.0/ql[RHO];
      
    uxR = qr[VXn];
    uxL = ql[VXn];

    pR = qr[PRS];
    pL = ql[PRS];

    wR = pR*tauR;
    wL = pL*tauL;

    VR = tauR/gR;
    VL = tauL/gL;
    
    // 2c. First iteration outside of the loop

    #if EOS == IDEAL
    j2R = tauR*(hR*(gmmr - 2.0) + 1.0) / (gmmr*pR);
    j2L = tauL*(hL*(gmmr - 2.0) + 1.0) / (gmmr*pL);
    #elif EOS == TAUB
    a0  = (5.0*hL - 8.0*wL)/(2.0*hL - 5.0*wL);
    j2L = -tauL*(a0*wL + hL*(1.0 - a0))/(a0*pL);
      
    a0  = (5.0*hR - 8.0*wR)/(2.0*hR - 5.0*wR);
    j2R = -tauR*(a0*wR + hR*(1.0 - a0))/(a0*pR); 
    #endif

    gR1 = sqrt(VR*VR + (1.0 - uxR*uxR)*j2R);   //    RIGHT    
    wR  = VR*uxR + gR1;

    gL1 = -sqrt(VL*VL + (1.0 - uxL*uxL)*j2L);  //    LEFT    
    wL  = VL*uxL + gL1;
    
    duR = gR1/(hR*gR);
    duL = gL1/(hL*gL);
    
    p1  = ql[VXn] - qr[VXn] + duR*pR - duL*pL;
    p1 /= duR - duL;

    if (p1 < 0.0) {
      p1 = MIN(pR,pL);
    }
      
    double u1, dp;
    double uxL1, uxR1;

    // 2d. Begin Newton's iteration loop

    for (iter = 1; iter < TWO_SHOCK_MAX_ITER; iter++)  {

      TwoShock_Shock (tauL, uxL, pL, gL, VL, hL, p1, uxL1, duL, wL, gamma, -1);
      TwoShock_Shock (tauR, uxR, pR, gR, VR, hR, p1, uxR1, duR, wR, gamma,  1);

      dp  = (uxR1 - uxL1)/(duL - duR);
      if (-dp > p1) dp = -0.5*p1;
      p1 += dp;

      if (p1 < 0.0){
        p1 -= dp;
        p1 *= 0.5;
      }
      g_maxRiemannIter = MAX(g_maxRiemannIter,iter);
      if ( (fabs(dp) < TWO_SHOCK_TOL*p1) ) break;
    } 
    u1 = 0.5*(uxR1 + uxL1); 

    // 2e. Check possible failures and flag zones

    zone_fail = 0;
    if (isnan(u1) || iter == TWO_SHOCK_MAX_ITER) {
      u1 = 0.5*(uxL + uxR);
      p1 = 0.5*(pL  + pR);
      NFLX_LOOP(nv) qc[nv] = ql[nv];
      zone_fail = 1;
      ATOMIC_ERR_UPDATE(err_kernel,       ERR_RIEMANN_SOLVE_FAIL);
      KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_RIEMANN_SOLVE_FAIL);
    }

   // FLATTENING MULTID PART

    #if HLL_SWITCH == YES
    zone_fail += (flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL);
    #endif
  
    // 2f. No zone failure: proceed with two-shock Solver 
      
    if (zone_fail == 0){    

      Vstar[VXn] = u1;
      Vstar[PRS] = p1;

      // 2f1. Sample solution on  x/t = 0 axis

      if (u1 >= 0.0) {  

        // 2f1L. Sample solution LEFT of contact 

        gL1   = gL*hL/(gL*hL + (p1 - pL)*(VL + wL*uxL));
        Vstar[VXt] = ql[VXt]*gL1;
        Vstar[VXb] = ql[VXb]*gL1;

        VL1   = VL - (u1 - uxL)*wL;
        scrhL = Vstar[VX1]*Vstar[VX1] + Vstar[VX2]*Vstar[VX2] + Vstar[VX3]*Vstar[VX3];
        wl2   = 1.0/(1.0 - scrhL);
        gL1   = sqrt(wl2);

        Vstar[RHO] = MAX(TWO_SHOCK_SMALL_RHO,1.0/(VL1*gL1));

        // Shock or rarefaction ?
      
        #if INTERPOLATE_RAREFACTION == YES
        if (p1 < ql[PRS]){

          // Rarefaction tail 

          iside = -1;
          a2local = stateL.a2;

          vx   = ql[VXn];
          vt2  = ql[VXt]*ql[VXt] + ql[VXb]*ql[VXb];
          vel2 = vx*vx + vt2;

          sroot = a2local*(1.0 - vx*vx - vt2*a2local)*(1.0 - vel2);   // this is eta/gamma 
          
          if (sroot < 0.0){
          // printf ("! TwoShock_RarefactionSpeed(): sroot < 0 \n");
            ATOMIC_ERR_UPDATE(err_kernel,       ERR_RIEMANN_SOLVE_FAIL);
            KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_RIEMANN_SOLVE_FAIL);
            zone_fail = 1;
          }
          sroot = sqrt(sroot);   // this is eta/gamma 
          delta2   = 1.0 - vel2*a2local;
          am = ( (vx*(1.0 - a2local) + (double)iside*sroot)/delta2);

          // Rarefaction head

          iside = -1;

          SoundSpeed2 (state_star, run_config);
          a2local = state_star.a2;

          vx   = Vstar[VXn];
          vt2  = Vstar[VXt]*Vstar[VXt] + Vstar[VXb]*Vstar[VXb];
          vel2 = vx*vx + vt2;

          sroot = a2local*(1.0 - vx*vx - vt2*a2local)*(1.0 - vel2);   // this is eta/gamma 
          
          if (sroot < 0.0){
            //printf ("! TwoShock_RarefactionSpeed(): sroot < 0 \n");
            ATOMIC_ERR_UPDATE(err_kernel,       ERR_RIEMANN_SOLVE_FAIL);
            KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_RIEMANN_SOLVE_FAIL);
            zone_fail = 1;

          }
          sroot = sqrt(sroot);   //this is eta/gamma 
          delta2   = 1.0 - vel2*a2local;
          ap = ((vx*(1.0 - a2local) + (double)iside*sroot)/delta2);

          am  = MIN(am, ap);
        }else{        
          ap = am = VL/wL + uxL;                // -- Shock speed --  //
        }
        #else
        am = ap = VL/wL + uxL;                  //  -- Shock speed --  //
        #endif
        
        if (am >= 0.0) {                      // --  region L  -- 
          NFLX_LOOP(nv) qc[nv] = ql[nv];
        }else if (ap <= 0.0) {                // -- region L1 -- 
          NFLX_LOOP(nv) qc[nv] = Vstar[nv];
        }else{    //  Solution is inside rarefaction fan, --> interpolate  
          NFLX_LOOP(nv) {
            qc[nv] = (am*Vstar[nv] - ap*ql[nv])/(am - ap);
          }
        }
            
      } else {   

        // 2f1R. Sample solution RIGHT of contact 

        gR1 = gR*hR/(gR*hR + (p1 - pR)*(VR + wR*uxR));
        Vstar[VXt] = qr[VXt]*gR1;
        Vstar[VXb] = qr[VXb]*gR1;

        scrhR = Vstar[VX1]*Vstar[VX1] + Vstar[VX2]*Vstar[VX2] + Vstar[VX3]*Vstar[VX3];
        wl2 = 1.0/(1.0 - scrhR);
        gR1 = sqrt(wl2);

        VR1        = VR - (u1 - uxR)*wR;
        Vstar[RHO] = MAX(TWO_SHOCK_SMALL_RHO,1.0/(VR1*gR1));
        
        // Shock or rarefaction ?

        #if INTERPOLATE_RAREFACTION == YES
        if (p1 < qr[PRS]){

          // Rarefaction tail 
            
          iside = 1;
          SoundSpeed2 (state_star, run_config);
          a2local = state_star.a2;

          vx   = Vstar[VXn];
          vt2  = Vstar[VXt]*Vstar[VXt] + Vstar[VXb]*Vstar[VXb];
          vel2 = vx*vx + vt2;

          sroot = a2local*(1.0 - vx*vx - vt2*a2local)*(1.0 - vel2);   // this is eta/gamma 
          
          if (sroot < 0.0){
            //printf ("! TwoShock_RarefactionSpeed(): sroot < 0 \n");
            ATOMIC_ERR_UPDATE(err_kernel,       ERR_RIEMANN_SOLVE_FAIL);
            KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_RIEMANN_SOLVE_FAIL);
            zone_fail = 1;
          }
          sroot = sqrt(sroot);   //this is eta/gamma 
          delta2   = 1.0 - vel2*a2local;
          am = ( (vx*(1.0 - a2local) + (double)iside*sroot)/delta2);

          // Rarefaction head

          iside = 1;
          a2local = stateR.a2;

          vx   = qr[VXn];
          vt2  = qr[VXt]*qr[VXt] + qr[VXb]*qr[VXb];
          vel2 = vx*vx + vt2;

          sroot = a2local*(1.0 - vx*vx - vt2*a2local)*(1.0 - vel2);   // this is eta/gamma 
          
          if (sroot < 0.0){
            //printf ("! TwoShock_RarefactionSpeed(): sroot < 0 \n");
            ATOMIC_ERR_UPDATE(err_kernel,       ERR_RIEMANN_SOLVE_FAIL);
            KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_RIEMANN_SOLVE_FAIL);
            zone_fail = 1;

          }
          sroot = sqrt(sroot);   //this is eta/gamma 
          delta2   = 1.0 - vel2*a2local;
          ap = ( (vx*(1.0 - a2local) + (double)iside*sroot)/delta2);
          ap  = MAX(am, ap);
        }else{        
          am = ap = VR/wR + uxR;               // -- Shock speed -- //
        }
        #else
        am = ap = VR/wR + uxR;                 // -- Shock speed -- //
        #endif
            
        if (ap <= 0.0) {                        // -- Shock speed -- //
          NFLX_LOOP(nv) qc[nv] = qr[nv];
        }else if (am >= 0.0){              //   region R1   //
          NFLX_LOOP(nv) qc[nv] = Vstar[nv];
        }else{      //  Solution is inside rarefaction fan, --> interpolate  //
          NFLX_LOOP(nv) qc[nv] = (ap*Vstar[nv] - am*qr[nv])/(ap - am);
        }
      }

      // 2f2. Compute max wave speed

      PrimToCons (qc, uc, run_config);

      State stateC;
      STATE_SET(stateC, qc, uc, NULL, x1[ti], x2[tj], x3[tk]);
      SoundSpeed2(stateC, run_config);
      MaxSignalSpeed<dir> (stateC, cminC, cmaxC);
      cmax(i) = MAX(fabs(cminC), fabs(cmaxC));

      // 2f3. Compute flux from sampled solution
    
      Flux<dir>(stateC, f_twoshock, press_twoshock);

      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = f_twoshock[nv];
      press(i) = press_twoshock;

      #ifdef ARTIFICIAL_VISC
      a0 = ARTIFICIAL_VISC*(MAX(0.0, ql[VXn] - qr[VXn]));
      NFLX_LOOP(nv) {
        flux(i,nv) += a0*(ul[nv]- ur[nv]);
      }
      #endif      
    } // end if (zone_fail == 0)

    // 2g. Zone failure: replace flux with HLL flux.

    if (zone_fail > 0){       

      HLL_Speed<dir> (stateL, stateR, gL, gR);
      
      a0 = MAX(fabs(gR), fabs(gL));
      cmax(i) = a0;

      gL = MIN(0.0, gL);
      gR = MAX(0.0, gR);
      a0 = 1.0/(gR - gL);

      Flux<dir>(stateL, fl, pL);
      Flux<dir>(stateR, fr, pR);

      for (nv = 0; nv < NFLX; nv++){
        flux(i,nv) = (gL*gR*(ur[nv] - ul[nv]) + gR*fl[nv] - gL*fr[nv])*a0;
      }
      press(i) = (gR*pL - gL*pR)*a0;
    }

    // 2h. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*ql[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*qr[nv];
    }
    #endif

    // 2i. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel,       ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_NAN, ERR_QUIT); 
      }
    }

  }}} // End loop (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------
  
  #if WARNING_MESSAGES > 0 
  Where ("TwoShock_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  int quit_condition = (err_kernel & ERR_QUIT);

  #if FALLBACK_MODE != NO
  quit_condition *= (data->integration_mode != DEFAULT);
  Fallback_Trigger(data->err, data->flag, box);
  #endif

  if (quit_condition) QuitPLUTO("from TwoShock_Solver()",err_kernel);
  #endif

}

////////////////////////////////////////////////////////////////////////
void TwoShock_Shock (double tau0, double u0, double p0, double g0, 
                     double V0, double h0, double p1, double &u1,
                     double &dudp, double &zeta, double gamma, int isweep)
/*!
  Compute post shock quantities  u1, dudp, zeta, for a given
  value of the post-shock pressure p1
*/
////////////////////////////////////////////////////////////////////////
{
  double a,b,c,dp;
  double d_htau1,g1;
  double j2;

////////////////////////////////////////////////////////////////////////
/* Use Taub Adiabat to find post-shock Enthalpy
        and mass flux; here j2 --> 1/(j)^2                            */
////////////////////////////////////////////////////////////////////////

  dp = p1 - p0;

  #if EOS == IDEAL
  double gmmr;
  double tau1, h1;

  gmmr = gamma/(gamma - 1.0);
  a = 1.0 - dp/(gmmr*p1);
  b = 1.0 - a;
  c = -h0*(h0 + tau0*dp);
   
  h1 = 0.5/a*(-b + sqrt(b*b - 4.0*a*c));
  tau1 = (h1 - 1.0)/(gmmr*p1);
  g1   = 2.0*h1*gmmr/(2.0*h1 - 1.0);

  j2  = h0*gmmr*tau0 + (h1*tau1 + h0*tau0)*(1.0/(h0 + h1) - 1.0);
  j2 /= gmmr*p1;

  d_htau1  = (h1*tau1 + h0*tau0 - g1*h1*tau1);
  d_htau1 /= (g1*p1 - dp);

  #elif EOS == TAUB
  double dx,da,db,dc;

  a = p0*(3.0*p1 + p0);
  b = -h0*h0*(3.0*p1 + 2.0*p0)
      -h0*tau0*(3.0*p1*p1 - 7.0*p1*p0 - 4.0*p0*p0) - dp;
  c = -h0*tau0*(h0*h0 + 2.0*h0*tau0*p1 + 2.0);
   
  dx = 2.0*c*dp/(-b + sqrt(b*b - 4.0*a*c*dp));// = [h\tau] 

  g1 = 2.0*c/(-b + sqrt(b*b - 4.0*a*c*dp));// = [h\tau]/[dp] 
  // h1 = sqrt(h0*h0 + (dx + 2.0*h0*tau0)*dp);

  j2 = -g1;

  da = 3.0*p0;
  db = -h0*h0*3.0 - h0*tau0*(6.0*p1 - 7.0*p0) - 1.0;
  dc = c - h0*tau0*dp*2.0*h0*tau0;

  d_htau1 = -(da*dx*dx + db*dx + dc)/(2.0*a*dx + b);
  #endif

  g1    = ((double)isweep)*sqrt(V0*V0 + (1.0 - u0*u0)*j2);
  zeta = (V0*u0 + g1) / (1.0 - u0*u0);

////////////////////////////////////////////////////////////////////////
/* Get post-shock velocity                                            */
////////////////////////////////////////////////////////////////////////

  b   = 1.0/(h0*g0 + (p1 - p0)*((zeta)*u0 + V0));
  u1 = (h0*g0*u0 + (zeta)*(p1 - p0)) * b;

////////////////////////////////////////////////////////////////////////
/* Get du/dp for next iteration                                       */
////////////////////////////////////////////////////////////////////////

  a     = -0.5*(d_htau1 + j2) / g1;
  dudp = ((zeta) + a - (u1)*((zeta)*u0 + V0 + a*u0)) * b;
}

#undef  TWO_SHOCK_MAX_ITER
#undef  TWO_SHOCK_SMALL_PRS
#undef  TWO_SHOCK_SMALL_RHO
#undef  TWO_SHOCK_TOL 
#undef  INTERPOLATE_RAREFACTION 
