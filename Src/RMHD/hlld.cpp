////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLLD Riemann solver for the relativistic MHD equations.

  Solve the Riemann problem for the relativistic MHD equations
  using the HLLD solver of Mignone, Ugliano & Bodo (2009).
  The solver requires the solution of a nonlinear equation (Eq. [48])
  and the maximum number of iterations before convergence can be
  set using the macro HLLD_MAX_ITER (default 20).
  A physical relevant solution is accepted if it satisfies
  the constraints by Eq. [54], implemented in the function HLLD_Fstar().
  If this is not the case or if other failure conditions occur during
  the iteration loop, we switch to the HLL Riemann solver.

  The macro ::COUNT_FAILURES can be set to YES if one wishes to count
  how many times the solver fails.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2.
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
    - "A five wave Harte-Lax-van Leer Riemann solver for relativistic
       magnetohydrodynamics" Mignone et al, MNRAS (2009) 393,1141.

  \authors A. Mignone  (andrea.mignone@unito.it)  \n 
           M. Rossazza (marco.rossazza@unito.it)  \n
           A. Suriano  (alessio.suriano@unito.it)
           M. Bugli    (matteo.bugli@cea.fr)
  \date    Jun 09, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define HLLD_ITERATIONS   -1

#define MAX_ITER  20
#define COUNT_FAILURES  NO /**< When set to YES, count number of failures and
        write the count to "hlld_fails.dat" (works in parallel as well).
        The number of failures is normalized to the total number of zones:
        \c nf/nz where \c nz is incremented during the main loop and \c nf
        is incremented when a failure occur. */

typedef struct RIEMANN_STATE{
  int    fail;
  double vx, vy, vz;
  double Bx, By, Bz;
  double Kx, Ky, Kz, K2;
  double w, sw, p, rho;
  double u[NVAR], R[NVAR];
  double S, Sa;
} Riemann_State;

template<int dir>
double HLLD_Fstar (Riemann_State *, Riemann_State *, double, double, double *, Grid *);

template<int dir>
void HLLD_GetCState (Riemann_State *, Riemann_State *, double, double, double *, Grid *);

template<int dir>
int HLLD_GetRiemannState (Riemann_State *, double, double, int, Grid *);

template<int dir>
void HLLD_GetAState (Riemann_State *, double, double, Grid *);

template<int dir>
void HLLD_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box);

#pragma acc routine seq
double HLLD_TotalPressure (double *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void HLLD_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) HLLD_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) HLLD_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) HLLD_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLLD_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
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

  #if (defined STAGGERED_MHD) && (CT_EMF_AVERAGE == UCT_GFORCE)
  QuitPLUTO("! HLLD Solver : not compatible with UCT_GFORCE", 1);
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  constexpr int MXn = MX1+dir;
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
  #pragma acc parallel loop collapse(2) present(data, run_config, grid, Dts)\
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  #pragma acc loop
  for (i = ibeg; i <= iend; i++){

    // 1a. Compute offsets and assign local pointers

    long int offset = ni*(j + nj*k);
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

    double Uhll[NVAR];
    double Fhll[NVAR];
    double Vhll[NVAR];
    double Uc[NVAR];

    // 2. Loop over normal direction

    // #pragma acc loop private(Uhll[:NVAR], Fhll[:NVAR], Vhll[:NVAR], Uc[:NVAR])   
    //  for (i = ibeg; i <= iend; i++){

    int ti, tj, tk;
    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl, a2L;
    double vr[NVAR], ur[NVAR], fr[NVAR], Sr, pr, a2R;
    double scrh;
    State stateL, stateR, stateHLL;

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    #pragma acc loop seq
    for (nv = 0; nv < NVAR; nv++) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    data->err(ti,tj,tk) = 0;   // reset error 

    #if PHYSICS == GRMHD
    Metric G = DEFAULT_METRIC;
    if (!SingularMetric(x1[ti],x2[tj],x3[tk])){
      ComputeLocalMetric(x1[ti], x2[tj], x3[tk], &G, run_config);
      ComputeSqrtDet(&G);
      ComputeGcon(&G);
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

    // 2b. Compute sound speed, cons variables, fluxes
    //     and HLL signal speeds.

    STATE_SET(stateL,    vl,    ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR,    vr,    ur, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateHLL, Vhll, Uhll, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    a2L = stateL.a2;
    a2R = stateR.a2;
    #if PHYSICS != GRMHD
    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);
    #else
    Flux<dir> (stateL, fl, pl, run_config, &G);
    Flux<dir> (stateR, fr, pr, run_config, &G);

    HLL_Speed<dir> (stateL, stateR, Sl, Sr, &G);
    #endif

    // 2c. Get max propagation speed & sound speed squared

    cmax(i) = MAX(fabs(Sl), fabs(Sr));

    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    // 2d. Solve equation iteratively
    //     Compute flux directly if SL > 0 or SR < 0

    if (Sl >= 0.0){

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fl[nv];
      press(i) = pl;

    }else if (Sr <= 0.0){

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fr[nv];
      press(i) = pr;

    }else{

      // 2e. Compute HLLD Flux when  SL < 0 < SR

      int switch_to_hll = 0;
      double p0, f0;
      double dS_1 = 1.0/(Sr - Sl);
      double p, f, dp, Sc;

      // 2e1. Build the HLL average state & flux
      //      Entropy is necessary for next guess 

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        Uhll[nv]  = Sr*ur[nv] - Sl*ul[nv] + fl[nv] - fr[nv];
        Uhll[nv] *= dS_1;

        Fhll[nv]  = Sr*fl[nv] - Sl*fr[nv] + Sl*Sr*(ur[nv] - ul[nv]);
        Fhll[nv] *= dS_1;
      }
      #if ENTROPY_SWITCH
      {
        double vxR = vr[VXn];
        double vxL = vl[VXn];
        nv = ENTR;
        Uhll[nv]  = (Sr - vxR)*ur[nv] - (Sl - vxL)*ul[nv];
        Uhll[nv] *= dS_1;

        Fhll[nv]  =   Sr*ul[nv]*vl[VXn] - Sl*ur[nv]*vr[VXn] 
                    + Sl*Sr*(ur[nv] - ul[nv]);
        Fhll[nv] *= dS_1;
      }
      #endif

      // 2e2. Switch to HLL at strong shocks 

      #if HLL_SWITCH == YES
      if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
        switch_to_hll = 1;
      }
      #endif

      // 2e3. Build the right hand sides
      //      (lambda*U - F) for L/R states

      Riemann_State PaL, PaR;
      PaL.S = Sl;
      PaR.S = Sr;
      double Bx    = Uhll[BXn];
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        PaL.R[nv] = Sl*ul[nv] - fl[nv];
        PaR.R[nv] = Sr*ur[nv] - fr[nv];
      }
      #if RMHD_REDUCED_ENERGY == YES
      PaL.R[ENG] += PaL.R[RHO];
      PaR.R[ENG] += PaR.R[RHO];
      #endif

      // 2e4. Provide initial guess

      scrh = MAX(pl, pr);
      if (Bx*Bx/scrh < 0.01) { // Try the B->0 limit 
        double a,b,c;

        a = Sr - Sl;
        b = PaR.R[ENG] - PaL.R[ENG] + Sr*PaL.R[MXn] - Sl*PaR.R[MXn];
        c = PaL.R[MXn]*PaR.R[ENG] - PaR.R[MXn]*PaL.R[ENG];
        scrh = b*b - 4.0*a*c;
        scrh = MAX(scrh,0.0);
        p0 = 0.5*(- b + sqrt(scrh))*dS_1;

      }else{  // Use HLL average 

        // double Uav[NVAR];
        // #pragma acc loop
        // NFLX_LOOP(nv) Uav[nv]  = 0.5*(ur[nv] + ul[nv]  + fl[nv] - fr[nv]);
        // uint32_t local_err = ConsToPrim(Uav, Vhll, &data->flag(i,j,k));
        // uint32_t local_err = ConsToPrim(Uhll, Vhll, run_config, &data->flag(ti,tj,tk));

        stateHLL.flag = &data->flag(ti,tj,tk);
        #if PHYSICS != GRMHD
        uint32_t local_err = ConsToPrim(stateHLL, run_config);
        #else
        uint32_t local_err = ConsToPrim(stateHLL, run_config, &G);
        #endif
        if (local_err) {
          ATOMIC_ERR_UPDATE(err_kernel, local_err);
          KERNEL_ERR_UPDATE(data->err(ti,tj,tk), local_err);
        }

        if (local_err) switch_to_hll = 1;
        p0    = HLLD_TotalPressure(Vhll);
      }

      // 2e5. Check if guess is physically ok.
      //      If not, switch to HLL.

      f0 = HLLD_Fstar<dir> (&PaL, &PaR, p0, Bx, &Sc, grid);
      if (isnan(f0)|| PaL.fail) switch_to_hll = 1;

      // 2e6. Start root solver

      if (fabs(f0) > 1.e-12 && !switch_to_hll){

        p  = 1.025*p0; f  = f0;
        for (int k1 = 1; k1 < MAX_ITER; k1++){
          f  = HLLD_Fstar<dir> (&PaL, &PaR, p, Bx, &Sc, grid);
          if (f != f  || PaL.fail || (k1 > 7) || (fabs(f) > fabs(f0) && k1 > 4)) {
            switch_to_hll = 1;
            break;
          }

          dp = (p - p0)/(f - f0)*f;
          p0 = p; f0 = f;
          p -= dp;
          if (p < 0.0) p = 1.e-6;
          if (fabs(dp) < 1.e-5*p || fabs(f) < 1.e-6) break;
        }
      } else p = p0;

      // #if HLLD_DEBUG == YES
      // HLLD_PrintWhatsWrong(&PaL, &PaR, phll, phll, p0Bx, p, vL, vR);
      // #endif

      // Too many iterations ? --> switch to HLL

      if (PaL.fail) switch_to_hll = 1;

      // 2e7. Comput final flux in the * region

      if (!switch_to_hll){   // -- Solution should be reliable -- 

        // g_maxRiemannIter = MAX(g_maxRiemannIter, k);

        #ifdef GLM_MHD
        PaL.u[PSI_GLM] = PaR.u[PSI_GLM] = Uc[PSI_GLM] = ul[PSI_GLM];
        #endif

        if (PaL.Sa >= -1.e-6){

          HLLD_GetAState<dir> (&PaL, p, Bx, grid);
          #if RMHD_REDUCED_ENERGY == YES
          PaL.u[ENG] -= PaL.u[RHO];
          #endif

          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++){
            flux(i,nv) = fl[nv] + Sl*(PaL.u[nv] - ul[nv]);
          }
          press(i) = pl;

        }else if (PaR.Sa <= 1.e-6){

          HLLD_GetAState<dir> (&PaR, p, Bx, grid);
          #if RMHD_REDUCED_ENERGY == YES
          PaR.u[ENG] -= PaR.u[RHO];
          #endif

          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++){
            flux(i,nv) = fr[nv] + Sr*(PaR.u[nv] - ur[nv]);
          }
          press(i) = pr;

        }else{

          HLLD_GetCState<dir> (&PaL, &PaR, p, Bx, Uc, grid);

          if (Sc > 0.0){

            #if RMHD_REDUCED_ENERGY == YES
            PaL.u[ENG] -= PaL.u[RHO];
            Uc[ENG]    -= Uc[RHO];
            #endif
            #pragma acc loop seq
            for (nv = 0; nv < NFLX; nv++){
              flux(i,nv) = fl[nv] + Sl*(PaL.u[nv] - ul[nv])
                                      + PaL.Sa*(Uc[nv] - PaL.u[nv]);
            }
            press(i) = pl;

          }else{

            #if RMHD_REDUCED_ENERGY == YES
            PaR.u[ENG] -= PaR.u[RHO];
            Uc[ENG]    -= Uc[RHO];
            #endif
            #pragma acc loop seq
            for (nv = 0; nv < NFLX; nv++){
              flux(i,nv) = fr[nv] + Sr*(PaR.u[nv] - ur[nv])
                                      + PaR.Sa*(Uc[nv] - PaR.u[nv]);
            }
            press(i) = pr;
          }
        }  // PaL.Sa, Par.Sa 

      } else {   // Switch to HLL flux

        #if COUNT_FAILURES == YES
        totfail += 1.0;
        #endif

        for (nv = 0; nv < NFLX; nv++) flux(i,nv) = Fhll[nv];
        press(i)  = (Sr*pl - Sl*pr)*dS_1;

      } // End if !switch_to_hll 
    } // End if on SL, SR -- */

    // 2f. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2g. Compute Powell / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = (flux(i,RHO) > 0.0 ? vl[BXn]:vr[BXn]);
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

    // 2i. Store emf-flux components for CT 

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
  Where ("HLLD_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from HLLD_Solver()",err_kernel);
  #endif

}
#undef MAX_ITER
#undef COUNT_FAILURES

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
double HLLD_Fstar (Riemann_State *PaL, Riemann_State *PaR, double p, 
                   double Bx, double *Sc, Grid *grid)
/*!
  Return the velocity difference across the contact mode as a
  function of the total pressure p.
*/
////////////////////////////////////////////////////////////////////////
{
  int    success = 1;
  double dK, Bxc, Byc, Bzc;
  double fun;
  double vxcL, KLBc;
  double vxcR, KRBc;

  // sBx = Bx > 0.0 ? 1.0:-1.0;

  success *= HLLD_GetRiemannState<dir> (PaL, p, Bx, -1, grid);
  success *= HLLD_GetRiemannState<dir> (PaR, p, Bx,  1, grid);

  // Compute B from average state

  dK  = PaR->Kx - PaL->Kx + 1.e-12;

  Bxc = Bx*dK;

  Byc =   PaR->By*(PaR->Kx - PaR->vx)
        - PaL->By*(PaL->Kx - PaL->vx)
        + Bx*(PaR->vy - PaL->vy);

  Bzc =   PaR->Bz*(PaR->Kx - PaR->vx)
        - PaL->Bz*(PaL->Kx - PaL->vx)
        + Bx*(PaR->vz - PaL->vz);

  KLBc = PaL->Kx*Bxc + PaL->Ky*Byc + PaL->Kz*Bzc;
  KRBc = PaR->Kx*Bxc + PaR->Ky*Byc + PaR->Kz*Bzc;

  vxcL = PaL->Kx - dK*Bx*(1.0 - PaL->K2)/(PaL->sw*dK - KLBc);
  vxcR = PaR->Kx - dK*Bx*(1.0 - PaR->K2)/(PaR->sw*dK - KRBc);

  PaL->Sa = PaL->Kx;
  PaR->Sa = PaR->Kx;
  *Sc      = 0.5*(vxcL + vxcR);
  fun     = vxcL - vxcR;   // Essentially, Eq. [48]

/*
fun = dK*(1.0 - Bx*(  (1.0 - PaR->K2)/(PaR->sw*dK - KRBc)
                    - (1.0 - PaL->K2)/(PaL->sw*dK - KLBc)) );
*/

  // Check if state makes physically sense 

  success  = (vxcL - PaL->Kx)  > -1.e-6;
  success *= (PaR->Kx  - vxcR) > -1.e-6;

  success *= (PaL->S - PaL->vx) < 0.0;
  success *= (PaR->S - PaR->vx) > 0.0;

  success *= (PaR->w - p) > 0.0;
  success *= (PaL->w - p) > 0.0;
  success *= (PaL->Sa - PaL->S)  > -1.e-6;
  success *= (PaR->S  - PaR->Sa) > -1.e-6;

  PaL->fail = !success;

/*
scrh  = (1.0 - PaR->K2)*(PaL->sw*dK - KLBc);
scrh -= (1.0 - PaL->K2)*(PaR->sw*dK - KRBc);

PaL->fun1 = (PaR->sw*dK - KRBc)*(PaL->sw*dK - KLBc) - Bx*scrh;
PaL->fun2 = scrh;
PaL->denL = (PaL->sw*dK - KLBc);
PaL->denR = (PaR->sw*dK - KRBc);
*/
  return (fun);
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
int HLLD_GetRiemannState (Riemann_State *Pv, double p, double Bx, int side, 
                          Grid *grid)
/*!
  Express the sweep behind a wave as function of the total
  pressure p and the right hand side on the other side of the wave.
  
  On output, return 1 if succesful, 0 if w < 0 is encountered.
  side = -1 : means left
  side =  1 : means right
*/
////////////////////////////////////////////////////////////////////////
{
  double A, C, G, X, s;
  double vx, vy, vz, scrh, S, *R;
  
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);

  S = Pv->S;
  R = Pv->R;

  A = R[MXn] + p*(1.0 - S*S) - S*R[ENG];         // Eq. [26] 
  G = R[BXt]*R[BXt] + R[BXb]*R[BXb];             // Eq. [27] 
  C = R[BXt]*R[MXt] + R[BXb]*R[MXb];             // Eq. [28] 
  X = Bx*(A*S*Bx + C) - (A + G)*(S*p + R[ENG]);  // Eq. [30] 

  // Compute the numerators of Eqs. [23,23,45] 

  vx = ( Bx*(A*Bx + C*S) - (R[MXn] + p)*(G + A) );

  vy = ( - (A + G - Bx*Bx*(1.0 - S*S))*R[MXt]
         + R[BXt]*(C + Bx*(S*R[MXn] - R[ENG])) );

  vz = ( - (A + G - Bx*Bx*(1.0 - S*S))*R[MXb]
         + R[BXb]*(C + Bx*(S*R[MXn] - R[ENG])) );

  scrh = vx*R[MXn] + vy*R[MXt] + vz*R[MXb];
  scrh = X*R[ENG] - scrh;

  Pv->w = p + scrh/(X*S - vx);  // Eq [31] 

  if (Pv->w < 0.0) return(0);  // Failure

  Pv->vx = vx/X;
  Pv->vy = vy/X;
  Pv->vz = vz/X;

  // -- original Eq. [21] --
/*
  Pv->Bx = Bx;
  Pv->By = (R[BXt]*X - Bx*vy)/(X*S - vx);
  Pv->Bz = (R[BXb]*X - Bx*vz)/(X*S - vx);
*/
/* --------------------------------------------------------------------- */
/*! When computing Bx, By and Bz, we use Eq. [21] with
    vx, vy and vz replaced by by Eq. [23,24,25].
    The following short MAPLE script is used to verify the
    correctness.
 \code
   restart;
   A    := R[mx] - lambda*R[E] + p*(1 - lambda^2):
   G    := R[By]*R[By] + R[Bz]*R[Bz]:
   C    := R[my]*R[By] + R[mz]*R[Bz]:
   Q    := - A - G + B[x]*B[x]*(1-lambda^2):
   X    := B[x]*(A*lambda*B[x]+ C) - (A + G)*(lambda*p + R[E]):
   v[x] := (B[x]*(A*B[x] + lambda*C) - (A + G)*(p + R[mx]))/X:
   v[y] := (Q*R[my] + R[By]*(C + B[x]*(lambda*R[mx] - R[E])))/X:
   B[y] := (R[By] - B[x]*v[y])/(lambda - v[x]):
   simplify(B[y]);
 \endcode
*/
// --------------------------------------------------------------------- 

  Pv->Bx = Bx;
  Pv->By = -(R[BXt]*(S*p + R[ENG]) - Bx*R[MXt])/A;
  Pv->Bz = -(R[BXb]*(S*p + R[ENG]) - Bx*R[MXb])/A;

  s  = Bx > 0.0 ? 1.0:-1.0;
  if (side < 0) s *= -1.0;

  Pv->sw = s*sqrt(Pv->w);

  scrh = 1.0/(S*p +  R[ENG] + Bx*Pv->sw);
  Pv->Kx = scrh*(R[MXn] + p + R[BXn]*Pv->sw);
  Pv->Ky = scrh*(R[MXt]     + R[BXt]*Pv->sw);
  Pv->Kz = scrh*(R[MXb]     + R[BXb]*Pv->sw);

  Pv->K2 = Pv->Kx*Pv->Kx + Pv->Ky*Pv->Ky + Pv->Kz*Pv->Kz;

  return(1); // -- success -- 
}
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
template<int dir>
void HLLD_GetAState (Riemann_State *Pa,  double p, double Bx, Grid *grid)
/*!
  Compute states aL and aR behind fast waves.
*/
////////////////////////////////////////////////////////////////////////
{
  double vB, *ua, *R, S;
  double scrh;

  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);

  ua = Pa->u;
  S  = Pa->S;
  R  = Pa->R;

  scrh = 1.0/(S - Pa->vx);

  ua[RHO] =  R[RHO]*scrh;
  ua[BXn] =  Bx;
  ua[BXt] = (R[BXt] - Bx*Pa->vy)*scrh;
  ua[BXb] = (R[BXb] - Bx*Pa->vz)*scrh;

  vB     = Pa->vx*ua[BXn] + Pa->vy*ua[BXt] + Pa->vz*ua[BXb];
  ua[ENG] = (R[ENG] + p*Pa->vx - vB*Bx)*scrh;

  ua[MXn] = (ua[ENG] + p)*Pa->vx - vB*ua[BXn];
  ua[MXt] = (ua[ENG] + p)*Pa->vy - vB*ua[BXt];
  ua[MXb] = (ua[ENG] + p)*Pa->vz - vB*ua[BXb];
  
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
template<int dir>
void HLLD_GetCState (Riemann_State *PaL, Riemann_State *PaR, double p,
                     double Bx, double *Uc, Grid *grid)
/*!
  Compute states cL and cR across contact mode.
*/
////////////////////////////////////////////////////////////////////////
{
  double dK, *ua;
  double vxcL, vycL, vzcL, KLBc;
  double vxcR, vycR, vzcR, KRBc;
  double vxc, vyc, vzc, vBc;
  double Bxc, Byc, Bzc, Sa, vxa;
  double scrhL, scrhR;

  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);

  HLLD_GetAState<dir> (PaL, p, Bx, grid);
  HLLD_GetAState<dir> (PaR, p, Bx, grid);
  dK = (PaR->Kx - PaL->Kx) + 1.e-12;

  Bxc = Bx*dK;

  Byc =   PaR->By*(PaR->Kx - PaR->vx)    // Eq. [45] 
        - PaL->By*(PaL->Kx - PaL->vx)
        + Bx*(PaR->vy - PaL->vy);

  Bzc =   PaR->Bz*(PaR->Kx - PaR->vx)    // Eq.[45] 
        - PaL->Bz*(PaL->Kx - PaL->vx)
        + Bx*(PaR->vz - PaL->vz);

  Bxc  = Bx;
  Byc /= dK;
  Bzc /= dK;

  Uc[BXn] = Bxc;
  Uc[BXt] = Byc;
  Uc[BXb] = Bzc;

  KLBc = PaL->Kx*Bxc + PaL->Ky*Byc + PaL->Kz*Bzc;
  KRBc = PaR->Kx*Bxc + PaR->Ky*Byc + PaR->Kz*Bzc;

  scrhL = (1.0 - PaL->K2)/(PaL->sw - KLBc);
  scrhR = (1.0 - PaR->K2)/(PaR->sw - KRBc);

  vxcL = PaL->Kx - Uc[BXn]*scrhL;     // Eq [47] 
  vxcR = PaR->Kx - Uc[BXn]*scrhR;

  vycL = PaL->Ky - Uc[BXt]*scrhL;
  vycR = PaR->Ky - Uc[BXt]*scrhR;

  vzcL = PaL->Kz - Uc[BXb]*scrhL;
  vzcR = PaR->Kz - Uc[BXb]*scrhR;

  vxc = 0.5*(vxcL + vxcR);
  vyc = 0.5*(vycL + vycR);
  vzc = 0.5*(vzcL + vzcR);

  if (vxc > 0.0) {
    HLLD_GetAState<dir> (PaL, p, Bx, grid);
    ua  = PaL->u;
    Sa  = PaL->Sa;
    vxa = PaL->vx;
  } else {
    HLLD_GetAState<dir> (PaR, p, Bx, grid);
    ua  = PaR->u;
    Sa  = PaR->Sa;
    vxa = PaR->vx;
  }

  vBc = vxc*Uc[BXn] + vyc*Uc[BXt] + vzc*Uc[BXb];

  Uc[RHO] = ua[RHO]*(Sa - vxa)/(Sa - vxc);                      // Eq [32] 
  Uc[ENG] = (Sa*ua[ENG] - ua[MXn] + p*vxc - vBc*Bx)/(Sa - vxc); // Eq [33] 

  Uc[MXn] = (Uc[ENG] + p)*vxc - vBc*Bx;                         // Eq [34] 
  Uc[MXt] = (Uc[ENG] + p)*vyc - vBc*Uc[BXt];
  Uc[MXb] = (Uc[ENG] + p)*vzc - vBc*Uc[BXb];
}

#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
double HLLD_TotalPressure (double *v)
/*!
  Compute total pressure
*/
////////////////////////////////////////////////////////////////////////
{
  double vel2, Bmag2, vB;
  double pt;

  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  vel2  = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  vB    = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  pt    = v[PRS] + 0.5*(Bmag2*(1.0 - vel2) + vB*vB);
  return(pt);
}
#pragma omp end declare target
