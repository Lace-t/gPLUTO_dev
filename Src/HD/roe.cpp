////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Roe Riemann solver for the HD equations.

  Solve the Riemann problem for the Euler equations of gasdynamics
  using the standard Roe solver with local characteristic decomposition.
  Eigenvectors are identical to the ones given in the book by Toro
  and can also be derived from the maple script "eigenv.maple"
  in Src/HD/.
  The solver can be used for adiabatic and isothermal hydrodynamics.

  The macro ROE_AVERAGE specifies how the averaging process is done:
    - ROE_AVERAGE == YES use Roe average (default);
    - ROE_AVERAGE == NO  use arithmetic average.

  The macro CHECK_ROE_MATRIX can be used to verify that the
  characteristic decomposition reproduces the Roe matrix.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2.
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
   -   "Riemann Solver and Numerical Methods for Fluid Dynamics"
        by E.F. Toro (Chapter 11)

  \authors A. Mignone  (andrea.mignone@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it)\n
           V. Berta    (vittoria.berta@unito.it)\n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define sqrt_1_2  (0.70710678118654752440)
#define ROE_AVERAGE       YES

template<int dir>
void Roe_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
void Roe_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  if      (dir == IDIR) Roe_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) Roe_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) Roe_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
  Solve the Riemann problem between L/R states using a
  Rusanov-Lax Friedrichs flux.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void Roe_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int      i, j, k, nv;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int      ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  #if EOS == IDEAL
  double gmm = data->run_config->gamma;
  double gmm1 = gmm - 1.0;
  double gmm1_inv  = 1.0/gmm1;
  #endif
  double delta     = 1.e-7;

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  RunConfig &run_config = *(data->run_config);

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);

  // Cycle on directional indexes to preserve sweep ordering

  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

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
    Ary1D press(&data->sweep->press[offset],ni);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    int switch_to_hll = 0;

    double vl[NVAR], ul[NVAR], fL[NVAR], Sl, pl;   
    double vr[NVAR], ur[NVAR], fR[NVAR], Sr, pr;
    double Rc[NVAR];
    double um[NVAR], dv[NVAR];
    double a, a2;

    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR)  {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR)  {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR)  {ti = j; tj = k; tk = i;}

    NFLUID_LOOP(nv) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    data->err(ti,tj,tk) = 0;   // reset error 

    double alambda, lambda, eta;

    // 2a. Compute sound speed, cons variables and fluxes.

    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    Flux<dir> (stateL, fL, pl);
    Flux<dir> (stateR, fR, pr);

    // 2b. Enable switch to HLL in proximity of strong shock

    #if HLL_SWITCH == YES
    if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
      switch_to_hll = 1;
    }
    #endif

    if(!switch_to_hll){ // Use ROE

      // ----------------------------------------
      // 2c. Compute fluxes using Roe solver
      // -----------------------------------------

      #pragma acc loop seq
      NFLX_LOOP(nv) dv[nv] = vr[nv] - vl[nv];

      #if ROE_AVERAGE == YES
      double s = sqrt(vr[RHO]/vl[RHO]);
      um[RHO]  = vl[RHO]*s;
      s        = 1.0/(1.0 + s);
      double c = 1.0 - s;

      um[VX1] = s*vl[VX1] + c*vr[VX1];
      um[VX2] = s*vl[VX2] + c*vr[VX2];
      um[VX3] = s*vl[VX3] + c*vr[VX3];

      #if EOS == IDEAL
      double vel2 = um[VX1]*um[VX1] + um[VX2]*um[VX2] + um[VX3]*um[VX3];

      double hl = 0.5*(vl[VX1]*vl[VX1] + vl[VX2]*vl[VX2] + vl[VX3]*vl[VX3]);
      hl += a2L*gmm1_inv;

      double hr = 0.5*(vr[VX1]*vr[VX1] + vr[VX2]*vr[VX2] + vr[VX3]*vr[VX3]);
      hr += a2R*gmm1_inv;

      double h = s*hl + c*hr;

      // 2c1. The following should be  equivalent to
      // 
      //      scrh = dv[VX1]*dv[VX1] + dv[VX2]*dv[VX2] + dv[VX3]*dv[VX3];
      //      a2 = s*a2L + c*a2R + 0.5*gmm1*s*c*scrh;
      //
      //      and therefore always positive.

      a2 = gmm1*(h - 0.5*vel2);
      a  = sqrt(a2);
      #endif // EOS == IDEAL 
      #else
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) um[nv] = 0.5*(vl[nv] + vr[nv]);
      #if EOS == IDEAL
      a2   = gmm*um[PRS]/um[RHO];
      a    = sqrt(a2);

      double vel2 = um[VX1]*um[VX1] + um[VX2]*um[VX2] + um[VX3]*um[VX3];
      double h    = 0.5*vel2 + a2/gmm1;
      #endif // EOS == IDEAL 
      #endif // ROE_AVERAGE == YES/NO 

      #if EOS == ISOTHERMAL
      a2 = 0.5*(a2L + a2R);
      a  = sqrt(a2);
      #endif

      // 2c2. Define non-zero components of conservative eigenvectors Rc,
      //      eigenvalues (lambda), wave strenght eta = L.du and
      //      Roe numerical flux

      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) = 0.5*(fL[nv] + fR[nv]);

      press(i) = 0.5*(pl + pr);

      // Lambda = u - c_s

      lambda = um[VXn] - a;

      #if EOS == IDEAL
      eta = 0.5/a2*(dv[PRS] - dv[VXn]*um[RHO]*a);
      #elif EOS == ISOTHERMAL
      eta = 0.5*(dv[RHO] - um[RHO]*dv[VXn]/a);
      #endif

      Rc[RHO] = 1.0;
      Rc[MXn] = um[VXn] - a;
      Rc[MXt] = um[VXt];
      Rc[MXb] = um[VXb];
      #if EOS == IDEAL
      Rc[ENG] = h - um[VXn]*a;
      #endif

      Sl = lambda;
      alambda = fabs(lambda);
      if (alambda <= delta) alambda = 0.5*lambda*lambda/delta + 0.5*delta;
      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) -= 0.5*alambda*eta*Rc[nv];

      // Lambda = u + c_s

      lambda = um[VXn] + a;
      #if EOS == IDEAL
      eta    = 0.5/a2*(dv[PRS] + dv[VXn]*um[RHO]*a);
      #elif EOS == ISOTHERMAL
      eta = 0.5*(dv[RHO] + um[RHO]*dv[VXn]/a);
      #endif

      Rc[RHO] = 1.0;
      Rc[MXn] = um[VXn] + a;
      Rc[MXt] = um[VXt];
      Rc[MXb] = um[VXb];
      #if EOS == IDEAL
      Rc[ENG] = h + um[VXn]*a;
      #endif

      Sr = lambda;
      alambda = fabs(lambda);
      if (alambda <= delta) alambda = 0.5*lambda*lambda/delta + 0.5*delta;
      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) -= 0.5*alambda*eta*Rc[nv];

      // Lambda = u  (1)

      #if EOS == IDEAL
      lambda  = um[VXn];
      eta     = dv[RHO] - dv[PRS]/a2;
      Rc[RHO] = 1.0;
      Rc[MX1] = um[VX1];
      Rc[MX2] = um[VX2];
      Rc[MX3] = um[VX3];
      Rc[ENG] = 0.5*vel2;
      #endif

      #if EOS == IDEAL
      alambda = fabs(lambda);
      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) -= 0.5*alambda*eta*Rc[nv];
      #endif

      // Lambda = u  (2)

      lambda  = um[VXn];
      eta     = um[RHO]*dv[VXt];
      Rc[RHO] = 0.0;
      Rc[MXn] = 0.0;
      Rc[MXt] = 1.0;
      Rc[MXb] = 0.0;
      #if EOS == IDEAL
      Rc[ENG] = um[VXt];
      #endif

      alambda = fabs(lambda);
      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) -= 0.5*alambda*eta*Rc[nv];

      // Lambda = u  (3)

      lambda  = um[VXn];
      eta     = um[RHO]*dv[VXb];
      Rc[RHO] = 0.0;
      Rc[MXn] = 0.0;
      Rc[MXt] = 0.0;
      Rc[MXb] = 1.0;
      #if EOS == IDEAL
      Rc[ENG] = um[VXb];
      #endif

      alambda = fabs(lambda);
      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) -= 0.5*alambda*eta*Rc[nv];

      // 2c3. Maximum eigenvalue

      cmax(i) = fabs(um[VXn]) + a;
      cs2(i)  = a2;

      #if DIMENSIONS > 1

      // 2c4. Use the HLL flux function if the interface
      //      lies within a strong shock.
      //      The effect of this switch is visible
      //      in the Mach reflection test.

      double err = vr[VXn] - vl[VXn];
      #if EOS == IDEAL
      double scrh  = fabs(vl[PRS] - vr[PRS]);
      scrh /= MIN(vl[PRS],vr[PRS]);
      #elif EOS == ISOTHERMAL
      double scrh  = fabs(vl[RHO] - vr[RHO]);
      scrh /= MIN(vl[RHO],vr[RHO]);
      scrh *= a*a;
      #endif

      if (scrh > 0.5 && err < -1.0e-12){   // -- tunable parameter -- 
        double bmin = MIN(0.0, Sl);
        double bmax = MAX(0.0, Sr);
        double scrh1 = 1.0/(bmax - bmin);

        #pragma acc loop seq
        for (nv = 0; nv < NFLX; nv++){
          flux(i,nv)  = bmin*bmax*(ur[nv] - ul[nv])
                      + bmax*fL[nv] - bmin*fR[nv];
          flux(i,nv) *= scrh1;
        }
        press(i) = (bmax*pl - bmin*pr)*scrh1;
      }
      #endif  // DIMENSIONS > 1 

    } else { // Use HLL 

      // ----------------------------------------
      // 2d. Compute flux using HLL Solver
      // ----------------------------------------

      // HLL switching function as in Quirk (1994).
      // Since the problem is related to multidimensional
      // pathologies, it works in more than 1-D only.
      // Use the HLL flux function if the interface lies
      // within a strong shock.
      // The effect of this switch is visible in the
      // Mach reflection test.

      HLL_Speed<dir> (stateL, stateR, Sl, Sr);

      double scrh = MAX(fabs(Sl), fabs(Sr));
      cmax(i) = scrh;
      Sl  = MIN(0.0, Sl);
      Sr  = MAX(0.0, Sr);

      scrh = 1.0 / (Sr - Sl);
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        flux(i,nv)  = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fL[nv] - Sl*fR[nv];
        flux(i,nv) *= scrh;
      }
      press(i) = (Sr*pl - Sl*pr)*scrh;

    } // if(!switch_to_hll) 

    // 2e. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

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

  }}} // End Domain loop (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("Roe_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from Roe_Solver()",err_kernel);
  #endif
}
#undef ROE_AVERAGE
