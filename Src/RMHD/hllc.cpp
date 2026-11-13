////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Implement the HLLC Riemann solver for relativistic MHD.

  Solve the Riemann problem for the relativistic MHD (RMHD) equations
  using the HLLC solver of Mignone & Bodo (2006).

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
    - "An HLLC Riemann solver for relativistic flows - II.
       Magnetohydrodynamics", Mignone and Bodo, MNRAS (2006) 368, 1040.

  \authors A. Mignone  (andrea.mignone@unito.it)  \n
           M. Rossazza (marco.rossazza@unito.it)
           M. Bugli    (matteo.bugli@cea.fr)
  \date    Jun 09, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define BX_MIN  1.e-6

template<int dir>
void HLLC_Solver (Data *, timeStep *Dts, Grid *grid, RBox *box);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void HLLC_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) HLLC_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) HLLC_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) HLLC_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLLC_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
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
  QuitPLUTO ("! HLLC_Solver: not compatible with UCT_GFORCE or UCT_HLLD \n", 1);
  #endif
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  constexpr int VXn = VX1+dir;
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
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

    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    double Uhll[NVAR];
    double Fhll[NVAR];
    double usl[NVAR];
    double usr[NVAR];

    int  switch_to_hll = 0;
    int ti, tj, tk;
    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl, a2L;
    double vr[NVAR], ur[NVAR], fr[NVAR], Sr, pr, a2R;
    State stateL, stateR;

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

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

    // 2b. compute sound speed, cons variables, fluxes
    //     and HLL signal speeds.

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    a2L = stateL.a2;
    a2R = stateR.a2;
    #if PHYSICS != GRMHD
    HLL_Speed<dir> (stateL, stateR, Sl, Sr);
    #else
    HLL_Speed<dir> (stateL, stateR, Sl, Sr, &G);
    #endif

    // 2c. Get max propagation speed & 
    //     sound speed squared

    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);
    cmax(i) = MAX(fabs(Sl), fabs(Sr));
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    // 2d. Compute HLL state and flux

    #if PHYSICS != GRMHD
    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);
    #else
    Flux<dir> (stateL, fl, pl, run_config, &G);
    Flux<dir> (stateR, fr, pr, run_config, &G);
    #endif

    double scrh = 1.0/(Sr - Sl);

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) {
      Uhll[nv] = Sr*ur[nv] - Sl*ul[nv] + fl[nv] - fr[nv];
      Uhll[nv] *= scrh;

      Fhll[nv]  = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fl[nv] - Sl*fr[nv];
      Fhll[nv] *= scrh;
    }
    // Uhll[MXn] += (pl - pr)*scrh;
    // Fhll[MXn] += (Sr*pl - Sl*pr)*scrh;

    // 2e. Solve Riemann problem with HLLC

    if (Sl >= 0.0){      // Region L 

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fl[nv];
      press(i) = pl;

    }else if (Sr <= 0.0){   // Region R 

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fr[nv];
      press(i) = pr;

    }else{     // Region * 

      // 2e1. Region SL < 0 < SR   

      #if HLL_SWITCH == YES
      if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
        switch_to_hll = 1;
      }
      #endif

      // A failed HLLC construction must not reach sqrt or star-state divisions.
      // Keep the existing HLL output branch (including its pressure convention).
      double vxs = 0.0;
      if (!switch_to_hll) {
      double vxl = vl[VXn];
      double vxr = vr[VXn];

      double Bx  = Uhll[BXn];
      double Bys = Uhll[BXt];
      double Bzs = Uhll[BXb];

      // ! Attention: derivation below is done without Reduced energy !

      #if RMHD_REDUCED_ENERGY == YES
      Uhll[ENG] += Uhll[RHO];
      Fhll[ENG] += Fhll[RHO];
      #endif

      // 2e2. Obtain coefficients of quadratic 
      //      equation and obtain v* (=vxs)

      double a,b,c, BtFBt, Bt2, FBt2;
      if (fabs(Bx) < BX_MIN) {  // Bx -> 0 Limit, Eq. [47]  
        BtFBt = Bt2 = FBt2 = 0.0;
        a  = Fhll[ENG];
        b  = - (Fhll[MXn] + Uhll[ENG]);
        c  = Uhll[MXn];
      }else{
        BtFBt = Uhll[BXt]*Fhll[BXt] + Uhll[BXb]*Fhll[BXb];
        Bt2   = Uhll[BXt]*Uhll[BXt] + Uhll[BXb]*Uhll[BXb];
        FBt2  = Fhll[BXt]*Fhll[BXt] + Fhll[BXb]*Fhll[BXb];
        a  = Fhll[ENG] - BtFBt;
        b  = Bt2 + FBt2 - (Fhll[MXn] + Uhll[ENG]);
        c  = Uhll[MXn] - BtFBt;
      }

      // 2e3. Solve quadratic, get v*

      // Preserve the original root on regular interfaces. A negative
      // discriminant (even a tiny one) safely selects HLL rather than
      // inventing a real root. Degenerate/nonfinite arithmetic does likewise.
      const double b2 = b*b;
      if (!isfinite(a) || !isfinite(b) || !isfinite(c) ||
          !isfinite(b2) || !(b2 > 0.0)) {
        switch_to_hll = 1;
      } else {
        const double discriminant = 1.0 - 4.0*a*c/b2;
        if (!isfinite(discriminant) || discriminant < 0.0) {
          switch_to_hll = 1;
        } else {
          scrh = 1.0 + sqrt(discriminant);
          const double denominator = b*scrh;
          if (!isfinite(denominator) || denominator == 0.0) {
            switch_to_hll = 1;
          } else {
            vxs = -2.0*c/denominator;
            if (!isfinite(vxs) || !(fabs(vxs) < 1.0) ||
                !(Sl < vxs && vxs < Sr)) switch_to_hll = 1;
          }
        }
      }

      if (!switch_to_hll) {

      // 2e4. Compute *L and *R states

      if (fabs(Bx) < BX_MIN) {

        // Degenerate case: the values of 
        // (vy,vz) are irrelevant in this case

        double alpha_l = (Sl - vxl)/(Sl - vxs);
        double alpha_r = (Sr - vxr)/(Sr - vxs);
        double ps;

        usl[RHO] = ul[RHO]*alpha_l;
        usr[RHO] = ur[RHO]*alpha_r;

        ps  = Fhll[MXn] - Fhll[ENG]*vxs;

        usl[ENG] = (Sl*ul[ENG] - fl[ENG] + ps*vxs)/(Sl - vxs);
        usr[ENG] = (Sr*ur[ENG] - fr[ENG] + ps*vxs)/(Sr - vxs);

        usl[MXn] = (usl[ENG] + ps)*vxs;
        usr[MXn] = (usr[ENG] + ps)*vxs;

        usl[MXt] = ul[MXt]*alpha_l;
        usr[MXt] = ur[MXt]*alpha_r;
        usl[MXb] = ul[MXb]*alpha_l;
        usr[MXb] = ur[MXb]*alpha_r;

        #if RMHD_REDUCED_ENERGY == YES
        usl[MXn] += usl[RHO]*vxs;
        usr[MXn] += usr[RHO]*vxs;
        #endif

        usl[BXn] = Bx;
        usr[BXn] = Bx;
        usl[BXt] = ul[BXt]*alpha_l;
        usr[BXt] = ur[BXt]*alpha_r;
        usl[BXb] = ul[BXb]*alpha_l;
        usr[BXb] = ur[BXb]*alpha_r;

      }else{

        // Non-degenerate case

        double vys = (Bys*vxs - Fhll[BXt])/Bx;
        double vzs = (Bzs*vxs - Fhll[BXb])/Bx;

        double gammas_2 = vxs*vxs + vys*vys + vzs*vzs;
        gammas_2 = 1.0 - gammas_2;

        // Avoid superluminal speed by 
        // switching to HLL solver

        if (!isfinite(gammas_2) || gammas_2 <= 0.0) switch_to_hll = 1;

        if (!switch_to_hll) {

        // Complete solver

        double vBs = vxs*Bx + vys*Bys + vzs*Bzs;

        double alpha_l = (Sl - vxl)/(Sl - vxs);
        double alpha_r = (Sr - vxr)/(Sr - vxs);
        double ps;

        usl[RHO] = ul[RHO]*alpha_l;
        usr[RHO] = ur[RHO]*alpha_r;

        ps = (BtFBt - Fhll[ENG])*vxs + Bx*Bx + Fhll[MXn] - FBt2;

        usl[ENG] = (Sl*ul[ENG] - fl[ENG] + ps*vxs - vBs*Bx)/(Sl - vxs);
        usr[ENG] = (Sr*ur[ENG] - fr[ENG] + ps*vxs - vBs*Bx)/(Sr - vxs);

        usl[MXn] = (usl[ENG] + ps)*vxs - vBs*Bx;
        usr[MXn] = (usr[ENG] + ps)*vxs - vBs*Bx;

        usl[MXt] =  (Sl*ul[MXt] - fl[MXt] - Bx*(Bys*gammas_2 + vBs*vys))
                    /(Sl - vxs);
        usr[MXt] =  (Sr*ur[MXt] - fr[MXt] - Bx*(Bys*gammas_2 + vBs*vys))
                    /(Sr - vxs);
        usl[MXb] =  (Sl*ul[MXb] - fl[MXb] - Bx*(Bzs*gammas_2 + vBs*vzs))
                    /(Sl - vxs);
        usr[MXb] =  (Sr*ur[MXb] - fr[MXb] - Bx*(Bzs*gammas_2 + vBs*vzs))
                    /(Sr - vxs);

        #if RMHD_REDUCED_ENERGY == YES
        usl[MXn] += usl[RHO]*vxs;
        usr[MXn] += usr[RHO]*vxs;
        #endif

        usl[BXn] = usr[BXn] = Bx;
        usl[BXt] = usr[BXt] = Bys;
        usl[BXb] = usr[BXb] = Bzs;
        } // finite, subluminal star velocity
      }

      #ifdef GLM_MHD
      usl[PSI_GLM] = usr[PSI_GLM] = ul[PSI_GLM];
      #endif

      if (!switch_to_hll) {
        if (!(usl[RHO] > 0.0) || !(usr[RHO] > 0.0)) switch_to_hll = 1;
        #pragma acc loop seq
        for (int n = 0; n < NFLX; ++n) {
          if (!isfinite(usl[n]) || !isfinite(usr[n])) switch_to_hll = 1;
        }
      }
      } // valid contact speed
      } // HLLC requested

      // 2e5. Check consistency condition
/*
      NFLX_LOOP(nv) {
        scrh  = (vxs - Sl)*usl[nv] + (Sr - vxs)*usr[nv];
        scrh -= Sr*ur[i][nv] - Sl*ul[i][nv] +
                fL[i][nv] - fR[i][nv];
        if (fabs(scrh/(Sr-Sl)) > 1.e-5){
          //printf (" Consistency condition violated\n");
          //printf ("%d %d  %12.6e\n",i,nv, scrh);
        }
      }
*/
      // 2e6. Compute flux

      if (!switch_to_hll){  // Use HLLC flux 

        if (vxs > 0.0) {
          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++) {
            flux(i,nv) = fl[nv] + Sl*(usl[nv] - ul[nv]);
          }
          press(i) = pl;
        }else {
          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++) {
            flux(i,nv) = fr[nv] + Sr*(usr[nv] - ur[nv]);
          }
          press(i) = pr;
        }

        // Finite star states can still overflow when forming the flux.
        #pragma acc loop seq
        for (int n = 0; n < NFLX; ++n) {
          if (!isfinite(flux(i,n))) switch_to_hll = 1;
        }
        if (!isfinite(press(i))) switch_to_hll = 1;
      }
      if (switch_to_hll) {   // Use HLL flux 

        // Note: need to re-calculate F(hll) since the previous 
        // array Fhll[] may have been modified with RMHD_REDUCED_ENERGY 

        scrh = 1.0/(Sr - Sl);
        #pragma acc loop seq
        for (nv = 0; nv < NFLX; nv++) {
          flux(i,nv)  = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fl[nv] - Sl*fr[nv];
          flux(i,nv) *= scrh;
        }
        press(i) = (Sr*pl - Sl*pr)*scrh;
      }
    }

    // 2f. Compute Powell's interface state

    #if DIVB_CONTROL == EIGHT_WAVES
    //      Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    //      Bn(i) = (flux(i,RHO) > 0.0 ? vl[BXn]:vr[BXn]);
    scrh = 1.0/(bR - bL);
    Bn(i) = (bR*vl[BXn] - bL*vr[BXn])*scrh;
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    // 2g. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2h. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (!isfinite(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

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
  Where ("HLLC_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from HLLC_Solver()",err_kernel);
  #endif

}
