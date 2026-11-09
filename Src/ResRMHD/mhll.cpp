////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHLL Riemann solver for resistive RMHD.

  Solve the Riemann problem for the resistive relativistic MHD
  using the Maxwell + HLL solver.
   
  \b Reference:
    - "A constrained transport method for the solution of the
       resistive relativistic MHD equations",
      Mignone et al., MNRAS (2019) 486, 4252.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Dec 18, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#define BX_MIN  1.e-6

template<int dir>
void MHLL_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
void MHLL_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  if      (dir == IDIR) MHLL_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) MHLL_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) MHLL_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
template<int dir>
void MHLL_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
/*
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
  double glm_ch  = run_config.glm_ch;
  double glm_ch2 = glm_ch*glm_ch;
  #endif
  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  ARRAY_STRIDES (grid, dir, ni, nj);

  constexpr int VXn  = VX1+dir;
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, EX1, EXn, EXt, EXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  // --------------------------------------------------------
  // 0. Compatibility check
  // --------------------------------------------------------

  #if   (defined STAGGERED_MHD) \
     && (CT_EMF_AVERAGE != CT_MAXWELL) && (CT_EMF_AVERAGE != ARITHMETIC)
  QuitPLUTO("! LF Solver(): EMF_AVERAGE not supported", 1);
  #endif

  // --------------------------------------------------------
  // 1. Loop over domain box and solve Riemann problem
  // in the direction given by dir.
  // --------------------------------------------------------
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts)\
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

      // ------------------------------------------------------
      // 2. Loop over normal direction
      // ------------------------------------------------------

  //#pragma acc loop private(WtL[:NVAR], WtR[:NVAR], HtL[:NVAR], \
                             HtR[:NVAR], Wt_hll[:NVAR], Ft_hll[:NVAR])
    
   



        // ------------------------------------------------------
        // 1a. Compute offsets and assign local pointers
        // ------------------------------------------------------

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary1D cs2(&data->sweep->cs2[offset],ni);
    Ary1D press(&data->sweep->press[offset],ni);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    double  WtL[NVAR], HtL[NVAR];
    double  WtR[NVAR], HtR[NVAR];
    double  Wt_hll[NVAR], Ft_hll[NVAR];
      int ti, tj, tk;
      double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl, a2L, hL;
      double ur[NVAR], vr[NVAR], fr[NVAR], Sr, pr, a2R, hR;
      State stateL, stateR;

      if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
      if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
      if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

      for (nv = 0; nv < NVAR; nv++) {
        vl[nv] = vL(i,nv);
        vr[nv] = vR(i,nv);
      }

      data->err(ti,tj,tk) = 0;   // reset error 

      // ----------------------------------------------------
      // 2a. GLM pre-conditioner
      // ----------------------------------------------------

      #ifdef GLM_MHD
      GLM_Solve (vl, vr, run_config, dir);
      #endif
      PrimToCons (vl, ul, run_config);
      PrimToCons (vr, ur, run_config);

      // ----------------------------------------------------
      // 2b. Compute sound speed, cons variables, fluxes, 
      // hydro HLL speeds
      // ----------------------------------------------------

      STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
      STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

      SoundSpeed2 (stateL, run_config);
      SoundSpeed2 (stateR, run_config);

      a2L = stateL.a2; hL = stateL.h;
      a2R = stateR.a2; hR = stateR.h;
      stateL.eta = stateR.eta = Resistive_eta (stateL, run_config);

      Flux<dir> (stateL, fl, pl, run_config);
      Flux<dir> (stateR, fr, pr, run_config);

      HLL_Speed_RHD<dir> (vl, a2L, vr, a2R, &Sl, &Sr);

      // ----------------------------------------------------
      // 2c. Get max propagation speed & sound speed squared
      // ----------------------------------------------------

      cmax(i) = 1.0;
      {
        double bL = MIN(Sl, 0.0);
        double bR = MAX(0.0, Sr);
        cs2(i) = (bR*a2R - bL*a2L)/(bR - bL);
      }

      // ------------------------------------------------------
      // 2d. Define hydro states W in the tL, tR regions
      // ------------------------------------------------------

      NFLX_LOOP(nv){
        WtL[nv] = 0.0;
        WtR[nv] = 0.0;
        HtL[nv] = 0.0;
      }

      WtL[RHO] = ul[RHO];
      WtR[RHO] = ur[RHO];

      {
      double wL = hL*vl[RHO];
      double gL = ul[RHO]/vl[RHO];
      double gL2 = gL*gL;
      double wR = hR*vr[RHO];
      double gR = ur[RHO]/vr[RHO];
      double gR2 = gR*gR;

      WtL[MX1] = wL*gL2*vl[VX1];
      WtL[MX2] = wL*gL2*vl[VX2];
      WtL[MX3] = wL*gL2*vl[VX3];

      WtR[MX1] = wR*gR2*vr[VX1];
      WtR[MX2] = wR*gR2*vr[VX2];
      WtR[MX3] = wR*gR2*vr[VX3];

      WtL[ENG] = wL*gL2 - vl[PRS];
      WtR[ENG] = wR*gR2 - vr[PRS];
      }

      // ------------------------------------------------------
      // 2e. Compute hydro fluxes in tL, tR regions
      // ------------------------------------------------------

      HtL[RHO] = WtL[RHO]*vl[VXn];
      HtR[RHO] = WtR[RHO]*vr[VXn];

      HtL[MX1] = WtL[MX1]*vl[VXn];
      HtL[MX2] = WtL[MX2]*vl[VXn];
      HtL[MX3] = WtL[MX3]*vl[VXn];

      HtR[MX1] = WtR[MX1]*vr[VXn];
      HtR[MX2] = WtR[MX2]*vr[VXn];
      HtR[MX3] = WtR[MX3]*vr[VXn];

      HtL[MXn] += vl[PRS];
      HtR[MXn] += vr[PRS];

      HtL[ENG] = WtL[MXn];
      HtR[ENG] = WtR[MXn];

      // ----------------------------------------------------
      // 2f. Define HLL states and fluxes in tilde region
      // ----------------------------------------------------

      {
      double scrh = 1.0/(Sr - Sl);
      NFLX_LOOP(nv){
        Wt_hll[nv]  = Sr*WtR[nv] - Sl*WtL[nv] + HtL[nv] - HtR[nv];
        Wt_hll[nv] *= scrh;

        Ft_hll[nv]  =   Sr*Sl*(WtR[nv] - WtL[nv])
                      + Sr*HtL[nv] - Sl*HtR[nv];
        Ft_hll[nv] *= scrh;
      }
      }

      // ----------------------------------------------------
      // 2g. Compute hydro flux
      // ----------------------------------------------------

      if (Sl > 0.0){
        NFLX_LOOP(nv) flux(i,nv) = HtL[nv];
      }else if (Sr < 0.0){
        NFLX_LOOP(nv) flux(i,nv) = HtR[nv];
      }else{
        NFLX_LOOP(nv) flux(i,nv) = Ft_hll[nv];
      }
      press(i) = 0.0;

      // 2h. Add EM contribution to obtain total flux

      {
      double Bxt = 0.5*(ul[BXn] + ur[BXn]);
      double Byt = 0.5*(ul[BXt] + ur[BXt]) + 0.5*(ur[EXb] - ul[EXb]);
      double Bzt = 0.5*(ul[BXb] + ur[BXb]) - 0.5*(ur[EXt] - ul[EXt]);

      double Ext = 0.5*(ul[EXn] + ur[EXn]);
      double Eyt = 0.5*(ul[EXt] + ur[EXt]) - 0.5*(ur[BXb] - ul[BXb]);
      double Ezt = 0.5*(ul[EXb] + ur[EXb]) + 0.5*(ur[BXt] - ul[BXt]);

      double Sxt = Eyt*Bzt - Ezt*Byt;
      // Syt = Ezt*Bxt - Ext*Bzt;
      // Szt = Ext*Byt - Eyt*Bxt;
      double Pem_t =   0.5*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
              + 0.5*(Ext*Ext + Eyt*Eyt + Ezt*Ezt);

      flux(i,MXn) += - Ext*Ext - Bxt*Bxt + Pem_t;
      flux(i,MXt) += - Eyt*Ext - Byt*Bxt;
      flux(i,MXb) += - Ezt*Ext - Bzt*Bxt;
      flux(i,ENG) +=   Sxt;

      flux(i,EXn) =  0.0;
      flux(i,EXt) =  Bzt;
      flux(i,EXb) = -Byt;

      flux(i,BXn) =  0.0;
      flux(i,BXt) = -Ezt;
      flux(i,BXb) =  Eyt;
      }

      #if DIVB_CONTROL == DIV_CLEANING
      flux(i,BXn)     = 0.5*(vl[PSI_GLM] + vr[PSI_GLM]);
      flux(i,PSI_GLM) = 0.5*(vl[BXn] + vr[BXn])*glm_ch2;
      #endif

      #if DIVE_CONTROL == DIV_CLEANING
      flux(i,PHI_GLM) = 0.5*(vl[EXn] + vr[EXn])*glm_ch2;
      flux(i,EXn)     = 0.5*(vl[PHI_GLM] + vr[PHI_GLM]);
      flux(i,CRG)     = 0.5*(fl[CRG] + fr[CRG]) - 0.5*(ur[CRG] - ul[CRG]);
      #endif

      // ----------------------------------------------------
      // 2i. Check for NaNs in flux
      // ----------------------------------------------------
 
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        if (isnan(flux(i,nv))) {
          ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
          KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
        }
      }

      // ----------------------------------------------------
      // 2j. Compute tracer fluxes
      // ----------------------------------------------------

      #if NSCL > 0
      if (flux(i,RHO) > 0.0){
        #pragma acc loop seq
        NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
      }else{
        #pragma acc loop seq
        NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
      }
      #endif

      // ------------------------------------------------------
      // 2k. Store emf-flux components for CT
      // ------------------------------------------------------

      #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
      #if CT_EMF_AVERAGE == ARITHMETIC 
      CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
      #elif CT_EMF_AVERAGE == CT_MAXWELL
      CT_Store_Maxwell<dir>(data, fl, fr, flux, i, ti, tj, tk);
      #endif
      #endif

    }  // End loop on normal direction 
  }} // End loop on transverse direction 

  // --------------------------------------------------------
  // 3. Error Check / Quit condition
  // --------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("MHLL_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from MHLL_Solver()",err_kernel);
  #endif

}
