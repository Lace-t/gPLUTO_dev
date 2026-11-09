////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLLC Riemann solver for MHD.

  Solve the Riemann problem for the adiabatic  MHD equations using
  a modified version of the  HLLC Riemann solver of Li (2005).
  The isothermal version has not been implemented yet.

  Our formulation differs from Li's original solver in the way
  transverse momenta are computed.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b References:
    - "An HLLC Riemann Solver for MHD", S. Li, JCP (2000) 203, 344

  \authors A. Mignone  (andrea.mignone@unito.it)   \n
           M. Rossazza (marco.rossazza@unito.it)   \n
           A. Suriano  (alessio.suriano@unito.it)  
  \date    Apr 1, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void HLLC_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void HLLC_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if     (dir == IDIR) HLLC_Solver<IDIR>(data, Dts, grid, box);
  else if(dir == JDIR) HLLC_Solver<JDIR>(data, Dts, grid, box);
  else if(dir == KDIR) HLLC_Solver<KDIR>(data, Dts, grid, box);
}

#if HAVE_ENERGY
////////////////////////////////////////////////////////////////////////
/*
 * 
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLLC_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int      nv, i, k, j;
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

  #if BACKGROUND_FIELD == YES
  QuitPLUTO ("! HLLC_Solver(): background field splitting not allowed\n",1);
  #endif

  #if (defined STAGGERED_MHD) 
  #if (CT_EMF_AVERAGE == UCT_GFORCE) || (CT_EMF_AVERAGE == UCT_HLLD)
  QuitPLUTO ("! HLLC_Solver(): not compatible with UCT_GFORCE or UCT_HLLD \n", 1);
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
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts) \
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
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
    Ary1D press(&data->sweep->press[offset],ni);

    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    int switch_to_hll = 0;

    double Uhll[NVAR];
    double Fhll[NVAR];
    double usl[NVAR];
    double usr[NVAR];

    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, pl;
    double vr[NVAR], ur[NVAR], fr[NVAR], Sr, pr;
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

    // 2b. Compute sound & HLL speeds
      
    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    // 2c. Get max propagation speed for dt comp.

    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);

    cmax(i) = MAX(fabs(Sl), fabs(Sr));    
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    // 2d. Compute Flux if SL > 0 or SR < 0
      
    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    if (Sl >= 0.0){
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        flux(i,nv) = fl[nv];
      }
      press(i) = pl;

    }else if (Sr <= 0.0){
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        flux(i,nv) = fr[nv];
      }
      press(i) = pr;

    }else{

      // 2e. Compute Flux when  SL < 0 < SR

      double scrh = 1.0/(Sr - Sl);

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        Uhll[nv]  = Sr*ur[nv] - Sl*ul[nv] + fl[nv] - fr[nv];
        Uhll[nv] *= scrh;

        Fhll[nv]  =   Sl*Sr*(ur[nv] - ul[nv])
                    + Sr*fl[nv] - Sl*fr[nv];
        Fhll[nv] *= scrh;
      }
      // Uhll[MXn] += (pl - pr)*scrh;
      // Fhll[MXn] += (Sr*pl - Sl*pr)*scrh;

      #if HLL_SWITCH == YES
      if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
        switch_to_hll = 1;
      }
      #endif
      if (!switch_to_hll){  // Use HLLC 

        double Bxs = Uhll[BXn];
        double Bys = Uhll[BXt];
        double Bzs = Uhll[BXb];

        double vxL = vl[VXn];
        double vxR = vr[VXn];

        // ---- normal velocity vx  ----  

        double vxs = Uhll[MXn]/Uhll[RHO];

        // 2e1. Compute states

        {
        double vBL = vl[VX1]*vl[BX1] + vl[VX2]*vl[BX2] + vl[VX3]*vl[BX3];
        double vBR = vr[VX1]*vr[BX1] + vr[VX2]*vr[BX2] + vr[VX3]*vr[BX3];
        double vBs = Uhll[BX1]*Uhll[MX1]
                  + Uhll[BX2]*Uhll[MX2]
                  + Uhll[BX3]*Uhll[MX3];

        vBs /= Uhll[RHO];

        usl[RHO] = ul[RHO]*(Sl - vxL)/(Sl - vxs);
        usr[RHO] = ur[RHO]*(Sr - vxR)/(Sr - vxs);

        double ps  = Fhll[MXn] + Bxs*Bxs - Fhll[RHO]*vxs;

        usl[ENG] = (ul[ENG]*(Sl - vxL) +
                    ps*vxs - pl*vxL - Bxs*vBs + vl[BXn]*vBL)/(Sl - vxs);
        usr[ENG] = (ur[ENG]*(Sr - vxR) +
                    ps*vxs - pr*vxR - Bxs*vBs + vr[BXn]*vBR)/(Sr - vxs);
        }

        usl[MXn] = usl[RHO]*vxs;
        usr[MXn] = usr[RHO]*vxs;

        usl[MXt] =   (ul[MXt]*(Sl - vxL)
                  - (Bxs*Bys - vl[BXn]*vl[BXt]))/(Sl - vxs);
        usr[MXt] =   (ur[MXt]*(Sr - vxR)
                  - (Bxs*Bys - vr[BXn]*vr[BXt]))/(Sr - vxs);

        usl[MXb] =   (ul[MXb]*(Sl - vxL)
                  - (Bxs*Bzs - vl[BXn]*vl[BXb]))/(Sl - vxs);
        usr[MXb] =   (ur[MXb]*(Sr - vxR)
                  - (Bxs*Bzs - vr[BXn]*vr[BXb]))/(Sr - vxs);

        usl[BXn] = usr[BXn] = Bxs;
        usl[BXt] = usr[BXt] = Bys;
        usl[BXb] = usr[BXb] = Bzs;

        #ifdef GLM_MHD
        usl[PSI_GLM] = usr[PSI_GLM] = vl[PSI_GLM];
        #endif

        // 2e2. Compute flux

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
      }else{ // Use HLL 

        scrh = 1.0/(Sr - Sl);
        #pragma acc loop seq
        for (nv = 0; nv < NFLX; nv++) flux(i,nv) = Fhll[nv];

        press(i) = (Sr*pl - Sl*pr)*scrh;
      } // if (!switch_to_hll) 
    }  // end if (SL < 0 SR > 0) 

    // 2f. Compute Powell's / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
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
      if (isnan(flux(i,nv))) {
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
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk);
    #endif
    #endif

    // 2l. Conservative Energy correction step [Experimental]

    #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
    double EgL, EgR;
    double fEgL = HD_EnergyFlux(vl, &EgL, dir, run_config);
    double fEgR = HD_EnergyFlux(vr, &EgR, dir, run_config);
    double scrh = 1.0 / (bR - bL);

    flux(i,ENG) = (bR*fEgL - bL*fEgR + bL*bR*(EgR - EgL))*scrh; 
    #elif CT_ENERGY_CORRECTION == 5
    flux(i,ENG) = (bR*fEgL - bL*fEgR)*scrh; 
    #endif

  }}} // End loop (i,j,k)

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
#elif EOS == ISOTHERMAL

////////////////////////////////////////////////////////////////////////
template<int dir>
void HLLC_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  QuitPLUTO("from HLLC_Solver() not implemented for Isothermal EOS",1);
}

#endif // end #if on EOS 
