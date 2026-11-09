////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLLD Riemann solver for the classical MHD equations.

  Solve the Riemann problem for the classical MHD equations with
  an ideal or isothermal equation of state.
  For the ideal case, we implement the four-states (of five-wave)
  HLLD solver of   Miyoshi & Kusano (2005).
  For the isothermal case, we use the three-states approximation
  described by Mignone (2007).

  The macro VERIFY_CONSISTENCY_CONDITION can be turned to YES to verify
  the correctness of the implementation.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b References:
    - "A muLti-state HLL approximate Riemann Solver for ideal MHD",
       Miyoshi, T., Kusano, K., JCP 2005
    - "A simple and accurate Riemann solver for isothermal MHD",
       Mignone, JCP (2007) 225, 1427.

  \authors A. Mignone  (andrea.mignone@unito.it)   \n  
           M. Rossazza (marco.rossazza@unito.it)   \n
           S. Truzzi   (stefano.truzzi@unito.it)   \n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////

#include"gpluto.hpp"

#define VERIFY_CONSISTENCY_CONDITION  NO

template<int dir>
void HLLD_Solver (Data *, timeStep *, Grid *, RBox *);

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

#if EOS == IDEAL || EOS == PVTE_LAW
////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLLD_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int  nv, i, k , j;
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  int  ni, nj;
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
  #pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts)\
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
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    int switch_to_hll = 0;
   
    double vl[NVAR], ul[NVAR], fl[NVAR], usl[NVAR], Sl, ptl, sqrL;
    double vr[NVAR], ur[NVAR], fr[NVAR], usr[NVAR], Sr, ptr, sqrR;
    double bgf[3];
    double Uhll[NVAR];

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

    // 2b. Compute background field, sound 
    //     and HLL speeds
 
    #if BACKGROUND_FIELD == YES
    BackgroundField (stateL, run_config); // stateL->bgf = stateR->bgf
    #endif

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);

    cmax(i) = MAX(fabs(Sl), fabs(Sr));    
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    // 2c. Compute fluxes and total pressure in the * region

    Flux<dir> (stateL, fl, ptl, run_config);
    Flux<dir> (stateR, fr, ptr, run_config);

    #if BACKGROUND_FIELD == YES
    double B0n = bgf[BXn-BX1];
    double B0t = bgf[BXt-BX1];
    double B0b = bgf[BXb-BX1];
    #endif

    double dul, dur, Bx, Bx1, sBx, pts, SM;
    {
    double scrh = 1.0/(Sr - Sl);
    Bx1  = Bx = (Sr*vr[BXn] - Sl*vl[BXn])*scrh;

    #if BACKGROUND_FIELD == YES
    Bx += B0n;   // Bx will be now the (normal) total field 
    #endif

    sBx  = (Bx > 0.0 ? 1.0 : -1.0);

    dul  = Sl - vl[VXn];
    dur  = Sr - vr[VXn];

    scrh = 1.0/(dur*ur[RHO] - dul*ul[RHO]);
    SM   = (dur*ur[MXn] - dul*ul[MXn] - ptr + ptl)*scrh;

    pts  =   dur*ur[RHO]*ptl - dul*ul[RHO]*ptr
           + vl[RHO]*vr[RHO]*dur*dul*(vr[VXn]- vl[VXn]);
    pts *= scrh;
    }

    // 2d. Set Alfven and contact waves velocities

    usl[RHO] = ul[RHO]*dul/(Sl - SM);
    usr[RHO] = ur[RHO]*dur/(Sr - SM);

    if (usl[RHO] < 0.0 || usr[RHO] < 0.0)  switch_to_hll = 1;

    sqrL = sqrt(usl[RHO]);
    sqrR = sqrt(usr[RHO]);

    double S1L = SM - fabs(Bx)/sqrL;
    double S1R = SM + fabs(Bx)/sqrR;

    // 2e. Compute Flux if SL > 0 or SR < 0

    if (Sl >= 0.0){                     // --- Region L ---
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fl[nv];
      press(i) = ptl;

    }else if (Sr <= 0.0) {              // --- Region R ---
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fr[nv];
      press(i) = ptr;

    } else {

      // ----------------------------------
      // 2f. Compute Flux when  SL < 0 < SR
      // ----------------------------------

      #if HLL_SWITCH == YES
      if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
        switch_to_hll = 1;
      }
      #endif

      if (!switch_to_hll){ // Use HLLD/HLLC 

        // 2f1. When S1L -> SL or S1R -> SR a degeneracy occurs.
        //      Although Miyoshi & Kusano say that no jump exists, we don't
        //      think this is actually true.
        //      Indeed, vy*, vz*, By*, Bz* cannot be solved independently.
        //      In this case we revert to the HLLC solver of Li (2005),  except
        //      for the term v.B in the region, which we compute in our own way.
        //      Note, that by comparing the expressions of Li (2005) and
        //      Miyoshi & Kusano (2005), the only change involves a
        //      re-definition of By* and Bz* in terms of By(HLL), Bz(HLL).
        
        int switch_to_hllc = 0;
        if ( (S1L - Sl) <  1.e-4*(SM - Sl) ) switch_to_hllc = 1;
        if ( (S1R - Sr) > -1.e-4*(Sr - SM) ) switch_to_hllc = 1;

        double scrhL, scrhR;

        // 2f2. Revert to hllc solver when |Sa-S| < tol

        if (!switch_to_hllc){ // Use HLLD 

          // 2f3. Compute states in the *L and *R regions.
          //      With bck field splitting, the intermediate
          //      mag. fields always contain the deviation B1.

          scrhL = (ul[RHO]*dul*dul - Bx*Bx)/(ul[RHO]*dul*(Sl - SM) - Bx*Bx);
          scrhR = (ur[RHO]*dur*dur - Bx*Bx)/(ur[RHO]*dur*(Sr - SM) - Bx*Bx);

          usl[BXn] = Bx1;
          usl[BXt] = ul[BXt]*scrhL;
          usl[BXb] = ul[BXb]*scrhL;
          #if BACKGROUND_FIELD == YES
          usl[BXt] += B0t*(scrhL - 1.0);    // Eq. [40] of         
          usl[BXb] += B0b*(scrhL - 1.0);    // Miyoshi etal (2010) 
          #endif

          usr[BXn] = Bx1;
          usr[BXt] = ur[BXt]*scrhR;
          usr[BXb] = ur[BXb]*scrhR;

          #if BACKGROUND_FIELD == YES
          usr[BXt] += B0t*(scrhR - 1.0);   // Eq. [40] of         
          usr[BXb] += B0b*(scrhR - 1.0);   // Miyoshi etal (2010) 
          #endif
        } else { // Use HLLC 
          double scrh = 1.0/(Sr - Sl);
          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++) {
            Uhll[nv]  = Sr*ur[nv] - Sl*ul[nv] + fl[nv] - fr[nv];
            Uhll[nv] *= scrh;
          }

          // In gpluto we add the press(i)ure term to Uhll[MXn]
          // that is not used here though                    

          usl[BXn] = usr[BXn] = Bx1;
          usl[BXt] = usr[BXt] = Uhll[BXt];
          usl[BXb] = usr[BXb] = Uhll[BXb];

          S1L = S1R = SM; // region ** should never be computed since 
                          // fluxes are given in terms of UL* and UR* 

        } // end if (switch_to_hllc) 

        scrhL = Bx/(ul[RHO]*dul);
        scrhR = Bx/(ur[RHO]*dur);

        double vsL = vl[VXt] - scrhL*(usl[BXt] - ul[BXt]);
        double vsR = vr[VXt] - scrhR*(usr[BXt] - ur[BXt]);

        double wsL = vl[VXb] - scrhL*(usl[BXb] - ul[BXb]);
        double wsR = vr[VXb] - scrhR*(usr[BXb] - ur[BXb]);

        usl[MXn] = usl[RHO]*SM;
        usr[MXn] = usr[RHO]*SM;

        usl[MXt] = usl[RHO]*vsL;
        usr[MXt] = usr[RHO]*vsR;

        usl[MXb] = usl[RHO]*wsL;
        usr[MXb] = usr[RHO]*wsR;

        scrhL     = vl[VXn]*Bx1 + vl[VXt]*ul[BXt] + vl[VXb]*ul[BXb];
        scrhL    -= SM*Bx1      +    vsL*usl[BXt] + wsL*usl[BXb];
        usl[ENG]  = dul*ul[ENG] - ptl*vl[VXn] + pts*SM + Bx*scrhL;
        usl[ENG] /= Sl - SM;

        scrhR     = vr[VXn]*Bx1 + vr[VXt]*ur[BXt] + vr[VXb]*ur[BXb];
        scrhR    -= SM*Bx1      +    vsR*usr[BXt] +    wsR*usr[BXb];
        usr[ENG]  = dur*ur[ENG] - ptr*vr[VXn] + pts*SM + Bx*scrhR;
        usr[ENG] /= Sr - SM;

        #ifdef GLM_MHD
        usl[PSI_GLM] = usr[PSI_GLM] = vl[PSI_GLM];
        #endif

        // 2f5. Compute flux when S1L > 0 or S1R < 0

        if (S1L >= 0.0){       //  --- Region L* --- 
          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++) {
            flux(i,nv) = fl[nv] + Sl*(usl[nv] - ul[nv]);
          }
          press(i) = ptl;

        }else if (S1R <= 0.0) {    //  --- Region R* ---
          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++) {
            flux(i,nv) = fr[nv] + Sr*(usr[nv] - ur[nv]);
          }
          press(i) = ptr;

        } else {   // -- This state exists only if B_x != 0 -- 

          // 2f6. Compute U**

          double ussl[NVAR], ussr[NVAR];

          ussl[RHO] = usl[RHO];
          ussr[RHO] = usr[RHO];

          double vss  = sqrL*vsL + sqrR*vsR + (usr[BXt] - usl[BXt])*sBx;
          vss /= sqrL + sqrR;

          double wss  = sqrL*wsL + sqrR*wsR + (usr[BXb] - usl[BXb])*sBx;
          wss /= sqrL + sqrR;

          ussl[MXn] = ussl[RHO]*SM;
          ussr[MXn] = ussr[RHO]*SM;

          ussl[MXt] = ussl[RHO]*vss;
          ussr[MXt] = ussr[RHO]*vss;

          ussl[MXb] = ussl[RHO]*wss;
          ussr[MXb] = ussr[RHO]*wss;

          ussl[BXn] = ussr[BXn] = Bx1;

          ussl[BXt]  = sqrL*usr[BXt] + sqrR*usl[BXt] + sqrL*sqrR*(vsR - vsL)*sBx;
          ussl[BXt] /= sqrL + sqrR;
          ussr[BXt]  = ussl[BXt];

          ussl[BXb]  = sqrL*usr[BXb] + sqrR*usl[BXb] + sqrL*sqrR*(wsR - wsL)*sBx;
          ussl[BXb] /= sqrL + sqrR;
          ussr[BXb]  = ussl[BXb];

          scrhL  = SM*Bx1 + vsL*usl [BXt] + wsL*usl [BXb];
          scrhL -= SM*Bx1 + vss*ussl[BXt] + wss*ussl[BXb];

          scrhR  = SM*Bx1 + vsR*usr [BXt] + wsR*usr [BXb];
          scrhR -= SM*Bx1 + vss*ussr[BXt] + wss*ussr[BXb];

          ussl[ENG] = usl[ENG] - sqrL*scrhL*sBx;
          ussr[ENG] = usr[ENG] + sqrR*scrhR*sBx;

          #ifdef GLM_MHD
          ussl[PSI_GLM] = ussr[PSI_GLM] = vl[PSI_GLM];
          #endif

          // 2f7. Flux in the **L or **R regions

          if (SM >= 0.0){           //  Region L**  
            #pragma acc loop seq
            for (nv = 0; nv < NFLX; nv++) {
              flux(i,nv) = fl[nv] + S1L*(ussl[nv]  - usl[nv])
                                    + Sl*(usl[nv] - ul[nv]);
            }
            press(i) = ptl;
          }else{                   // Region R**  
            #pragma acc loop seq
            for (nv = 0; nv < NFLX; nv++) {
              flux(i,nv) = fr[nv] + S1R*(ussr[nv]  - usr[nv])
                                    + Sr*(usr[nv] - ur[nv]);

            }
            press(i) = ptr;
          }
        }  // end if (S1L < 0 S1R > 0) 

      } else {  // USE HLL 

        double scrh = 1.0/(Sr - Sl);
        #pragma acc loop seq
        for (nv = 0; nv < NFLX; nv++) {
          flux(i,nv)  = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fl[nv] - Sl*fr[nv];
          flux(i,nv) *= scrh;
        }
        press(i) = (Sr*ptl - Sl*ptr)*scrh;

      } // End if (!switch_to_hll) 
    }  // end if (SL < 0, SR > 0) 

    // 2g. Compute Powell's / GLM interface states

    #if (DIVB_CONTROL == EIGHT_WAVES)
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
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
    #elif CT_EMF_AVERAGE == UCT_HLLD
    // if (switch_to_hll) CT_Store_UCTHll<dir>  (data, vl, vr, Sl, Sr, ti, tj, tk);
    // else               CT_Store_UCTHlld<dir> (data, ul, ur, vl, vr, Sl, Sr, ti, tj, tk);
    CT_Store_UCTHlld<dir> (data, ul, ur, vl, vr, Sl, Sr, ti, tj, tk);
    
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

#elif EOS == ISOTHERMAL

////////////////////////////////////////////////////////////////////////
/*!
  Solve Riemann problem for the adiabatic MHD equations using the
  four-state HLLD Riemann solver of Miyoshi & Kusano (2005).
  
  \param[in,out] sweep    pointer to Sweep structure
  \param[in]     beg     initial grid index
  \param[out]    end     final grid index
  \param[out]    cmax    1D array of maximum characteristic speeds
  \param[in]     grid    pointer to array of Grid structures.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLLD_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int  nv, i, k , j;
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  int ni, nj;
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
  QuitPLUTO("! HLLD_Solver(): not compatible with UCT_GFORCE", 1);
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
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;

    double usl[NVAR];
    double usr[NVAR];
    double usc[NVAR];

    double vl[NVAR], ul[NVAR], fl[NVAR], Sl, ptl;
    double vr[NVAR], ur[NVAR], fr[NVAR], Sr, ptr;
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

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. GLM pre-conditioner

    #ifdef GLM_MHD
    GLM_Solve (vl, vr, run_config,dir);
    #endif
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    // 2b. - compute sound speed, cons variables and fluxes
    //     - compute HLL speeds.

    STATE_SET(stateL, vl, ul, bgf, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, bgf, x1[ti], x2[tj], x3[tk]);

    #if BACKGROUND_FIELD == YES
    BackgroundField (stateL, run_config); // stateL->bgf = stateR->bgf
    #endif

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    Flux<dir> (stateL, fl, ptl, run_config);
    Flux<dir> (stateR, fr, ptr, run_config);

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    // 2c. Set Alfven and contact waves velocities

    #if BACKGROUND_FIELD == YES
    double B0n = bgf[BXn-BX1];
    double B0t = bgf[BXt-BX1];
    double B0b = bgf[BXb-BX1];
    #endif

    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);
    double scrh  = MAX(fabs(Sl), fabs(Sr));
    cmax(i) = scrh;
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    scrh = 1.0/(Sr - Sl);
    double dul  = Sl - vl[VXn];
    double dur  = Sr - vr[VXn];

    double Bx = (Sr*vr[BXn] - Sl*vl[BXn])*scrh;
    #if BACKGROUND_FIELD == YES
    Bx += B0n;   // total field 
    #endif

    double rho       = (ur[RHO]*dur - ul[RHO]*dul)*scrh;
    flux(i,RHO) = (Sl*ur[RHO]*dur - Sr*ul[RHO]*dul)*scrh;

    double sqrho = sqrt(rho);

    double SM  = flux(i,RHO)/rho;
    double S1L = SM - fabs(Bx)/sqrho;
    double S1R = SM + fabs(Bx)/sqrho;

    // 2d. Compute Flux if SL > 0 or SR < 0

    if (Sl >= 0.0){                     // -- Region L --
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fl[nv];
      press(i) = ptl;

    }else if (Sr <= 0.0) {              // -- Region R  --
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fr[nv];
      press(i) = ptr;

    } else {

      // 2e. Compute Flux when  SL < 0 < SR

      // 2e1. When S1L -> SL or S1R -> SR a degeneracy occurs.
      //      Although Miyoshi & Kusano say that no jump exists, we don't
      //      think this is actually true.
      //      Indeed, vy*, vz*, By*, Bz* cannot be solved independently.
      //      In this case we revert to the HLLC solver of Li (2005),  except
      //      for the term v.B in the region, which we compute in our own way.
      //      Note, that by comparing the expressions of Li (2005) and
      //      Miyoshi & Kusano (2005), the only change involves a
      //      re-definition of By* and Bz* in terms of By(HLL), Bz(HLL).   

      int switch_to_hll = 0;
      if ( (S1L - Sl) <  1.e-4*(SM - Sl) ) switch_to_hll = 1;
      if ( (S1R - Sr) > -1.e-4*(Sr - SM) ) switch_to_hll = 1;

      if (!switch_to_hll){ // Use HLLD 

        flux(i,MXn) = (  Sr*fl[MXn] - Sl*fr[MXn] 
                       + Sr*Sl*(ur[MXn] - ul[MXn]))*scrh;

        press(i) = (Sr*ptl - Sl*ptr)*scrh;
        #ifdef GLM_MHD
        flux(i,BXn)     = fl[BXn];
        flux(i,PSI_GLM) = fl[PSI_GLM];
        #else
        flux(i,BXn) = Sr*Sl*(ur[BXn] - ul[BXn])*scrh;
        #endif

        // 2e2. Compute states in the *L and *R regions.
        //      With bck field splitting, the intermediate
        //      magnetic fields always contain the deviation B1. 

        {
        double scrhL = 1.0/((Sl - S1L)*(Sl - S1R));
        double scrhR = 1.0/((Sr - S1L)*(Sr - S1R));

        #if BACKGROUND_FIELD == YES
        usl[MXt] = rho*vL[VXt] - Bx*(ul[BXt]+B0t)*(SM - vL[VXn])*scrhL;
        usr[MXt] = rho*vR[VXt] - Bx*(ur[BXt]+B0t)*(SM - vR[VXn])*scrhR;

        usl[MXb] = rho*vL[VXb] - Bx*(ul[BXb]+B0b)*(SM - vL[VXn])*scrhL;
        usr[MXb] = rho*vR[VXb] - Bx*(ur[BXb]+B0b)*(SM - vR[VXn])*scrhR;
        #else
        usl[MXt] = rho*vl[VXt] - Bx*ul[BXt]*(SM - vl[VXn])*scrhL;
        usr[MXt] = rho*vr[VXt] - Bx*ur[BXt]*(SM - vr[VXn])*scrhR;

        usl[MXb] = rho*vl[VXb] - Bx*ul[BXb]*(SM - vl[VXn])*scrhL;
        usr[MXb] = rho*vr[VXb] - Bx*ur[BXb]*(SM - vr[VXn])*scrhR;
        #endif

        #if BACKGROUND_FIELD == YES
        usl[BXt] = (ul[BXt]+B0t)/rho*(ul[RHO]*dul*dul - Bx*Bx)*scrhL-B0t;
        usr[BXt] = (ur[BXt]+B0t)/rho*(ur[RHO]*dur*dur - Bx*Bx)*scrhR-B0t;

        usl[BXb] = (ul[BXb]+B0b)/rho*(ul[RHO]*dul*dul - Bx*Bx)*scrhL-B0b;
        usr[BXb] = (ur[BXb]+B0b)/rho*(ur[RHO]*dur*dur - Bx*Bx)*scrhR-B0b;
        #else
        usl[BXt] = ul[BXt]/rho*(ul[RHO]*dul*dul - Bx*Bx)*scrhL;
        usr[BXt] = ur[BXt]/rho*(ur[RHO]*dur*dur - Bx*Bx)*scrhR;

        usl[BXb] = ul[BXb]/rho*(ul[RHO]*dul*dul - Bx*Bx)*scrhL;
        usr[BXb] = ur[BXb]/rho*(ur[RHO]*dur*dur - Bx*Bx)*scrhR;
        #endif
        }

        // 2e3. Compute flux when S1L > 0 or S1R < 0

        if (S1L >= 0.0){       //  ----  Region L*  ---- 

          flux(i,MXt) = fl[MXt] + Sl*(usl[MXt] - ul[MXt]);
          flux(i,MXb) = fl[MXb] + Sl*(usl[MXb] - ul[MXb]);

          flux(i,BXt) = fl[BXt] + Sl*(usl[BXt] - ul[BXt]);
          flux(i,BXb) = fl[BXb] + Sl*(usl[BXb] - ul[BXb]);

        }else if (S1R <= 0.0) {    //  ----  Region R*  ---- 

          flux(i,MXt) = fr[MXt] + Sr*(usr[MXt] - ur[MXt]);
          flux(i,MXb) = fr[MXb] + Sr*(usr[MXb] - ur[MXb]);

          flux(i,BXt) = fr[BXt] + Sr*(usr[BXt] - ur[BXt]);
          flux(i,BXb) = fr[BXb] + Sr*(usr[BXb] - ur[BXb]);

        } else {

          // 2e4. Compute U**

          {
          double sBx = (Bx > 0.0 ? 1.0 : -1.0);

          usc[MXt] = 0.5*(usr[MXt] + usl[MXt] + (usr[BXt] - usl[BXt])*sBx*sqrho);
          usc[MXb] = 0.5*(usr[MXb] + usl[MXb] + (usr[BXb] - usl[BXb])*sBx*sqrho);

          usc[BXt] = 0.5*(usr[BXt] + usl[BXt] + (usr[MXt] - usl[MXt])*sBx/sqrho);
          usc[BXb] = 0.5*(usr[BXb] + usl[BXb] + (usr[MXb] - usl[MXb])*sBx/sqrho);
          }

          flux(i,MXt) = usc[MXt]*SM - Bx*usc[BXt];
          flux(i,MXb) = usc[MXb]*SM - Bx*usc[BXb];
          #if BACKGROUND_FIELD == YES
          flux(i,MXt) -= Bx1*B0t;
          flux(i,MXb) -= Bx1*B0b;
          #endif

          flux(i,BXt) = usc[BXt]*SM - Bx*usc[MXt]/rho;
          flux(i,BXb) = usc[BXb]*SM - Bx*usc[MXb]/rho;
          #if BACKGROUND_FIELD == YES
          flux(i,BXt) += B0t*SM;
          flux(i,BXb) += B0b*SM;
          #endif

        } // end if (S1L < 0 S1R > 0) 

        // 2e5. Revert do HLL if  |Sa-S| < tol

      } else { // Use HLL 

        #pragma acc loop seq
        NFLX_LOOP(nv) {
          flux(i,nv) = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fl[nv] - Sl*fr[nv];
          flux(i,nv) *= scrh;
        }
        press(i) = (Sr*ptl - Sl*ptr)*scrh;

      } // end if (!switch_to_hll) 

    }  // end if (SL < 0 SR > 0) 

    // 2f. Compute tracer fluxes

    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }

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
    #elif CT_EMF_AVERAGE == UCT_HLLD
    CT_Store_UCTHlld<dir> (data, ul, ur, vl, vr, 
                            Sl, Sr, ti, tj, tk);
    #endif
    #endif

    // 2h. Compute Powell's interface state

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

  }}}  // End for i, j, k

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
#endif
