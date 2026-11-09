////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  High-order reconstructions (MP5 / WENO)

  Compute interface states using:
  - MP5 5th-order monotonicity-preserving
    schemes of Suresh & Huynh (1997); 
  - WENO-Z (Borges et al. 2008).
  - WENO-3 ().
  - LimO3.

  The reconstruction workflow can be divided into 4 subtasks:

  1. Transform to 4vel
  2. Apply reconstruction,  apply order reduction
  3. Adjust L/R states in case of unphysical values, normal field assignment
  4. Transform to 3vel

  \author A. Mignone  (andrea.mignone@unito.it)
          V. Berta    (vittoria.berta@unito.it)
          M. Rossazza (marco.rossazza@unito.it)

  \date    Dec 17, 2025

  \b References
     - "Accurate Monotonicity-Preserving Schemes 
        with Runge-Kutta Time Stepping"
        Suresh, JCP (1997) 136, 83
     - "An improved weighted essentially non-oscillatory scheme 
        for hyperbolic conservation laws"
        Borges et al, JCP (2008) 227, 3191
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "States/states.hpp"

// Explicit template instantiations
template void Reconstruct<IDIR> (Data *, Grid *, RBox *);
template void Reconstruct<JDIR> (Data *, Grid *, RBox *);
template void Reconstruct<KDIR> (Data *, Grid *, RBox *);

#if CHAR_LIMITING == NO
////////////////////////////////////////////////////////////////////////
/*!
  Compute states using fifth-order reconstruction in primitive
  variables.
  
  \return This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void Reconstruct (Data *data, Grid *grid, RBox *box)
{
  int i, j, k, nv;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;

  int  ni, nj;

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  CT_ASSIGN(constexpr int BXn = BX1 + dir;,
            constexpr int EXn = EX1 + dir;)
  #endif

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  ibeg--; // Reconstruct one point to the left
  iend++; // Reconstruct one point to the right

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("MP5_Reconstruct",GREEN);

  // 1a. Reflect / symmetrized profiles when zones are 
  //     flagged with INTERNAL_BOUNDARY.

  #if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
  #endif

  // 1b. Transform to four-velocity if required.

  #if RECONSTRUCT_4VEL
  #if (RECONSTRUCTION == LimO3) || (RECONSTRUCTION == WENO3)
  int stncl = 1;    // reconstruction stencil
  #elif (RECONSTRUCTION == WENOZ) || (RECONSTRUCTION == MP5)
  int stncl = 2;    // reconstruction stencil
  #endif

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data) 
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);

    #pragma acc loop vector
    for (i = ibeg-stncl; i <= iend+stncl; i++) ConvertTo4vel(v, i);
  }}
  #endif

  // ------------------------------------------------------
  // 2. Main Reconstruction Loop (primitive variables)
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data, grid) 
  for (nv = 0; nv < NVAR; nv++) {
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  
    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
   
    double dvp, dvm, dv_lim;
    double cp, cm;
    double w[5];

    // 2a. Copy variables from the stencil region
    //     for reconstruction

    #if (RECONSTRUCTION == MP5) || (RECONSTRUCTION == WENOZ)
    w[4] = v(i+2,nv);
    w[3] = v(i+1,nv);
    w[2] = v(i,  nv);
    w[1] = v(i-1,nv);
    w[0] = v(i-2,nv);
    #elif (RECONSTRUCTION == LimO3) || (RECONSTRUCTION == WENO3)
    w[3] = v(i+1,nv);
    w[2] = v(i,  nv);
    w[1] = v(i-1,nv);
    #endif

    cp = cm = 2.0;

    // 2b. Reduce order in flagged zones

    if (flag(i) & FLAG_FLAT) {
      vL(i+1,nv) = vR(i,nv) = v(i,nv);
    } else if (flag(i) & FLAG_MINMOD){
      dvp = v(i+1,nv) - v(i,nv);
      dvm = v(i,nv)   - v(i-1,nv);
      dv_lim    = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,nv) = v(i,nv) + 0.5*dv_lim;
      vR(i,nv)   = v(i,nv) - 0.5*dv_lim;
    }else
    #if (defined HIGH_ORDER) && (HO_ORDER_REDUCTION == WENO3)
      if (flag(i) & FLAG_HO_LAP_LIMITER){
        vL(i+1,nv) = WENO3_States (&w[2], +1, grid->dx[dir][i]);
        vR(i,nv)   = WENO3_States (&w[2], -1, grid->dx[dir][i]);
      }else
    #endif

    // 2c. Apply high-order reconstruction

    #if RECONSTRUCTION == LimO3
    {
    vL(i+1,nv) = LimO3_States (&w[2], +1, grid->dx[dir][i]);
    vR(i,nv)   = LimO3_States (&w[2], -1, grid->dx[dir][i]);
    }
    #elif RECONSTRUCTION == WENO3
    {
    vL(i+1,nv) = WENO3_States (&w[2], +1, grid->dx[dir][i]);
    vR(i,nv)   = WENO3_States (&w[2], -1, grid->dx[dir][i]);
    }
    #elif RECONSTRUCTION == WENOZ
    {
    vL(i+1,nv) = WENOZ_States (&w[2], +1);
    vR(i,nv)   = WENOZ_States (&w[2], -1);
    }
    #elif RECONSTRUCTION == MP5
    {
    vL(i+1,nv) = MP5_States (&w[2], +1);
    vR(i,nv)   = MP5_States (&w[2], -1);
    }
    #else
      #error Invalid RECONSTRUCTION in mp5_states
    #endif
  }}}} // End loop (nv, k, j, i) main reconstruction loop 

  // ------------------------------------------------------
  // 3. Adjustment Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  
    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*offset;

    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D En(&data->sweep->En[offset],ni);
    #endif
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    double dvp, dvm, dv_lim;

    // 3a. Correct for negative press. or dens.

    if (vL(i+1,RHO) < 0.0 || vR(i,RHO) < 0.0){
      dvp = v(i+1,RHO) - v(i,RHO);
      dvm = v(i,RHO)   - v(i-1,RHO);
      dv_lim      = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,RHO) = v(i,RHO) + 0.5*dv_lim;
      vR(i,RHO)   = v(i,RHO) - 0.5*dv_lim;
    }
    #if HAVE_ENERGY
    if (vL(i+1,PRS) < 0.0 || vR(i,PRS) < 0.0){
      dvp = v(i+1,PRS) - v(i,PRS);
      dvm = v(i,PRS)   - v(i-1,PRS);
      dv_lim    = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,PRS) = v(i,PRS) + 0.5*dv_lim;
      vR(i,PRS)   = v(i,PRS) - 0.5*dv_lim;
    }
    #endif
    #if RADIATION
    if (vL(i+1,ENR) < 0.0 || vR(i,ENR) < 0.0){
      dvp = v(i+1,ENR) - v(i,ENR);
      dvm = v(i,ENR)   - v(i-1,ENR);
      dv_lim    = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,ENR) = v(i,ENR) + 0.5*dv_lim;
      vR(i,ENR)   = v(i,ENR) - 0.5*dv_lim;
    }
    #endif

    // 3b. Apply relativistic velocity limiter, 
    //     convert L/R states to three-velocity 

    #if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
    VelocityLimiter (v, vL, vR, i);
    #endif
    #if RADIATION
    RadFluxLimFlatten (v, vL, vR, i);
    #endif

    #if RECONSTRUCT_4VEL
    ConvertTo3vel(vL, i+1);
    ConvertTo3vel(vR, i);
    #endif

    // 3c. Assign face-centered magnetic field

    // #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    // CT_ASSIGN(vL(i,BXn) = vR(i,BXn) = Bn(i);    ,
    //           vL(i,EXn) = vR(i,EXn) = En(i);)
    // #endif

    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    CT_ASSIGN(vL(i+1,BXn) = Bn(i+1);  
              vR(i,BXn)   = Bn(i);    ,
              vL(i+1,EXn) = En(i+1);
              vR(i,EXn)   = En(i);)
    #endif
  }}} // End loop (i,j,k)

  // 4. Convert cell-center back to three-velocity 

  #if RECONSTRUCT_4VEL
  #pragma acc parallel loop collapse(2) present(data) 
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    
    #pragma acc loop vector
    for (i = ibeg-stncl; i <= iend+stncl; i++) ConvertTo3vel(v, i);
  }} 
  #endif

  mynvtxstop_();

}
#endif

#if CHAR_LIMITING == YES
////////////////////////////////////////////////////////////////////////
/*!
  Compute states using fifth-order reconstruction in characteristic 
  variables.
  
  \return This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void Reconstruct (Data *data, Grid *grid, RBox *box)
{
  int i, j, k, nv, l; 
  int  ni, nj;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  CT_ASSIGN(constexpr int BXn = BX1 + dir;,
            constexpr int EXn = EX1 + dir;)
  #endif

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  RunConfig &run_config = *(data->run_config);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  ibeg--; // Reconstruct one point to the left
  iend++; // Reconstruct one point to the right

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("MP5_Reconstruct",GREEN);

  // 1. Reflect / symmetrized profiles when zones are 
  //     flagged with INTERNAL_BOUNDARY.

  #if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
  #endif

  // ------------------------------------------------------
  // 2. Main Reconstruction Loop (characteristic variables)
  // ------------------------------------------------------ 

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    
    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*offset;

    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D En(&data->sweep->En[offset],ni);
    #endif
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    double L[NFLX][NFLX], R[NFLX][NFLX];
    double vm[NVAR], vc[NVAR], vp[NVAR];

    int    ti, tj, tk;
    double dvp, dvm, dv, dv_lim;
    double cp, cm;
    double w[5][NVAR], w1[5];
    double wp[NVAR], wm[NVAR];
    double lambda[NFLX];
    double dp, dm;

    State state;

    cp = cm = 2.0;

    if (dir == IDIR)  {ti = i; tj = j; tk = k;}
    if (dir == JDIR)  {ti = j; tj = i; tk = k;}
    if (dir == KDIR)  {ti = j; tj = k; tk = i;}

    // 2a. Reduce order if zone has been flagged

    if (flag(i) & FLAG_FLAT) {
      #pragma acc loop
      NVAR_LOOP(nv) vL(i+1,nv) = vR(i,nv) = v(i,nv);
    }else if (flag(i) & FLAG_MINMOD) {
      #pragma acc loop
      NVAR_LOOP(nv){
        dvp = v(i,nv)   - v(i-1,nv);
        dvm = v(i+1,nv) - v(i,nv);
        dv  = MINMOD_LIMITER(dvp, dvm);
        vL(i+1,nv) = v(i,nv) + 0.5*dv;
        vR(i,nv)   = v(i,nv) - 0.5*dv;
      }
    } else {

      // 2a. Copy variables from the stencil region
      //    for reconstruction

      #if (RECONSTRUCTION == MP5) || (RECONSTRUCTION == WENOZ)
      double vpp[NVAR], vmm[NVAR];

      #pragma acc loop seq
      NVAR_LOOP(nv){  
        vpp[nv] = v(i+2,nv);
        vp[nv]  = v(i+1,nv);
        vc[nv]  = v(i,nv);
        vm[nv]  = v(i-1,nv);
        vmm[nv] = v(i-2,nv);
      }

      // 2b. Obtain characteristic variables by
      //     projecting on left eigenvectors.
   
      STATE_SET (state, vc, NULL, NULL, x1[ti], x2[tj], x3[tk]);
      PrimEigenvectors<dir> (state, run_config, L, R, lambda);

      PrimToChar<dir> (L, vpp, w[4]);
      PrimToChar<dir> (L, vp,  w[3]);
      PrimToChar<dir> (L, vc,  w[2]);  // '2' = stencil center
      PrimToChar<dir> (L, vm,  w[1]);
      PrimToChar<dir> (L, vmm, w[0]);
      #elif (RECONSTRUCTION == LimO3) || (RECONSTRUCTION == WENO3)
      #pragma acc loop seq
      NVAR_LOOP(nv){  
        vp[nv]  = v(i+1,nv);
        vc[nv]  = v(i,nv);
        vm[nv]  = v(i-1,nv);
      }

      STATE_SET (state, vc, NULL, NULL, x1[ti], x2[tj], x3[tk]);
      PrimEigenvectors<dir> (state, run_config, L, R, lambda);

      PrimToChar<dir> (L, vp,  w[3]);
      PrimToChar<dir> (L, vc,  w[2]);  // '2' = stencil center
      PrimToChar<dir> (L, vm,  w[1]);
      #endif

      // 2c. Reconstruct characteristic variables

      #pragma acc loop
      NVAR_LOOP(nv){

        #if (RECONSTRUCTION == MP5) || (RECONSTRUCTION == WENOZ)
        w1[4] = w[4][nv];
        w1[3] = w[3][nv];
        w1[2] = w[2][nv];
        w1[1] = w[1][nv];
        w1[0] = w[0][nv];
        #elif (RECONSTRUCTION == LimO3) || (RECONSTRUCTION == WENO3)
        w1[3] = w[3][nv];
        w1[2] = w[2][nv];
        w1[1] = w[1][nv];
        #endif

        #if (defined HIGH_ORDER) && (HO_ORDER_REDUCTION == WENO3)
        if (flag(i) & FLAG_HO_LAP_LIMITER){
          wp[nv] = WENO3_States (&w1[2], +1, grid->dx[dir][i]);
          wm[nv] = WENO3_States (&w1[2], -1, grid->dx[dir][i]);
        }else
        #endif

        #if RECONSTRUCTION == LimO3
        {
        wp[nv] = LimO3_States (&w1[2], +1, grid->dx[dir][i]);
        wm[nv] = LimO3_States (&w1[2], -1, grid->dx[dir][i]);
        }
        #elif RECONSTRUCTION == WENO3
        {
        wp[nv] = WENO3_States (&w1[2], +1, grid->dx[dir][i]);
        wm[nv] = WENO3_States (&w1[2], -1, grid->dx[dir][i]);
        }
        #elif RECONSTRUCTION == WENOZ
        {
        wp[nv] = WENOZ_States (&w1[2], +1);
        wm[nv] = WENOZ_States (&w1[2], -1);
        }
        #elif RECONSTRUCTION == MP5
        {
        wp[nv] = MP5_States (&w1[2], +1);
        wm[nv] = MP5_States (&w1[2], -1);
        }
        #else
          #error Invalid RECONSTRUCTION in mp5_states
        #endif

      // -- Lower order if flagged -- 

        // if (flag(i) & FLAG_MINMOD){
        //   dvp    = w1[3] - w1[2];
        //   dvm    = w1[2] - w1[1];
        //   dv_lim = MINMOD_LIMITER(dvp, dvm);
        //   wp[nv] = w1[2] + 0.5*dv_lim;
        //   wm[nv] = w1[2] - 0.5*dv_lim;
        // }

      }  // End NVAR_LOOP(nv) on characteristics 

      // 2d. Obtain L/R states in primitive variables

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        dp = dm = 0.0;
        #pragma acc loop seq
        for (l = 0; l < NFLX; l++){
          dp += wp[l]*R[nv][l];
          dm += wm[l]*R[nv][l];
        }
        vL(i+1,nv) = dp;
        vR(i,nv)   = dm;
      }

      #if NSCL > 0
      #pragma acc loop
      NSCL_LOOP(nv){
        vL(i+1,nv) = wp[nv];
        vR(i,nv)   = wm[nv];
      }
      #endif

      #if RADIATION
      #pragma acc loop
      NRAD_LOOP(nv){
        vL(i+1,nv) = wp[nv];
        vR(i,nv)   = wm[nv];
      }
      #endif
    } // End if full reconstruction 

    // 2e. Correct for negative pressures or densities
    
    if (vL(i+1,RHO) < 0.0 || vR(i,RHO) < 0.0){
      dvp = v(i+1,RHO) - v(i,RHO);
      dvm = v(i,RHO)   - v(i-1,RHO);
      dv_lim      = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,RHO) = v(i,RHO) + 0.5*dv_lim;
      vR(i,RHO)   = v(i,RHO) - 0.5*dv_lim;
    }
    #ifdef ENG
    if (vL(i+1,PRS) < 0.0 || vR(i,PRS) < 0.0){
      dvp = v(i+1,PRS) - v(i,PRS);
      dvm = v(i,PRS)   - v(i-1,PRS);
      dv_lim    = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,PRS) = v(i,PRS) + 0.5*dv_lim;
      vR(i,PRS)   = v(i,PRS) - 0.5*dv_lim;
    }
    #endif
    #if RADIATION
    if (vL(i+1,ENR) < 0.0 || vR(i,ENR) < 0.0){
      dvp = v(i+1,ENR) - v(i,ENR);
      dvm = v(i,ENR)   - v(i-1,ENR);
      dv_lim    = MINMOD_LIMITER (dvp, dvm);
      vL(i+1,ENR) = v(i,ENR) + 0.5*dv_lim;
      vR(i,ENR)   = v(i,ENR) - 0.5*dv_lim;
    }
    #endif

    // 2f. Apply relativistic Limiter to avoid
    //     super-luminal velocities.

    #if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
    VelocityLimiter (v, vL, vR, i);
    #endif
    #if RADIATION
    RadFluxLimFlatten (v, vL, vR, i);
    #endif

    // 2g. Assign face-centered magnetic field

    // #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    // CT_ASSIGN(vL(i,BXn) = vR(i,BXn) = Bn(i);    ,
    //           vL(i,EXn) = vR(i,EXn) = En(i);)
    // #endif

    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    CT_ASSIGN(vL(i+1,BXn) = Bn(i+1);  
              vR(i,BXn)   = Bn(i);    ,
              vL(i+1,EXn) = En(i+1);
              vR(i,EXn)   = En(i);)
    #endif

  }}} // End loop on (i,j,k)

  mynvtxstop_();

}
#endif
