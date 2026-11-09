////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Piecewise linear reconstruction.

  Compute interface states using piecewise linear reconstruction
  inside each zone.
  Reconstruction is performed in primitive variable when
   <tt> CHAR_LIMITING == NO </tt>) or characteristic variables
  when <tt> CHAR_LIMITING == YES </tt>.

  The convention used throughout is that \b vL and \b vR are left and
  right states with respect to the \e interface, while \b vm and \b vp
  refer to the cell \e center:

                      vL(i)-> <-vR(i)
       |----------*----------|----------*----------|
                (i-1)      i-1/2        i

  The default setting applies VanLeer limiter to all variables and 
  MC limiter to tracers.

  A stencil of 3 zones is required for all limiters.

  \b References
     - "High-order conservative reconstruction schemes for finite
        volume methods in cylindrical and spherical coordinates"
        A. Mignone, JCP (2014), 270, 784.

  \authors A. Mignone  (andrea.mignone@unito.it) \n
           M. Rossazza (marco.rossazza@unito.it)
  \date    Dec 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// Explicit template instantiations
template void Reconstruct<IDIR> (Data *, Grid *, RBox *);
template void Reconstruct<JDIR> (Data *, Grid *, RBox *);
template void Reconstruct<KDIR> (Data *, Grid *, RBox *);

#if CHAR_LIMITING == NO
////////////////////////////////////////////////////////////////////////
/*!
  Perform 1D linear reconstruction of primitive variables
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
  constexpr int BXn = BX1 + dir;
  #endif
  #if PHYSICS == ResRMHD && DIVE_CONTROL == CONSTRAINED_TRANSPORT
  constexpr int EXn = EX1 + dir;
  #endif

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  ibeg--; // Reconstruct one point to the left
  iend++; // Reconstruct one point to the right

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("PLM_Reconstruct",GREEN);

  // 1a. Reflect / symmetrized profiles when zones are 
  //     flagged with INTERNAL_BOUNDARY.

#if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
#endif

  // 1b. Transform to four-velocity if required.

  #if RECONSTRUCT_4VEL
  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data) 
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);

    #pragma acc loop vector
    for (i = 0; i < ni; i++) ConvertTo4vel(v, i);
  }}
  #endif

  // ------------------------------------------------------
  // 2. Main Reconstruction Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (nv = 0; nv < NVAR; nv++){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*offset;

    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){
      double dvp, dvm, dv;
      double cp = 2.0, cm = 2.0;

      // 2a. Compute forward (dvp) and backward (dvm) differences

      dvp = v(i+1,nv) - v(i,nv);
      dvm = v(i,nv)   - v(i-1,nv);

      // 2b. Apply limiter / reduce order if needed

      if (flag(i) & FLAG_FLAT) {
        dv = 0.0;
      }else if (flag(i) & FLAG_MINMOD) {
        dv = MINMOD_LIMITER (dvp, dvm);
      }else{
        SET_LIMITER(dv, dvp, dvm, cp, cm);
      }

      // 2c. Construct left and right states

      vL(i+1,nv) = v(i,nv) + 0.5*dv;
      vR(i,nv)   = v(i,nv) - 0.5*dv;

    }  // End loop on i
  }}}  // End loop on j,k,nv 

  // ------------------------------------------------------
  // 3. Adjustment Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D En(&data->sweep->En[offset],ni);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    // 3a. Relativistic velocity limiter 
    //     or transform back to three-velocity 

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

    // 3b. Assign face-centered magnetic field

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

  }}} // End Loop on (i,j,k)

  // ------------------------------------------------------
  // 4. Transform back to three-velocity 
  // ------------------------------------------------------ */

  #if RECONSTRUCT_4VEL
  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data) 
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    
    #pragma acc loop vector
    for (i = 0; i < ni; i++) ConvertTo3vel(v, i);
  }} 
  #endif

  mynvtxstop_();
}
#endif

#if CHAR_LIMITING == YES
////////////////////////////////////////////////////////////////////////
/*!
  Perform 1D linear reconstruction of primitive variables
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void Reconstruct (Data *data, Grid *grid, RBox *box)
{
  int i, j, k, nv, l;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int  ni, nj;

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  constexpr int BXn = BX1 + dir;
  #endif
  #if PHYSICS == ResRMHD && DIVE_CONTROL == CONSTRAINED_TRANSPORT
  constexpr int EXn = EX1 + dir;
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

  mynvtxstart_("PLM_Reconstruct",GREEN);

  // 1. Reflect / symmetrized profiles when zones are 
  //    flagged with INTERNAL_BOUNDARY.

  #if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
  #endif

  // ------------------------------------------------------
  // 2. Main Reconstruction Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) private(nv, l)
  #pragma acc parallel loop collapse(3) present(data, data->run_config, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      int    ti, tj, tk;

      long int offset = ni*(j + nj*k);
      long int offset1 = NVAR*ni*(j + nj*k);

      Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
      Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
      #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
      Ary1D Bn(&data->sweep->Bn[offset],ni);
      Ary1D En(&data->sweep->En[offset],ni);
      #endif
      Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

      double L[NFLX][NFLX], R[NFLX][NFLX];
      double dv_lim[NVAR], d2v, vc[NVAR];
      double dw_lim[NVAR], dwp[NVAR], dwm[NVAR];
      double dvp[NVAR], dvm[NVAR];
      double lambda[NFLX];
      double dp, dm, dc;
      double cp, cm;

      State state;

      if (dir == IDIR)  {ti = i; tj = j; tk = k;}
      if (dir == JDIR)  {ti = j; tj = i; tk = k;}
      if (dir == KDIR)  {ti = j; tj = k; tk = i;}

      // 2a. Project forward, backward (and centered) 
      //     undivided differences of prim. variables 
      //     along characteristics, dw(k) = L(k).dv

      cp = cm = 2.0;
      dp = dm = 0.5;

      #pragma acc loop seq
      NVAR_LOOP(nv) vc[nv]  = v(i,nv);

      STATE_SET (state, vc, NULL, NULL, x1[ti], x2[tj], x3[tk]);
      PrimEigenvectors<dir> (state, run_config, L, R, lambda);

      #pragma acc loop seq
      NVAR_LOOP(nv) {
        dvp[nv] = v(i+1,nv) - v(i,nv);
        dvm[nv] = v(i,nv)   - v(i-1,nv);
      }

      PrimToChar<dir> (L, dvm, dwm);
      PrimToChar<dir> (L, dvp, dwp);

      // 2b. Apply slope limiter to characteristic differences 
      //     for nv < NFLX, dw = Limiter(dwp, dwm).

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        SET_LIMITER(dw_lim[nv], dwp[nv], dwm[nv], cp, cm);
      }

      // 2c. Reduce order if necessary

      if (flag(i) & FLAG_FLAT) {
        #pragma acc loop
        NVAR_LOOP(nv) dw_lim[nv] = 0.0;
      }else if (flag(i) & FLAG_MINMOD) {
        #pragma acc loop
        NVAR_LOOP(nv) {
          SET_MM_LIMITER(dw_lim[nv], dwp[nv], dwm[nv], cp, cm);
        }
      }

      // 2d. Project limited slopes in characteristic vars on 
      //     right eigenvectors to obtain primitive slopes: 
      //
      //     dv = \sum dw.R
      //
      //     Also, enforce monotonicity in primitive variables 
      //     as well.

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        dc = 0.0;
        #pragma acc loop seq
        for (l = 0; l < NFLX; l++) dc += dw_lim[l]*R[nv][l];
        
        if (dvp[nv]*dvm[nv] > 0.0){
          d2v        = ABS_MIN(cp*dvp[nv], cm*dvm[nv]);
          dv_lim[nv] = MINMOD_LIMITER(d2v, dc);
        }else dv_lim[nv] = 0.0;
      }

      // 2e. Repeat construction for passive scalars (tracers),
      //     since the primitive variable is the same as the 
      //     characteristic one. We use the MC limiter.

      #if NSCL > 0
      #pragma acc loop seq
      NSCL_LOOP(nv){
        SET_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
      }
      #endif

      #if RADIATION
      #pragma acc loop seq
      NRAD_LOOP(nv){
        SET_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
      }
      #endif

      // 2f. Build L/R states

      #pragma acc loop seq
      NVAR_LOOP(nv){
        vL(i+1,nv) = v(i,nv) + dv_lim[nv]*dp;
        vR(i,nv)   = v(i,nv) - dv_lim[nv]*dm;
      }

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

    } // End loop on normal direction 
  }} // End loop on (j,k)

mynvtxstop_();

}
#endif
