////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Piecewise parabolic reconstruction.

  Compute interface states using piecewise parabolic reconstruction 
  inside each zone.
  Reconstruction is performed in primitive variable (when 
   <tt>CHAR_LIMITING == NO</tt>) or characteristic variables 
  (when <tt>CHAR_LIMITING == YES</tt>).
  The reconstruction process follows the following steps:
  - compute unlimited 4-th or 5-th order interface values from within 
    the zone (\f${\rm vp} = v_{i,+},\; {\rm vm} = v_{i,-} \f$).
    
  - interface states must lie between adjacent cell averages. 
    This is achieved by working on the increments,
    \f[
      \left\{\begin{array}{lll}
           \min(v_{i+1},v_i) \le v_{i,+} \le \max(v_{i+1},v_i)
           & \qquad \Longrightarrow\qquad &
           \delta v_{i,+} \to {\rm minmod}(\delta v_{i,+}, \Delta v_{i+\HALF})
       \\[0.5em]
           \min(v_{i-1},v_i) \le v_{i,-} \le \max(v_{i-1},v_i)
           & \qquad \Longrightarrow\qquad &
           \delta v_{i,-} \to {\rm minmod}(\delta v_{i,-}, \Delta v_{i+\HALF})
       \end{array}\right.
    \f]
    where \f$\delta v_{i,\pm} = v_{i,\pm} - v_i\f$.
    This step is equivalent to Eq. [45] of Mignone (JCP, 2014).
    The previous operation replaces, starting with PLUTO 4.1, the 
    conventional van-Leer limiting used in the original PPM formulation.
    
  - the parabola does not take any extrema inside the zone, see 
    Eq. [46] in Mignone (JCP, 2014).

  When interpolation is carried out in characteristic variables, we 
  first compute \c vp and \c vm (unlimited) in primitive variables and 
  then project the increments in characteristic space:
  \f[
      \delta w_{i,\pm} = \vec{l}_i\cdot\delta v_{i,\pm}
  \f]
  At this point, Eq. [45] (written in terms of increments) is 
  imposed on characteristic variables while Eq. [46] can be imposed
  in characteristic variables (<tt>PRIMITIVE_LIM == 0</tt>), 
  primitive (<tt>PRIMITIVE_LIM == 1</tt>) or both 
  (<tt>PRIMITIVE_LIM == 2</tt>).
  
  \author A. Mignone  (andrea.mignone@unito.it)
          M. Rossazza (marco.rossazza@unito.it)

  \date    Dec 17, 2025

  \b References
     - "High-order conservative reconstruction schemes for finite
        volume methods in cylindrical and spherical coordinates"
        A. Mignone, JCP (2014), 270, 784.
     - "Accurate monotonicity- and extrema-preserving methods
        through adaptive nonlinear hybridizations"
        Rider, Greenough and Kamm, JCP (2007) 225, 1827
     - "A limiter for PPM that preserved accuracy at smooth extrema"
        Colella and Sekora, JCP (2008) 227, 7069
     
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
  Perform 1D parabolic reconstruction of primitive variables
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void Reconstruct (Data *data, Grid *grid, RBox *box)
{
  int i, j, k, nv;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int ni,nj;

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  CT_ASSIGN(constexpr int BXn = BX1 + dir;,
            constexpr int EXn = EX1 + dir;)
  #endif

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  ibeg--; // Reconstruct one point to the left
  iend++; // Reconstruct one point to the right

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("PPM_Reconstruct",GREEN);

  // 1a. Reflect / symmetrized profiles when zones are 
  //     flagged with INTERNAL_BOUNDARY.

  #if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
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
    long int offset1 = NVAR*offset;

    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
   
    double dvm, dvp, dv, vp, vm;
    double cm = 2.0;
    double cp = 2.0;

    // 2a. Reduce order if necessary
       
    if (flag(i) & FLAG_FLAT) {
      vR(i,nv) = vL(i+1,nv) = v(i,nv);
    } else if (flag(i) & FLAG_MINMOD) {
      dvp = v(i+1,nv) - v(i,nv);
      dvm = v(i,nv)   - v(i-1,nv);
      dv  = MINMOD_LIMITER(dvp, dvm);
      vL(i+1,nv) = v(i,nv) + 0.5*dv;
      vR(i,nv)   = v(i,nv) - 0.5*dv;
    } else {

      // 2b. Initialize the left and right states 
      //      to unique interface value.

      vp  =  (-v(i-1,nv) + 7.0*(v(i,nv) + v(i+1,nv)) - v(i+2,nv))/12.0;
      vm  =  (-v(i-2,nv) + 7.0*(v(i-1,nv) + v(i,nv)) - v(i+1,nv))/12.0;

      // 2c. Constraint values to lie between adjacent cell averages.
      //      Same as imposing 
      //      MIN[v(i),v(i+1)] < vp < MAX[v(i), v(i+1)]

      dv  = v(i+1,nv) - v(i,nv);
      dvp = vp - v(i,nv);                       
      vp  = v(i,nv) + MINMOD_LIMITER(dvp, dv); 

      dv  = v(i,nv) - v(i-1,nv);
      dvm = vm - v(i-1,nv);
      vm  = v(i-1,nv) + MINMOD_LIMITER(dvm, dv);

      // 2d. Apply parabolic limiter

      dvp = vp - v(i,nv);
      dvm = vm - v(i,nv);

      if (dvp*dvm >= 0.0) dvp = dvm = 0.0;
      else{
        if      (fabs(dvp) >= cm*fabs(dvm)) dvp = -cm*dvm;
        else if (fabs(dvm) >= cp*fabs(dvp)) dvm = -cp*dvp;
      }

      //2e. Construct left and right states

      vL(i+1,nv) = v(i,nv) + dvp;
      vR(i,nv)   = v(i,nv) + dvm;
    }
  }}}} 

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

    // 3a. Relativistic Limiter

    #if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
    VelocityLimiter (v, vL, vR, i);
    #endif
    #if RADIATION
    RadFluxLimFlatten (v, vL, vR, i);
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
  }}}

  mynvtxstop_();
}
#endif


#if CHAR_LIMITING == YES
#define PARABOLIC_LIM  1
////////////////////////////////////////////////////////////////////////
/*!
  Perform 1D parabolic reconstruction of characteristic variables
*/
////////////////////////////////////////////////////////////////////////
template<int dir> 
void Reconstruct (Data *data, Grid *grid, RBox *box)
{
  int i, j, k, nv;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int ni,nj;
  
  #if NVAR != NFLX
  constexpr int VXn = VX1 + dir;
  #endif
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

    long int offset  = ni*(j + nj*k);
    long int offset1 = NVAR*offset;

    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D En(&data->sweep->En[offset],ni);
    #endif

    double L[NFLX][NFLX], R[NFLX][NFLX];
    double vppm4_p[NVAR], vppm4_m[NVAR];
    double dvF_p[NVAR], dvF_m[NVAR];
    int    ti, tj, tk;
    double dvp[NVAR], dvm[NVAR], dv;
    double cm = 2.0;
    double cp = 2.0;
    double dwp[NVAR], dwm[NVAR];
    double dwp1[NVAR], dwm1[NVAR];
    double lambda[NVAR], vc[NVAR];
    double dp, dm;

    State state;

    if (dir == IDIR)  {ti = i; tj = j; tk = k;}
    if (dir == JDIR)  {ti = j; tj = i; tk = k;}
    if (dir == KDIR)  {ti = j; tj = k; tk = i;}

    // 2a. Reduce order in flagged zones
      
    if (flag(i) & FLAG_FLAT) {
      //#pragma acc loop
      NVAR_LOOP(nv){
          vR(i,nv) = vL(i+1,nv) = v(i,nv);
      }
    } else if (flag(i) & FLAG_MINMOD) {
      //#pragma acc loop
      for (nv = 0; nv < NVAR; nv++){
        dvp[nv] = v(i+1,nv) - v(i,nv);
        dvm[nv] = v(i,nv)   - v(i-1,nv);
        dv  = MINMOD_LIMITER(dvp[nv], dvm[nv]);
        vL(i+1,nv) = v(i,nv) + 0.5*dv;
        vR(i,nv)   = v(i,nv) - 0.5*dv;
      }
    } else{

      // 2b. Compute unlimited interface values
      //     in primitive variables

      NVAR_LOOP(nv){
        vc[nv]    = v(i,nv);
        dvF_p[nv] = v(i+1,nv) - v(i,nv);
        dvF_m[nv] = v(i,nv) - v(i-1,nv);
        vppm4_p[nv] = (-v(i-1,nv) + 7.0*(v(i,nv) + v(i+1,nv)) - v(i+2,nv))/12.0;
        vppm4_m[nv] = (-v(i-2,nv) + 7.0*(v(i-1,nv) + v(i,nv)) - v(i+1,nv))/12.0;
      }

      STATE_SET (state, vc, NULL, NULL, x1[ti], x2[tj], x3[tk]);
      PrimEigenvectors<dir> (state, run_config, L, R, lambda);
            
      #if NVAR != NFLX
      for (nv = NFLX; nv < NVAR; nv++) lambda[nv] = v(i,VXn);
      #endif

      
      // 2c. Project unlimited increments (vp - v) and 
      //     (vm - v) along characteristics and apply limiter.

      NVAR_LOOP(nv){
        dvp[nv] = vppm4_p[nv] - v(i,nv);
        dvm[nv] = vppm4_m[nv] - v(i,nv);
      }

      // Project to characteristic spcae 

      PrimToChar<dir> (L, dvm, dwm);
      PrimToChar<dir> (L, dvp, dwp);

      PrimToChar<dir> (L, dvF_m, dwm1);
      PrimToChar<dir> (L, dvF_p, dwp1);

      // 2d. Enforce monotonicity by constraining interface 
      //     values to lie between adjacent cell averages. 
      //     Same as Eq. [45] in Mignone JCP (2014) 270, 
      //     rewritten by subtracting Q_i from both sides.

      NVAR_LOOP(nv){
        dwp[nv] = MINMOD_LIMITER(dwp[nv],  dwp1[nv]);
        dwm[nv] = MINMOD_LIMITER(dwm[nv], -dwm1[nv]);
      }

      // cm and cp are coefficients used in the parabolic limiter 

      #if PARABOLIC_LIM == 0 || PARABOLIC_LIM == 2
      NVAR_LOOP(nv){
        if (dwp[nv]*dwm[nv] >= 0.0) dwm[nv] = dwp[nv] = 0.0;
        else{
          if      (fabs(dwp[nv]) >= cm*fabs(dwm[nv])) dwp[nv] = -cm*dwm[nv];
          else if (fabs(dwm[nv]) >= cp*fabs(dwp[nv])) dwm[nv] = -cp*dwp[nv];
        }
      }
      #endif

      // 2e. Project limited differences in char. variables 
      //     on right eigenvectors to obtain the corresponding 
      //     primitive quantities,  dv = \sum dw.R

      for (nv = 0; nv < NFLX; nv++) {
        dp = dm = 0.0;
        for (int l = 0; l < NFLX; l++){
          dp += dwp[l]*R[nv][l];
          dm += dwm[l]*R[nv][l];
        }
        dvp[nv] = dp;
        dvm[nv] = dm;
      }
      #if NVAR != NFLX
      for (nv = NFLX; nv < NVAR; nv++){
        dvp[nv] = dwp[nv];
        dvm[nv] = dwm[nv];
      }
      #endif

      // 2f. Build L/R states in prim. variables and apply 
      //     parabolic limiter to primitive variables if 
      //     required.

      NVAR_LOOP(nv) {
        #if PARABOLIC_LIM == 1 || PARABOLIC_LIM == 2
        if (dvp[nv]*dvm[nv] >= 0.0) dvp[nv] = dvm[nv] = 0.0;
        else {
          if      (fabs(dvp[nv]) >= cm*fabs(dvm[nv])) dvp[nv] = -cm*dvm[nv];
          else if (fabs(dvm[nv]) >= cp*fabs(dvp[nv])) dvm[nv] = -cp*dvp[nv];
        }
        #endif
        vL(i+1,nv) = v(i,nv) + dvp[nv];
        vR(i,nv)   = v(i,nv) + dvm[nv];
      }

      // 2g. Make sure that left and right states are
      //     physically admissible. If not, use linear 
      //     reconstruction on density and pressure.

      if (vL(i+1,RHO) < 0.0 || vR(i,RHO) < 0.0) {
        dvp[RHO] = 0.5*(MINMOD_LIMITER(dvF_p[RHO], dvF_m[RHO]));
        dvm[RHO] = - dvp[RHO];
        vL(i+1,RHO) = v(i,RHO) + dvp[RHO];
        vR(i,RHO)   = v(i,RHO) + dvm[RHO];
      }
      #if HAVE_ENERGY
      if (vL(i+1,PRS) < 0.0 || vR(i,PRS) < 0.0) {
        dvp[PRS] = 0.5*(MINMOD_LIMITER(dvF_p[PRS], dvF_m[PRS]));
        dvm[PRS] = - dvp[PRS];
        vL(i+1,PRS) = v(i,PRS) + dvp[PRS];
        vR(i,PRS)   = v(i,PRS) + dvm[PRS];
      }
      #endif
      #if RADIATION
      if (vL(i+1,ENR) < 0.0 || vR(i,ENR) < 0.0) {
        dvp[ENR] = 0.5*(MINMOD_LIMITER(dvF_p[ENR], dvF_m[ENR]));
        dvm[ENR] = - dvp[ENR];
        vL(i+1,ENR) = v(i,ENR) + dvp[ENR];
        vR(i,ENR)   = v(i,ENR) + dvm[ENR];
      }
      #endif
    } // End block full reconstruction 
          
    #if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD)
    VelocityLimiter (v, vL, vR, i);
    #endif
    #if RADIATION
    RadFluxLimFlatten (v, vL, vR, i);
    #endif

   // 2h. Assign face-centered magnetic field

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
  }}}

}
#undef PARABOLIC_LIM
#endif  // CHAR_LIMITING == YES
