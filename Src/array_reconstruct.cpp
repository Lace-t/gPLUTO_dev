////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Reconstruct an array Q(i,j,k) along a given direction, 
  defining left and right values at the interfaces qL and qR.
  qL and qR are unique memory portion of \c data->sweep->vL and 
  \c data->sweep->vR respectively, where the offset is determined by the \c
  variable nv.
  In fact, all consecutive calls to ArrayReconstruct<>() have different 
  values of nv. 
  In this way we recycle the structures sweep used in Reconstruct().

  \authors A. Mignone  (andrea.mignone@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it)
  \date    Apr 20, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "States/states.hpp"

// Explicit template instantiations
template void ArrayReconstruct<IDIR> (const Data *, const Ary3D &, int,
                                      int, Grid *, RBox *, const int);
template void ArrayReconstruct<JDIR> (const Data *, const Ary3D &, int,
                                      int, Grid *, RBox *, const int);
template void ArrayReconstruct<KDIR> (const Data *, const Ary3D &, int,
                                      int, Grid *, RBox *, const int);

////////////////////////////////////////////////////////////////////////
/*!
  \param [in] data       Pointer to data structure
  \param [in] Q          The 3D Array used to be reconstructed
  \param [in] flag_dir   The variable displacement with respect to the 
                         center (e.g. JDIR means the Q is located at j-1/2).
                         Used to compute the flag in a consisten way.
  \param [in] lim        The slope limiter (only for LINEAR reconstruction)
  \param [in] grid       Pointer to the grid structure
  \param [in] box        The box over which reconstruction is needed
  \param [in] nv         The memory location used when indexing the 
  sweep array, allowing more arrays to be 
  reconstructed
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void ArrayReconstruct (const Data *data, const Ary3D &Q, int flag_dir,
                       int lim, Grid *grid, RBox *box, const int nv)
{
  int i,j,k,ni,nj;
  double dqp, dqm, dq;
  int rec0, rec;
  int ntot;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  ARRAY_STRIDES (grid, dir, ni, nj);

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];

  // 0a. Copy 3D array onto 1D array (x1-direction)

  if (dir == IDIR) {
    ntot = nx1_tot;
    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = 0; i < ntot; i++){
    
      long int offset1 = ni*(j + nj*k);
      long int offset  = offset1 + size_tot*nv;
      Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
      Ary1D q(&data->sweep->vC[offset],ni);
      
      q(i) = Q(i,j,k);

      if (flag_dir==IDIR) flag(i) = data->flag(i,j,k);
      if (flag_dir==JDIR) flag(i) = data->flag(i,j,k) | data->flag(i,j-1,k);
      if (flag_dir==KDIR) flag(i) = data->flag(i,j,k) | data->flag(i,j,k-1);
    }}}
  }
 
  // 0b. Copy 3D array onto 1D array (x2-direction)

  if (dir == JDIR) {
    ntot = nx2_tot;
    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data)
    
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = 0; j < ntot; j++){
    
      long int offset1 = ni*(i + nj*k);
      long int offset  = offset1 + size_tot*nv;
      Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
      Ary1D q(&data->sweep->vC[offset],ni);
     
      q(j) = Q(i,j,k);

      if (flag_dir==IDIR) flag(j) = data->flag(i,j,k) | data->flag(i-1,j,k);
      if (flag_dir==JDIR) flag(j) = data->flag(i,j,k);
      if (flag_dir==KDIR) flag(j) = data->flag(i,j,k) | data->flag(i,j,k-1);
    }}}
  }

  // 0c. Copy 3D array onto 1D array (x3-direction)

  if (dir == KDIR) {
    ntot = nx3_tot;
    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = 0; k < ntot; k++){

      long int offset1 = ni*(i + nj*j);
      long int offset  = offset1 + size_tot*nv;
      Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
      Ary1D q(&data->sweep->vC[offset],ni);
      
      q(k) = Q(i,j,k);
      
      if (flag_dir==IDIR) flag(k) = data->flag(i,j,k) | data->flag(i-1,j,k);
      if (flag_dir==JDIR) flag(k) = data->flag(i,j,k) | data->flag(i,j-1,k);
      if (flag_dir==KDIR) flag(k) = data->flag(i,j,k);
    }}}
  }

  // 1. Cycle on directional indexes to preserve sweep ordering
  //    In the following kernels, k,j and i correspond
  //    to binormal, tangent and normal directions         

  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);

  // 1a. Flat reconstruction
  
  rec  = lim;
  rec0 = (rec/10)*10;
  if (rec0 == FLAT) {
    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = 1; i < ntot-1; i++){
    
      long int offset = ni*(j + nj*k) + size_tot*nv;
      
      Ary1D q(&data->sweep->vC[offset],ni);
      Ary1D qL(&data->sweep->vL[offset],ni);
      Ary1D qR(&data->sweep->vR[offset],ni);

      qL(i+1) = q(i);
      qR(i)   = q(i);
    }}}
  }

  // 1b. Linear reconstruction
  
  rec  = lim;
  rec0 = (rec/10)*10;
  if (rec0 == LINEAR) {
    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = 1; i < ntot-1; i++){

      long int offset1 = ni*(j + nj*k);
      long int offset  = offset1 + size_tot*nv;

      Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
      Ary1D q(&data->sweep->vC[offset],ni);
      Ary1D qL(&data->sweep->vL[offset],ni);
      Ary1D qR(&data->sweep->vR[offset],ni);
          
      #if SHOCK_FLATTENING == MULTID
      if (flag(i) & FLAG_FLAT) {
        qL(i+1) = qR(i) = q(i);
        continue;
      }
      #endif

      dqp = q(i+1) - q(i);
      dqm = q(i)   - q(i-1);

      if      (rec == FLAT_LIM)    dq = 0.0;
      else if (rec == MINMOD_LIM)  dq = MINMOD_LIMITER(dqp, dqm);
      else if (rec == VANLEER_LIM) dq = VANLEER_LIMITER(dqp, dqm);
      else if (rec == MC_LIM)      dq = MC_LIMITER(dqp, dqm);
      else{
        SET_LIMITER (dq, dqp, dqm, 2.0, 2.0);
      }
      qL(i+1) = q(i) + 0.5*dq;
      qR(i)   = q(i) - 0.5*dq;   

    }}} // End loop (i,j,k)
  } // End if (rec0 == LINEAR)

  // 1c. PARABOLIC reconstruction

  #if RECONSTRUCTION == PARABOLIC
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = 1; i < ntot-1; i++){
    long int offset1 = ni*(j + nj*k);
    long int offset  = offset1 + size_tot*nv;

    Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
    Ary1D q(&data->sweep->vC[offset],ni);
    Ary1D qL(&data->sweep->vL[offset],ni);
    Ary1D qR(&data->sweep->vR[offset],ni);

    double qp, qm;
    
    #if SHOCK_FLATTENING == MULTID
    if (flag(i) & FLAG_FLAT) {
      qL(i+1) = qR(i) = q(i);
      continue;
    }else if (flag(i) & FLAG_MINMOD) {
      dqp = q(i+1) - q(i);
      dqm = q(i) - q(i-1);
      dq = MINMOD_LIMITER(dqp, dqm);
      qL(i+1) = q(i) + 0.5*dq;
      qR(i)   = q(i) - 0.5*dq;
      continue;
    }
    #endif

    qp = (    -q(i-1) + 5.0*q(i) + 2.0*q(i+1))/6.0;
    qm = ( 2.0*q(i-1) + 5.0*q(i) -     q(i+1))/6.0;
    
    dqp = qp - q(i);
    dqm = qm - q(i);

    dq  = q(i+1) - q(i);
    qp  = q(i) + MINMOD_LIMITER(dqp, dq);
    
    dq  = q(i) - q(i-1);
    qm  = q(i) + MINMOD_LIMITER(dqm, -dq);
    
    dqp = qp - q(i);
    dqm = qm - q(i);

    if (dqp*dqm >= 0.0) dqp = dqm = 0.0;
    else{
      if      (fabs(dqp) >= 2.0*fabs(dqm)) dqp = -2.0*dqm;
      else if (fabs(dqm) >= 2.0*fabs(dqp)) dqm = -2.0*dqp;
    }
    qL(i+1) = qp = q(i) + dqp; 
    qR(i)   = qm = q(i) + dqm;
  }}}
  #endif // RECONSTRUCTION == PARABOLIC

  // 1d. LIMO3 & WENO3 reconstruction (3 point stencil)

  #if (RECONSTRUCTION == LimO3) || (RECONSTRUCTION == WENO3)
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = 1; i < ntot-1; i++){

    long int offset1 = ni*(j + nj*k);
    long int offset  = offset1 + size_tot*nv;

    Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
    Ary1D q(&data->sweep->vC[offset],ni);
    Ary1D qL(&data->sweep->vL[offset],ni);
    Ary1D qR(&data->sweep->vR[offset],ni);
  
    #if SHOCK_FLATTENING == MULTID
    if (flag(i) & FLAG_FLAT) {
      qL(i+1) = q(i);
      qR(i)   = q(i);
      continue;
    }else if (flag(i) & FLAG_MINMOD) {
      dqp = q(i+1) - q(i);
      dqm = q(i) - q(i-1);
      dq = MINMOD_LIMITER(dqp, dqm);
      qL(i+1) = q(i) + 0.5*dq;
      qR(i)   = q(i) - 0.5*dq;
      continue;
    }
    #endif

    #if HO_ORDER_REDUCTION == LINEAR
    if (flag(i) & FLAG_HO_LAP_LIMITER){
      dqp = q(i+1) - q(i);
      dqm = q(i) - q(i-1);
      dq  = VANLEER_LIMITER(dqp, dqm);
      qL(i+1) = q(i) + 0.5*dq;
      qR(i)   = q(i) - 0.5*dq;
    }else
    #endif
    {
      #if RECONSTRUCTION == LimO3
      qL(i+1) = LimO3_States (&q(i), +1, grid->dx[dir][i]);
      qR(i)   = LimO3_States (&q(i), -1, grid->dx[dir][i]);
      #elif RECONSTRUCTION == WENO3
      qL(i+1) = WENO3_States (&q(i), +1, grid->dx[dir][i]);
      qR(i)   = WENO3_States (&q(i), -1, grid->dx[dir][i]);
      #endif
    }   
  }}}
  #endif // (RECONSTRUCTION == LimO3) || (RECONSTRUCTION == WENO3)

  // 1e. WENOZ & MP5 reconstructions (5 point stencil)

  #if (RECONSTRUCTION == WENOZ) || (RECONSTRUCTION == MP5)
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = 2; i < ntot-2; i++){
    long int offset1 = ni*(j + nj*k);
    long int offset  = offset1 + size_tot*nv;

    Ary1Duint32_t flag(&data->sweep->flag[offset1],ni);
    Ary1D q(&data->sweep->vC[offset],ni);
    Ary1D qL(&data->sweep->vL[offset],ni);
    Ary1D qR(&data->sweep->vR[offset],ni);

    #if SHOCK_FLATTENING == MULTID
    if (flag(i) & FLAG_FLAT) {
      qL(i+1) = q(i);
      qR(i)   = q(i);
      continue;
    }else if (flag(i) & FLAG_MINMOD) {
      dqp = q(i+1) - q(i);
      dqm = q(i) - q(i-1);
      dq = MINMOD_LIMITER(dqp, dqm);
      qL(i+1) = q(i) + 0.5*dq;
      qR(i)   = q(i) - 0.5*dq;
      continue;
    }
    #endif

    #if HO_ORDER_REDUCTION == WENO3 
    if (flag(i) & FLAG_HO_LAP_LIMITER){
      qL(i+1) = WENO3_States (&q(i), +1, grid->dx[dir][i]);
      qR(i)   = WENO3_States (&q(i), -1, grid->dx[dir][i]);
    }else
    #elif HO_ORDER_REDUCTION == LINEAR
    if (flag(i) & FLAG_HO_LAP_LIMITER){
      dqp = q(i+1) - q(i);
      dqm = q(i) - q(i-1);
      dq  = VANLEER_LIMITER(dqp, dqm);
      qL(i+1) = q(i) + 0.5*dq;
      qR(i)   = q(i) - 0.5*dq;
    }else
    #endif
    {
      #if RECONSTRUCTION == WENOZ
      qL(i+1) = WENOZ_States (&q(i), +1);
      qR(i)   = WENOZ_States (&q(i), -1);
      #elif RECONSTRUCTION == MP5
      qL(i+1) = MP5_States (&q(i), +1);
      qR(i)   = MP5_States (&q(i), -1); 
      #endif
    }   
  }}}
  #endif // (RECONSTRUCTION == WENOZ) || (RECONSTRUCTION == MP5)
  
  return;
}
