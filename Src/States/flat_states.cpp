////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Piecewise flat reconstruction.
      
  Note: this function is replicated twice since it is needed during
  fallback approaches

  \author A. Mignone  (andrea.mignone@unito.it)
          M. Rossazza (marco.rossazza@unito.it)
  \date   Dec 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if RECONSTRUCTION == FLAT
// Explicit template instantiations
template void Reconstruct<IDIR> (Data *, Grid *, RBox *);
template void Reconstruct<JDIR> (Data *, Grid *, RBox *);
template void Reconstruct<KDIR> (Data *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
template<int dir>
void Reconstruct (Data *data, Grid *grid, RBox *box)
/*!
*/
////////////////////////////////////////////////////////////////////////
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

  // 1a. Reflect / symmetrized profiles when zones are 
  //     flagged with INTERNAL_BOUNDARY.

  #if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
  #endif

  // ------------------------------------------------------
  // 2. Main Reconstruction Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  for (nv = 0; nv < NVAR; nv++){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    long int offset = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset],ni,NVAR);

    vL(i+1,nv) = v(i,nv);
    vR(i,nv)   = v(i,nv);
  }}}}

  // ------------------------------------------------------
  // 3. Adjustment Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){ 
    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D En(&data->sweep->En[offset],ni);

    // 3a. Assign face-centered magnetic field
       
    // #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    // vL(i,BXn) = vR(i,BXn) = Bn(i);
    // #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    // vL(i,EXn) = vR(i,EXn) = En(i);
    // #endif
    // #endif

    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    CT_ASSIGN(vL(i+1,BXn) = Bn(i+1);  
              vR(i,BXn)   = Bn(i);    ,
              vL(i+1,EXn) = En(i+1);
              vR(i,EXn)   = En(i);)
    #endif

  }}}
}
#endif // RECONSTRUCTION == FLAT


// Explicit template instantiations
template void ReconstructFlat<IDIR> (Data *, Grid *, RBox *);
template void ReconstructFlat<JDIR> (Data *, Grid *, RBox *);
template void ReconstructFlat<KDIR> (Data *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
template<int dir>
void ReconstructFlat (Data *data, Grid *grid, RBox *box)
/*!
*/
////////////////////////////////////////////////////////////////////////
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

  // 1a. Reflect / symmetrized profiles when zones are 
  //     flagged with INTERNAL_BOUNDARY.

  #if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
  InternalBoundaryReflect<dir>(data, grid, box);
  #endif

  // ------------------------------------------------------
  // 2. Main Reconstruction Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  for (nv = 0; nv < NVAR; nv++){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    long int offset = NVAR*ni*(j + nj*k);
    Ary2D v(&data->sweep->vC[offset],ni,NVAR);
    Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset],ni,NVAR);

    vL(i+1,nv) = v(i,nv);
    vR(i,nv)   = v(i,nv);
  }}}}

  // ------------------------------------------------------
  // 3. Adjustment Loop
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){ 
    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D En(&data->sweep->En[offset],ni);

    // 3a. Assign face-centered magnetic field
       
    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    CT_ASSIGN(vL(i+1,BXn) = Bn(i+1);  
              vR(i,BXn)   = Bn(i);    ,
              vL(i+1,EXn) = En(i+1);
              vR(i,EXn)   = En(i);)
    #endif

  }}}
}
