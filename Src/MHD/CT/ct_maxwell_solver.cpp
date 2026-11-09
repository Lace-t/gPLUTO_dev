////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Maxwell solver for the Resistive RMHD equations.

  \author  A. Mignone (andrea.mignone@unito.it)\n
           V. Berta   (vittoria.berta@unito.it)
           
  \date    Nov 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if PHYSICS == ResRMHD || PHYSICS == GRMHD
////////////////////////////////////////////////////////////////////////
void CT_MaxwellSolver (const Data *data, const EMF *emf, Grid *grid, RBox *box)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
#if DIMENSIONS > 1

  int i, j, k;
  int ni, nj;

#ifdef HIGH_ORDER
  #if HO_AVERAGE < 3
  int ns = 0;
  #else
  int ns = 1;
  #endif
#endif

  int recE = RECONSTRUCTION;
  int recB = RECONSTRUCTION;

#ifdef HIGH_ORDER 
  const Ary3D *Vs = data->Vps;
#else
  const Ary3D *Vs = data->Vs;
#endif

  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  int ibeg = box->ibeg, iend = box->iend+INCLUDE_IDIR;
  int jbeg = box->jbeg, jend = box->jend+INCLUDE_JDIR;
  int kbeg = box->kbeg, kend = box->kend+INCLUDE_KDIR;
  long int size_tot = nx3_tot*nx2_tot*nx1_tot;

  auto &Ex1e = emf->Ee[IDIR];
  auto &Ex2e = emf->Ee[JDIR];
  auto &Ex3e = emf->Ee[KDIR];

  auto &Bx1e = emf->Be[IDIR];
  auto &Bx2e = emf->Be[JDIR];
  auto &Bx3e = emf->Be[KDIR];

  RBox box_rec;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);

  // ------------------------------------------------------
  // 0. Lower reconstruction order 
  //    if failsafe is enabled
  // ------------------------------------------------------

  #if FALLBACK == YES
  if (data->integration_mode == FALLBACK) {
    recE = LINEAR;
    recB = LINEAR;
  }
  #endif



  mynvtxstart_("CT_MaxwellSolver",YELLOW);

#if INCLUDE_IDIR

  ARRAY_STRIDES (grid, IDIR, ni, nj);

  // ------------------------------------------------------
  // X1a. Reconstruct along the x1-direction:
  //
  //      o Ez(j+1/2) and By(j+1/2) --> (i+1/2, j+1/2, k)
  //      o Bz(j+1/2) and Ey(j+1/2) --> (i+1/2, j+1/2, k)
  //      o Initialize Ex3 and Bx3 at (i+1/2, j+1/2, k)
  // ------------------------------------------------------
  
  #if INCLUDE_JDIR
  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                     iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg,                     jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg - ns - INCLUDE_KDIR, kend = box->kend + INCLUDE_KDIR + ns;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ArrayReconstruct<IDIR> (data, emf->e[KDIR][JDIR], JDIR, recE, grid, &box_rec, 0); 
  ArrayReconstruct<IDIR> (data, Vs[BX2s],           JDIR, recB, grid, &box_rec, 1); 
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  ArrayReconstruct<IDIR> (data, emf->B[KDIR][JDIR], JDIR, recE, grid, &box_rec, 2); 
  ArrayReconstruct<IDIR> (data, Vs[EX2s],           JDIR, recB, grid, &box_rec, 3); 
  #endif

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf, grid, data) 
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    long int offset  = ni*(j + nj*k);
    long int offset1 = ni*(j + nj*k) + size_tot;
    Ary1D EzL(&data->sweep->vL[offset],ni);
    Ary1D EzR(&data->sweep->vR[offset],ni);
    Ary1D ByL(&data->sweep->vL[offset1],ni);
    Ary1D ByR(&data->sweep->vR[offset1],ni);
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    long int offset2 = ni*(j + nj*k) + size_tot*2;
    long int offset3 = ni*(j + nj*k) + size_tot*3;
    Ary1D BzL(&data->sweep->vL[offset2],ni);
    Ary1D BzR(&data->sweep->vR[offset2],ni);
    Ary1D EyL(&data->sweep->vL[offset3],ni);
    Ary1D EyR(&data->sweep->vR[offset3],ni);
    #endif

    Ex3e(i,j,k) = 0.25*(EzL(i) + EzR(i)) + 0.5*(ByR(i) - ByL(i));
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    Bx3e(i,j,k) = 0.25*(BzL(i) + BzR(i)) - 0.5*(EyR(i) - EyL(i));
    #endif

    #if CT_ENERGY_CORRECTION > 1
    emf->Bx2ez(i,j,k)    = 0.5*(ByL(i) + ByR(i));
    #endif
  }}}
  #endif  // INCLUDE_JDIR

  // ------------------------------------------------------
  // X1b. Reconstruct along the x1-direction:
  //
  //      o Ey(k+1/2) and Bz(k+1/2) --> (i+1/2, j, k+1/2)
  //      o By(k+1/2) and Ez(k+1/2) --> (i+1/2, j, k+1/2)
  //      o Initialize Ex2 and Bx2 at (i+1/2, j, k+1/2)
  // ------------------------------------------------------

  #if INCLUDE_KDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                     iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg - ns - INCLUDE_JDIR, jend = box->jend + INCLUDE_JDIR + ns;
  kbeg = box->kbeg,                     kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ArrayReconstruct<IDIR> (data, emf->e[JDIR][KDIR], KDIR, recE, grid, &box_rec, 0); 
  ArrayReconstruct<IDIR> (data, Vs[BX3s],           KDIR, recB, grid, &box_rec, 1); 
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  ArrayReconstruct<IDIR> (data, emf->B[JDIR][KDIR], KDIR, recE, grid, &box_rec, 2); 
  ArrayReconstruct<IDIR> (data, Vs[EX3s],           KDIR, recB, grid, &box_rec, 3); 
  #endif

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf, grid, data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    long int offset = ni*(j + nj*k);
    long int offset1 = ni*(j + nj*k) + size_tot;
    Ary1D EyL(&data->sweep->vL[offset],ni);
    Ary1D EyR(&data->sweep->vR[offset],ni);
    Ary1D BzL(&data->sweep->vL[offset1],ni);
    Ary1D BzR(&data->sweep->vR[offset1],ni);
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    long int offset2 = ni*(j + nj*k) + size_tot*2;
    long int offset3 = ni*(j + nj*k) + size_tot*3;
    Ary1D ByL(&data->sweep->vL[offset2],ni);
    Ary1D ByR(&data->sweep->vR[offset2],ni);
    Ary1D EzL(&data->sweep->vL[offset3],ni);
    Ary1D EzR(&data->sweep->vR[offset3],ni);
    #endif

    Ex2e(i,j,k) = 0.25*(EyL(i) + EyR(i)) - 0.5*(BzR(i) - BzL(i));
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    Bx2e(i,j,k) = 0.25*(ByL(i) + ByR(i)) + 0.5*(EzR(i) - EzL(i));
    #endif
  }}}
  #endif  // INCLUDE_KDIR 

#endif // INCLUDE_IDIR 

#if INCLUDE_JDIR
  
  ARRAY_STRIDES (grid, JDIR, ni, nj);

  // ------------------------------------------------------
  // X2a. Reconstruct along the x2-direction:
  //
  //      o Ex(k+1/2) and Bz(k+1/2) --> (i, j+1/2, k+1/2)
  //      o Bx(k+1/2) and Ez(k+1/2) --> (i, j+1/2, k+1/2)
  //      o Initialie Ex1 and Bx1 at (i, j+1/2, k+1/2)
  //
  //      o Ez(i+1/2) and Bx(i+1/2) --> (i+1/2, j+1/2, k)
  //      o Bz(i+1/2) and Ex(i+1/2) --> (i+1/2, j+1/2, k)
  // ------------------------------------------------------

  #if INCLUDE_KDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg - ns - INCLUDE_IDIR, iend = box->iend + INCLUDE_IDIR + ns;
  jbeg = box->jbeg,                     jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg,                     kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ArrayReconstruct<JDIR> (data, emf->e[IDIR][KDIR], KDIR, recE, grid, &box_rec, 0); 
  ArrayReconstruct<JDIR> (data, Vs[BX3s],           KDIR, recB, grid, &box_rec, 1); 
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  ArrayReconstruct<JDIR> (data, emf->B[IDIR][KDIR], KDIR, recE, grid, &box_rec, 2); 
  ArrayReconstruct<JDIR> (data, Vs[EX3s],           KDIR, recB, grid, &box_rec, 3); 
  #endif

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf, grid, data) 
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
  for (j = jbeg; j <= jend; j++){
    long int offset = ni*(i + nj*k);
    long int offset1 = ni*(i + nj*k) + size_tot;
    Ary1D ExL(&data->sweep->vL[offset],ni);
    Ary1D ExR(&data->sweep->vR[offset],ni);
    Ary1D BzL(&data->sweep->vL[offset1],ni);
    Ary1D BzR(&data->sweep->vR[offset1],ni);
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    long int offset2 = ni*(i + nj*k) + size_tot*2;
    long int offset3 = ni*(i + nj*k) + size_tot*3;
    Ary1D BxL(&data->sweep->vL[offset2],ni);
    Ary1D BxR(&data->sweep->vR[offset2],ni);
    Ary1D EzL(&data->sweep->vL[offset3],ni);
    Ary1D EzR(&data->sweep->vR[offset3],ni);
    #endif
    
    Ex1e(i,j,k) = 0.25*(ExL(j) + ExR(j)) + 0.5*(BzR(j) - BzL(j));
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT  
    Bx1e(i,j,k) = 0.25*(BxL(j) + BxR(j)) - 0.5*(EzR(j) - EzL(j));
    #endif
  }}}
  #endif // INCLUDE_KDIR 

  // ------------------------------------------------------
  // X2b. Reconstruct along the x2-direction:
  //
  //      o Ez(i+1/2) and Bx(i+1/2) --> (i+1/2, j+1/2, k)
  //      o Bz(i+1/2) and Ex(i+1/2) --> (i+1/2, j+1/2, k)
  //      o Complete Ex3 and Bx3 at (i+1/2,j+1/2, k)
  // ------------------------------------------------------

  #if INCLUDE_IDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                     iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg,                     jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg - ns - INCLUDE_KDIR, kend = box->kend + INCLUDE_KDIR + ns;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ArrayReconstruct<JDIR> (data, emf->e[KDIR][IDIR], IDIR, recE, grid, &box_rec, 0); 
  ArrayReconstruct<JDIR> (data, Vs[BX1s],           IDIR, recB, grid, &box_rec, 1); 
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  ArrayReconstruct<JDIR> (data, emf->B[KDIR][IDIR], IDIR, recE, grid, &box_rec, 2); 
  ArrayReconstruct<JDIR> (data, Vs[EX1s],           IDIR, recB, grid, &box_rec, 3); 
  #endif

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf, grid, data) 
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
  for (j = jbeg; j <= jend; j++){  
    long int offset = ni*(i + nj*k);
    long int offset1 = ni*(i + nj*k) + size_tot;
    Ary1D EzL(&data->sweep->vL[offset],ni);
    Ary1D EzR(&data->sweep->vR[offset],ni);
    Ary1D BxL(&data->sweep->vL[offset1],ni);
    Ary1D BxR(&data->sweep->vR[offset1],ni);
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    long int offset2 = ni*(i + nj*k) + size_tot*2;
    long int offset3 = ni*(i + nj*k) + size_tot*3;
    Ary1D BzL(&data->sweep->vL[offset2],ni);
    Ary1D BzR(&data->sweep->vR[offset2],ni);
    Ary1D ExL(&data->sweep->vL[offset3],ni);
    Ary1D ExR(&data->sweep->vR[offset3],ni);
    #endif
    
    Ex3e(i,j,k) += 0.25*(EzL(j) + EzR(j)) - 0.5*(BxR(j) - BxL(j));
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    Bx3e(i,j,k) += 0.25*(BzL(j) + BzR(j)) + 0.5*(ExR(j) - ExL(j));
    #endif

    #if CT_ENERGY_CORRECTION > 1
    emf->Bx1ez(i,j,k) = 0.5*(BxL(j) + BxR(j));
    #endif
  }}}
  #endif // INCLUDE_IDIR

#endif  // INCLUDE_JDIR

#if INCLUDE_KDIR
  
  ARRAY_STRIDES (grid, KDIR, ni, nj);

  // ------------------------------------------------------
  // X3a. Reconstruct along the x3-direction:
  //
  //      o Ey(i+1/2) and Bx(i+1/2) --> (i+1/2, j, k+1/2)
  //      o By(i+1/2) and Ex(i+1/2) --> (i+1/2, j, k+1/2)
  //      o Complete Ex2 at (i+1/2,j,k+1/2)
  //
  //      o Ex(j+1/2) and By(j+1/2) --> (i, j+1/2, k+1/2)
  //      o Bx(j+1/2) and Ey(j+1/2) --> (i, j+1/2, k+1/2)
  //  -------------------------------------------------------- 

  #if INCLUDE_IDIR
  
  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                     iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg - ns - INCLUDE_JDIR, jend = box->jend + INCLUDE_JDIR + ns;
  kbeg = box->kbeg,                     kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ArrayReconstruct<KDIR> (data, emf->e[JDIR][IDIR], IDIR, recE, grid, &box_rec, 0); 
  ArrayReconstruct<KDIR> (data, Vs[BX1s],           IDIR, recB, grid, &box_rec, 1); 
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  ArrayReconstruct<KDIR> (data, emf->B[JDIR][IDIR], IDIR, recE, grid, &box_rec, 2); 
  ArrayReconstruct<KDIR> (data, Vs[EX1s],           IDIR, recB, grid, &box_rec, 3); 
  #endif

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf, grid, data)
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  for (k = kbeg; k <= kend; k++){
  
    long int offset = ni*(i + nj*j);
    long int offset1 = ni*(i + nj*j) + size_tot;
    Ary1D EyL(&data->sweep->vL[offset],ni);
    Ary1D EyR(&data->sweep->vR[offset],ni);
    Ary1D BxL(&data->sweep->vL[offset1],ni);
    Ary1D BxR(&data->sweep->vR[offset1],ni);
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    long int offset2 = ni*(i + nj*j) + size_tot*2;
    long int offset3 = ni*(i + nj*j) + size_tot*3;
    Ary1D ByL(&data->sweep->vL[offset2],ni);
    Ary1D ByR(&data->sweep->vR[offset2],ni);
    Ary1D ExL(&data->sweep->vL[offset3],ni);
    Ary1D ExR(&data->sweep->vR[offset3],ni);
    #endif
    
    Ex2e(i,j,k) += 0.25*(EyL(k) + EyR(k)) + 0.5*(BxR(k) - BxL(k));
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    Bx2e(i,j,k) += 0.25*(ByL(k) + ByR(k)) - 0.5*(ExR(k) - ExL(k));
    #endif
  }}}
  #endif // INCLUDE_IDIR

  // ------------------------------------------------------
  // X3b. Reconstruct along the x3-direction:
  //
  //      o Ex(j+1/2) and By(j+1/2) --> (i, j+1/2, k+1/2)
  //      o Bx(j+1/2) and Ey(j+1/2) --> (i, j+1/2, k+1/2)
  //      o Complete Ex1 at (i, j+1/2, k+1/2)
  // ------------------------------------------------------

  #if INCLUDE_JDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg - ns - INCLUDE_IDIR, iend = box->iend + INCLUDE_IDIR + ns;
  jbeg = box->jbeg,                     jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg,                     kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ArrayReconstruct<KDIR> (data, emf->e[IDIR][JDIR], JDIR, recE, grid, &box_rec, 0); 
  ArrayReconstruct<KDIR> (data, Vs[BX2s],           JDIR, recB, grid, &box_rec, 1); 
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  ArrayReconstruct<KDIR> (data, emf->B[IDIR][JDIR], JDIR, recE, grid, &box_rec, 2); 
  ArrayReconstruct<KDIR> (data, Vs[EX2s],           JDIR, recB, grid, &box_rec, 3); 
  #endif

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf, grid, data)
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  for (k = kbeg; k <= kend; k++){   
    long int offset = ni*(i + nj*j);
    long int offset1 = ni*(i + nj*j) + size_tot;
    Ary1D ExL(&data->sweep->vL[offset],ni);
    Ary1D ExR(&data->sweep->vR[offset],ni);
    Ary1D ByL(&data->sweep->vL[offset1],ni);
    Ary1D ByR(&data->sweep->vR[offset1],ni);
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    long int offset2 = ni*(i + nj*j) + size_tot*2;
    long int offset3 = ni*(i + nj*j) + size_tot*3;
    Ary1D BxL(&data->sweep->vL[offset2],ni);
    Ary1D BxR(&data->sweep->vR[offset2],ni);
    Ary1D EyL(&data->sweep->vL[offset3],ni);
    Ary1D EyR(&data->sweep->vR[offset3],ni);
    #endif
 
    Ex1e(i,j,k) += 0.25*(ExL(k) + ExR(k)) - 0.5*(ByR(k) - ByL(k));
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    Bx1e(i,j,k) += 0.25*(BxL(k) + BxR(k)) + 0.5*(EyR(k) - EyL(k));
    #endif
  }}}
  #endif // INCLUDE_JDIR

#endif // INCLUDE_KDIR

  mynvtxstop_();

#endif
}
#endif // PHYSICS == ResRMHD 
