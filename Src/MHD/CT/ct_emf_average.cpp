////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Collects different EMF averaging schemes.

  \author  A. Mignone   (andrea.mignone@unito.it)\n
           M. Rossazza  (marco.rossazza@unito.it)
  \date    Nov 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Compute arithmetic average of EMF at cell edges by combining
  the four electric field values computed at zone faces
  as upwind Godunov fluxes into an edge-centered value.\n
  The face-centered EMF should have been stored by previous calls
  to CT_StoreEMF() during the one-dimensional sweeps.\n
  This function employs a simple arithmetic averaging of the
  face-centered electric field.
  
  \b References:
  - "A Staggered Mesh Algorithm Using High Order Godunov Fluxes to
  Ensure Solenoidal Magnetic Fields in Magnetohydrodynamic
  Simulations"\n
  Balsara \& Spicer, JCP (1999) 149, 270.
  
  \param [in]      Z1    pointer to EMF structure
  \param [in]      w     weighting factor
  \param [in]  Grid     w     weighting factor
  
  \return  This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
void CT_EMF_ArithmeticAverage (const EMF *Z1, const double w, RBox *box)
{
#if DIMENSIONS > 1

  int i, j, k;
  int di = 0, dj = 0, dk = 0;  // Extension zones for HIGH_ORDER 
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  #ifdef HIGH_ORDER 
  di = INCLUDE_IDIR;
  dj = INCLUDE_JDIR;
  dk = INCLUDE_KDIR;
  #endif

  #if CT_ENERGY_CORRECTION > 1
  QuitPLUTO("Conservative CT energy correction unavailable for Arithmetic EMF ",1);
  #endif

  mynvtxstart_("CT_EMF_ArithmeticAverage",YELLOW);

  #if INCLUDE_JDIR && INCLUDE_KDIR
  auto &Ex1e = Z1->Ee[IDIR];
  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3)
  for (k = kbeg;    k <= kend+1;  k++){
  for (j = jbeg;    j <= jend+1;  j++){
  for (i = ibeg-di; i <= iend+di; i++){
    Ex1e(i,j,k) = w*(  Z1->e[IDIR][KDIR](i,j,k) + Z1->e[IDIR][KDIR](i,j-1,k)
                     + Z1->e[IDIR][JDIR](i,j,k) + Z1->e[IDIR][JDIR](i,j,k-1));
  }}}
  #endif

  #if INCLUDE_IDIR && INCLUDE_KDIR
  auto &Ex2e = Z1->Ee[JDIR];
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3)
  for (k = kbeg;    k <= kend+1;  k++){
  for (j = jbeg-dj; j <= jend+dj; j++){
  for (i = ibeg;    i <= iend+1;  i++){
    Ex2e(i,j,k) = w*(  Z1->e[JDIR][IDIR](i,j,k) + Z1->e[JDIR][IDIR](i,j,k-1)
                     + Z1->e[JDIR][KDIR](i,j,k) + Z1->e[JDIR][KDIR](i-1,j,k));
  }}}
  #endif

  #if INCLUDE_IDIR && INCLUDE_JDIR
  auto &Ex3e = Z1->Ee[KDIR];
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3)
  for (k = kbeg-dk; k <= kend+dk; k++){
  for (j = jbeg;    j <= jend+1;  j++){
  for (i = ibeg;    i <= iend+1;  i++){
    Ex3e(i,j,k) = w*(  Z1->e[KDIR][IDIR](i,j,k) + Z1->e[KDIR][IDIR](i,j-1,k)
                     + Z1->e[KDIR][JDIR](i,j,k) + Z1->e[KDIR][JDIR](i-1,j,k));
  }}}
  #endif
  
  mynvtxstop_();

#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Add derivatives to the 4-point arithmetic average of magnetic fields.
  Obtain the electric field at corners.
  
  \b References:
  - "An unsplit Godunov method for ideal MHD via constrained transport"\n
  Gardiner & Stone, JCP (2005) 205, 509.
  See Eq. (41), (45) and (50).
  
  \param [in]      d     pointer to PLUTO Data structure
  \param [in]      emf   pointer to EMF structure
  \param [in]      grid  pointer to Grid structure
  
  \return  This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
void CT_EMF_IntegrateToCorner (const Data *data, const EMF *emf, RBox *box)

#define DEX_DYP(k,j,i) (emf->e[IDIR][JDIR](i,j+1,k) - data->Ex1(i,j,k))
#define DEX_DZP(k,j,i) (emf->e[IDIR][KDIR](i,j,k+1) - data->Ex1(i,j,k))

#define DEY_DXP(k,j,i) (emf->e[JDIR][IDIR](i+1,j,k) - data->Ex2(i,j,k))
#define DEY_DZP(k,j,i) (emf->e[JDIR][KDIR](i,j,k+1) - data->Ex2(i,j,k))

#define DEZ_DXP(k,j,i) (emf->e[KDIR][IDIR](i+1,j,k) - data->Ex3(i,j,k))
#define DEZ_DYP(k,j,i) (emf->e[KDIR][JDIR](i,j+1,k) - data->Ex3(i,j,k))

#define DEX_DYM(k,j,i) (data->Ex1(i,j,k) - emf->e[IDIR][JDIR](i,j,k))
#define DEX_DZM(k,j,i) (data->Ex1(i,j,k) - emf->e[IDIR][KDIR](i,j,k))

#define DEY_DXM(k,j,i) (data->Ex2(i,j,k) - emf->e[JDIR][IDIR](i,j,k))
#define DEY_DZM(k,j,i) (data->Ex2(i,j,k) - emf->e[JDIR][KDIR](i,j,k))

#define DEZ_DXM(k,j,i) (data->Ex3(i,j,k) - emf->e[KDIR][IDIR](i,j,k))
#define DEZ_DYM(k,j,i) (data->Ex3(i,j,k) - emf->e[KDIR][JDIR](i,j,k))
{
#if DIMENSIONS > 1
  int    i, j, k;
  signed char  sx, sy, sz;
  double wS, wN, wE, wW;

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  #if CT_ENERGY_CORRECTION > 1
  QuitPLUTO("Conservative CT energy correction unavailable for CT_Contact EMF ",1);
  #endif

  // --------------------------------------------
  // 1. Compute Ex1e
  // --------------------------------------------

  #if INCLUDE_JDIR && INCLUDE_KDIR
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, emf)
  for (k = kbeg; k <= kend + INCLUDE_KDIR; k++){
  for (j = jbeg; j <= jend + INCLUDE_JDIR; j++){
  for (i = ibeg; i <= iend; i++){

    // 1a. Span Y - Faces:    dEx/dz
    
    sy = emf->svy(i,j,k-1);
    sy = emf->sv[JDIR](i,j,k-1);

    if   (sy != 0) wS = sy > 0 ? 1.0:0.0;
    else           wS = 0.5;

    sy = emf->svy(i,j,k);
    sy = emf->sv[JDIR](i,j,k);

    if   (sy != 0) wN = sy > 0 ? 1.0:0.0;
    else           wN = 0.5;

    emf->Ee[IDIR](i,j,k)  +=   wS*DEX_DZP(k-1,j-1,i) + (1.0-wS)*DEX_DZP(k-1,j,i)
                             - wN*DEX_DZM(k,j-1,i)   - (1.0-wN)*DEX_DZM(k,j,i);  

    // 1b. Span Z - Faces:    dEx/dy

    sz = emf->svz(i,j-1,k);
    sz = emf->sv[KDIR](i,j-1,k);

    if   (sz != 0) wW = sz > 0 ? 1.0:0.0;
    else           wW = 0.5;

    sz = emf->svz(i,j,k);
    sz = emf->sv[KDIR](i,j,k);

    if   (sz != 0) wE = sz > 0 ? 1.0:0.0;
    else           wE = 0.5;

    emf->Ee[IDIR](i,j,k)  +=   wW*DEX_DYP(k-1,j-1,i) + (1.0-wW)*DEX_DYP(k,j-1,i)
                             - wE*DEX_DYM(k-1,j,i)   - (1.0-wE)*DEX_DYM(k,j,i);

  }}}
  #endif

  // --------------------------------------------
  // 2. Compute Ex2e
  // --------------------------------------------

  #if INCLUDE_IDIR && INCLUDE_KDIR
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, emf)
  for (k = kbeg; k <= kend + INCLUDE_KDIR; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend + INCLUDE_IDIR; i++){

    // 2a. Span Z - Faces:    dEy/dx

    sz = emf->svz(i-1,j,k);
    sz = emf->sv[KDIR](i-1,j,k);
    if   (sz != 0) wS = sz > 0 ? 1.0:0.0;
    else           wS = 0.5;

    sz = emf->svz(i,j,k);
    sz = emf->sv[KDIR](i,j,k);
    if   (sz != 0) wN = sz > 0 ? 1.0:0.0;
    else           wN = 0.5;

    emf->Ee[JDIR](i,j,k)  +=   wS*DEY_DXP(k-1,j,i-1) + (1.0-wS)*DEY_DXP(k,j,i-1)
                             - wN*DEY_DXM(k-1,j,i)   - (1.0-wN)*DEY_DXM(k,j,i);

    // 2b. Span X - Faces:    dEy/dz

    sx = emf->svx(i,j,k-1);
    sx = emf->sv[IDIR](i,j,k-1);

    if   (sx != 0) wW = sx > 0 ? 1.0:0.0;
    else           wW = 0.5;

    sx = emf->svx(i,j,k);
    sx = emf->sv[IDIR](i,j,k);

    if   (sx != 0) wE = sx > 0 ? 1.0:0.0;
    else           wE = 0.5;

    emf->Ee[JDIR](i,j,k)  +=   wW*DEY_DZP(k-1,j,i-1) + (1.0-wW)*DEY_DZP(k-1,j,i)
                             - wE*DEY_DZM(k,j,i-1)   - (1.0-wE)*DEY_DZM(k,j,i);

  }}}
  #endif

  // --------------------------------------------
  // 3. Compute Ex3e
  // --------------------------------------------

  #if INCLUDE_IDIR && INCLUDE_JDIR
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, emf)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend + INCLUDE_JDIR; j++){
  for (i = ibeg; i <= iend + INCLUDE_IDIR; i++){

    // 3a. Span X - Faces:    dEz/dy

    sx = emf->svx(i,j-1,k);
    sx = emf->sv[IDIR](i,j-1,k);

    if   (sx != 0) wS = sx > 0 ? 1.0:0.0;
    else           wS = 0.5;

    sx = emf->svx(i,j,k);
    sx = emf->sv[IDIR](i,j,k);

    if   (sx != 0) wN = sx > 0 ? 1.0:0.0;
    else           wN = 0.5;

    emf->Ee[KDIR](i,j,k)  +=   wS*DEZ_DYP(k,j-1,i-1) + (1.0-wS)*DEZ_DYP(k,j-1,i)
                            - wN*DEZ_DYM(k,j,i-1)   - (1.0-wN)*DEZ_DYM(k,j,i);

    // 3b. Span Y - Faces:    dEz/dx

    sy = emf->svy(i-1,j,k);
    sy = emf->sv[JDIR](i-1,j,k);
    if   (sy != 0) wW = sy > 0 ? 1.0:0.0;
    else           wW = 0.5;

    sy = emf->svy(i,j,k);
    sy = emf->sv[JDIR](i,j,k);

    if   (sy != 0) wE = sy > 0 ? 1.0:0.0;
    else           wE = 0.5;

    emf->Ee[KDIR](i,j,k)  +=   wW*DEZ_DXP(k,j-1,i-1) + (1.0-wW)*DEZ_DXP(k,j,i-1)
                             - wE*DEZ_DXM(k,j-1,i)   - (1.0-wE)*DEZ_DXM(k,j,i);

  }}}
  #endif

#endif
}
#undef DEX_DYP
#undef DEX_DZP
#undef DEY_DXP
#undef DEY_DZP
#undef DEZ_DXP
#undef DEZ_DYP
#undef DEX_DYM
#undef DEX_DZM
#undef DEY_DXM
#undef DEY_DZM
#undef DEZ_DXM
#undef DEZ_DYM

////////////////////////////////////////////////////////////////////////
/*!
  Reconstruct electric field at faces to the edges.
  The electric field interpolated from the 1D Riemann solver and
  defined as E = E(pnt) + 2*E(diff), that is, the sum of the smooth
  flux term plust \e twice the diffusive part.
*/
////////////////////////////////////////////////////////////////////////
void CT_EMF_Flux(const Data *data, const EMF *emf, Grid *grid, RBox *box)
{
#if DIMENSIONS > 1
  int i, j, k;
  int lim;
  int ni, nj;

  #if CT_ENERGY_CORRECTION > 1 
  QuitPLUTO("! Conservative CT energy correction unavailable for CT_EMF_Flux",1);
  #endif

#if RECONSTRUCTION == LINEAR
  #if LIMITER == DEFAULT
  lim = MC_LIM;
  #else
  lim = LIMITER;
  #endif
#else
  lim = RECONSTRUCTION;
#endif

  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];

  int ibeg = box->ibeg, iend = box->iend+INCLUDE_IDIR;
  int jbeg = box->jbeg, jend = box->jend+INCLUDE_JDIR;
  int kbeg = box->kbeg, kend = box->kend+INCLUDE_KDIR;

  #if (INCLUDE_KDIR) && (INCLUDE_JDIR)
  auto &Ex1e = emf->Ee[IDIR];
  #endif
  #if (INCLUDE_KDIR) && (INCLUDE_IDIR)
  auto &Ex2e = emf->Ee[JDIR];
  #endif
  #if (INCLUDE_JDIR) && (INCLUDE_IDIR)
  auto &Ex3e = emf->Ee[KDIR];
  #endif

  RBox box_rec;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);

  // ------------------------------------------------------
  // X1. Reconstruct along the x1-direction:
  //
  //     o Ezj(j-1/2) and Dy(j-1/2) --> (i-1/2, j-1/2, k)
  //     o Eyk(k-1/2) and Dz(k-1/2) --> (i-1/2, j, k-1/2)
  // ------------------------------------------------------
  
#if INCLUDE_IDIR
  ARRAY_STRIDES (grid, IDIR, ni, nj);

  #if INCLUDE_JDIR
  ArrayReconstruct<IDIR> (data, emf->e[KDIR][JDIR],     JDIR, lim, grid, &box_rec, 0); 
  ArrayReconstruct<IDIR> (data, emf->e_dff[KDIR][JDIR], JDIR, lim, grid, &box_rec, 1); 
  #endif

  #if INCLUDE_KDIR
  ArrayReconstruct<IDIR> (data, emf->e[JDIR][KDIR],     KDIR, lim, grid, &box_rec, 2); 
  ArrayReconstruct<IDIR> (data, emf->e_dff[JDIR][KDIR], KDIR, lim, grid, &box_rec, 3); 
  #endif

  #pragma acc parallel loop collapse(2) present(emf, grid)
  #pragma omp target teams distribute parallel for collapse(2)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){

    #if INCLUDE_JDIR
    long int offset = ni*(j + nj*k);
    long int offset1 = ni*(j + nj*k) + size_tot;
    Ary1D EzjL(&data->sweep->vL[offset],ni);
    Ary1D EzjR(&data->sweep->vR[offset],ni);
    Ary1D dzjL(&data->sweep->vL[offset1],ni);
    Ary1D dzjR(&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){
      Ex3e(i,j,k) = 0.25*(EzjL(i) + EzjR(i)) + 0.5*(dzjL(i) + dzjR(i));
    }
    #endif  // INCLUDE_JDIR 

    #if INCLUDE_KDIR
    long int offset2 = ni*(j + nj*k) + size_tot*2;
    long int offset3 = ni*(j + nj*k) + size_tot*3;
    Ary1D EykL(&data->sweep->vL[offset2],ni);
    Ary1D EykR(&data->sweep->vR[offset2],ni);
    Ary1D dykL(&data->sweep->vL[offset3],ni);
    Ary1D dykR(&data->sweep->vR[offset3],ni);

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){
      Ex2e(i,j,k) = 0.25*(EykL(i) + EykR(i)) + 0.5*(dykL(i) + dykR(i));
    }
    #endif  // INCLUDE_KDIR 

  }}
  #endif

  // ------------------------------------------------------
  // X2. Reconstruct along the x2-direction:
  //
  //     o Ezi(i-1/2) and Dx(i-1/2) --> (i-1/2, j-1/2, k)
  //     o Exk(k-1/2) and Dz(k-1/2) --> (i, j-1/2, k-1/2)
  // ------------------------------------------------------

  #if INCLUDE_JDIR

  ARRAY_STRIDES (grid, JDIR, ni, nj);

  #if INCLUDE_IDIR
  ArrayReconstruct<JDIR> (data, emf->e[KDIR][IDIR],     IDIR, lim, grid, &box_rec, 0); 
  ArrayReconstruct<JDIR> (data, emf->e_dff[KDIR][IDIR], IDIR, lim, grid, &box_rec, 1); 
  #endif

  #if INCLUDE_KDIR
  ArrayReconstruct<JDIR> (data, emf->e[IDIR][KDIR],     KDIR, lim, grid, &box_rec, 2); 
  ArrayReconstruct<JDIR> (data, emf->e_dff[IDIR][KDIR], KDIR, lim, grid, &box_rec, 3); 
  #endif
  
  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(emf, grid)
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
     
    #if INCLUDE_IDIR
    long int offset  = ni*(i + nj*k);
    long int offset1 = ni*(i + nj*k) + size_tot;
    Ary1D EziL(&data->sweep->vL[offset],ni);
    Ary1D EziR(&data->sweep->vR[offset],ni);
    Ary1D dziL(&data->sweep->vL[offset1],ni);
    Ary1D dziR(&data->sweep->vR[offset1],ni);
    #pragma acc loop vector
    for (j = jbeg; j <= jend; j++){
      Ex3e(i,j,k) += 0.25*(EziL(j) + EziR(j)) + 0.5*(dziL(j) + dziR(j));
    }
    #endif // INCLUDE_IDIR 

    #if INCLUDE_KDIR
    long int offset2 = ni*(i + nj*k) + size_tot*2;
    long int offset3 = ni*(i + nj*k) + size_tot*3;
    Ary1D ExkL(&data->sweep->vL[offset2],ni);
    Ary1D ExkR(&data->sweep->vR[offset2],ni);
    Ary1D dxkL(&data->sweep->vL[offset3],ni);
    Ary1D dxkR(&data->sweep->vR[offset3],ni);
    #pragma acc loop vector
    for (j = jbeg; j <= jend; j++){
      Ex1e(i,j,k) = 0.25*(ExkL(j) + ExkR(j)) + 0.5*(dxkL(j) + dxkR(j));
    }
    #endif  // INCLUDE_KDIR 

  }}
  #endif

  // ------------------------------------------------------
  // X3. Reconstruct along the x3-direction:
  //
  //     o Eyi(i-1/2) and Dx(i-1/2) --> (i-1/2, j, k-1/2)
  //     o Exj(j-1/2) and Dy(j-1/2) --> (i, j-1/2, k-1/2)
  // ------------------------------------------------------

#if INCLUDE_KDIR

  ARRAY_STRIDES (grid, KDIR, ni, nj);

  #if INCLUDE_IDIR
  ArrayReconstruct<KDIR> (data, emf->e[JDIR][IDIR],     IDIR, lim, grid, &box_rec, 0); 
  ArrayReconstruct<KDIR> (data, emf->e_dff[JDIR][IDIR], IDIR, lim, grid, &box_rec, 1); 
  #endif

  #if INCLUDE_JDIR
  ArrayReconstruct<KDIR> (data, emf->e[IDIR][JDIR],     JDIR, lim, grid, &box_rec, 2); 
  ArrayReconstruct<KDIR> (data, emf->e_dff[IDIR][JDIR], JDIR, lim, grid, &box_rec, 3); 
  #endif
  
  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(emf, grid)
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    #if INCLUDE_IDIR
    long int offset = ni*(i + nj*j);
    long int offset1 = ni*(i + nj*j) + size_tot;
    Ary1D EyiL(&data->sweep->vL[offset],ni);
    Ary1D EyiR(&data->sweep->vR[offset],ni);
    Ary1D dyiL(&data->sweep->vL[offset1],ni);
    Ary1D dyiR(&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (k = kbeg; k <= kend; k++){
     Ex2e(i,j,k) += 0.25*(EyiL(k) + EyiR(k)) + 0.5*(dyiL(k) + dyiR(k));
    }
    #endif // INCLUDE_IDIR 

    #if INCLUDE_JDIR
    long int offset2 = ni*(i + nj*j) + size_tot*2;
    long int offset3 = ni*(i + nj*j) + size_tot*3;
    Ary1D ExjL(&data->sweep->vL[offset2],ni);
    Ary1D ExjR(&data->sweep->vR[offset2],ni);
    Ary1D dxjL(&data->sweep->vL[offset3],ni);
    Ary1D dxjR(&data->sweep->vR[offset3],ni);

    #pragma acc loop vector
    for (k = kbeg; k <= kend; k++){
      Ex1e(i,j,k) += 0.25*(ExjL(k) + ExjR(k)) + 0.5*(dxjL(k) + dxjR(k));
    }
    #endif // INCLUDE_JDIR 

  }}
  #endif

#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Reconstruct the electric field to zone edge using the general
  UCT formalism,
  
  E = (aL*vL[i]*bL[i] + aR*vR[i]*bR[i]) +- (dR*bR - dL*bL)
  
  where the a and d coefficients are retrieved from the 1D Riemann solver
  at zone edges.
  This function is used by UCT_HLL, UCT_HLLD and UCT_GFORCE
*/
////////////////////////////////////////////////////////////////////////
void CT_EMF_Riemann2D(const Data *data, const EMF *emf, Grid *grid, RBox *box)
{
#if DIMENSIONS > 1

  int i, j, k;
  int ni, nj;
  int ibeg = box->ibeg, iend = box->iend+INCLUDE_IDIR;
  int jbeg = box->jbeg, jend = box->jend+INCLUDE_JDIR;
  int kbeg = box->kbeg, kend = box->kend+INCLUDE_KDIR;

  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];

  #ifdef HIGH_ORDER
  int ns = (HO_AVERAGE < 3 ? 1:2);
  const Ary3D *Vs = data->Vps;
  #else
  const Ary3D *Vs = data->Vs;
  #endif
  
  #if RECONSTRUCTION == LINEAR
  int recV = (LIMITER == DEFAULT ? VANLEER_LIM:LIMITER);
  int recB = MC_LIM;
  #else
  int recV = RECONSTRUCTION;
  int recB = RECONSTRUCTION;
  #endif
  
  #if BACKGROUND_FIELD == YES
  double *x1 = grid->xC[IDIR], *x1L = grid->xL[IDIR];
  double *x2 = grid->xC[JDIR], *x2L = grid->xL[JDIR];
  double *x3 = grid->xC[KDIR], *x3L = grid->xL[KDIR];
  RunConfig &run_config = *(data->run_config);
  #endif

  #if (INCLUDE_KDIR) && (INCLUDE_JDIR)
  auto &Ex1e = emf->Ee[IDIR];
  #endif
  #if (INCLUDE_KDIR) && (INCLUDE_IDIR)
  auto &Ex2e = emf->Ee[JDIR];
  #endif
  #if (INCLUDE_JDIR) && (INCLUDE_IDIR)
  auto &Ex3e = emf->Ee[KDIR];
  #endif

  RBox box_rec;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);

  mynvtxstart_("CT_EMF_Riemann2D_Allocation",YELLOW);

  // ------------------------------------------------------
  // X1a. Reconstruct along the x1-direction:
  //
  //      o By(j-1/2) and -vx(j-1/2) --> (i-1/2, j-1/2, k)
  //      o Initialize Ex3 at (i-1/2, j-1/2, k)
  // ------------------------------------------------------

#if INCLUDE_IDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                   iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg,                   jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg - ns*INCLUDE_KDIR, kend = box->kend + INCLUDE_KDIR*ns;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ARRAY_STRIDES (grid, IDIR, ni, nj);

  #if INCLUDE_JDIR
  ArrayReconstruct<IDIR> (data, Vs[BX2s],           JDIR, recB, grid, &box_rec, 0); 
  ArrayReconstruct<IDIR> (data, emf->e[KDIR][JDIR], JDIR, recV, grid, &box_rec, 1); 

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, emf, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset0 = ni*(j + nj*k);
    long int offset1 = offset0 + size_tot;

    Ary1D byL(&data->sweep->vL[offset0],ni);
    Ary1D byR(&data->sweep->vR[offset0],ni);
    Ary1D vL (&data->sweep->vL[offset1],ni);
    Ary1D vR (&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){

      #if BACKGROUND_FIELD == YES
      {
      double Bgf[3];
      State state;
      STATE_SET(state, NULL, NULL, Bgf, x1L[i], x2L[j], x3[k]);
      BackgroundField (state, run_config);
      byL(i) += state.Bgf[JDIR];
      byR(i) += state.Bgf[JDIR];
      }
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      double aL = 0.5;
      double aR = 0.5;
      #else
      double aL = 0.5*(emf->aL[IDIR](i,j-1,k) + emf->aL[IDIR](i,j,k));
      double aR = 0.5*(emf->aR[IDIR](i,j-1,k) + emf->aR[IDIR](i,j,k));
      #endif
      double dL = 0.5*(emf->dL[IDIR](i,j-1,k) + emf->dL[IDIR](i,j,k));
      double dR = 0.5*(emf->dR[IDIR](i,j-1,k) + emf->dR[IDIR](i,j,k));

      double phi = dR*byR(i) - dL*byL(i);
      Ex3e(i,j,k) = (aL*vL(i)*byL(i) + aR*vR(i)*byR(i)) + phi;

      #if CT_ENERGY_CORRECTION > 1
      emf->Bx2ez(i,j,k)    = aL*byL(i) + aR*byR(i);
      emf->Ex3e_phi(i,j,k) = phi;
emf->Ex3e_phix(i,j,k) = phi;
      #endif
    }
  }}
  #endif  // INCLUDE_JDIR 

  // ------------------------------------------------------
  // X1b. Reconstruct along the x1-direction:
  //
  //      o Bz(k-1/2) and vx(k-1/2) --> (i-1/2, j, k-1/2)
  //      o Initialize Ex2 at (i-1/2, j, k-1/2)
  // ------------------------------------------------------

  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                   iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg - ns*INCLUDE_JDIR, jend = box->jend + INCLUDE_JDIR*ns;
  kbeg = box->kbeg,                   kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  #if INCLUDE_KDIR
  ArrayReconstruct<IDIR> (data, Vs[BX3s],           KDIR, recB, grid, &box_rec, 0); 
  ArrayReconstruct<IDIR> (data, emf->e[JDIR][KDIR], KDIR, recV, grid, &box_rec, 1); 

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, emf, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset0 = ni*(j + nj*k);
    long int offset1 = offset0 + size_tot;

    Ary1D bzL(&data->sweep->vL[offset0],ni);
    Ary1D bzR(&data->sweep->vR[offset0],ni);
    Ary1D vL(&data->sweep->vL[offset1],ni);
    Ary1D vR(&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){

      #if BACKGROUND_FIELD == YES
      {
      double Bgf[3];
      State state;
      STATE_SET(state, NULL, NULL, Bgf, x1L[i], x2[j], x3L[k]);
      BackgroundField (state, run_config);
      bzL(i) += state.Bgf[KDIR];
      bzR(i) += state.Bgf[KDIR];
      }
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      double aL = 0.5;
      double aR = 0.5;
      #else
      double aL = 0.5*(emf->aL[IDIR](i,j,k-1) + emf->aL[IDIR](i,j,k));
      double aR = 0.5*(emf->aR[IDIR](i,j,k-1) + emf->aR[IDIR](i,j,k));
      #endif
      double dL = 0.5*(emf->dL[IDIR](i,j,k-1) + emf->dL[IDIR](i,j,k));
      double dR = 0.5*(emf->dR[IDIR](i,j,k-1) + emf->dR[IDIR](i,j,k));

      double phi = dR*bzR(i) - dL*bzL(i);
      Ex2e(i,j,k) = (aL*vL(i)*bzL(i) + aR*vR(i)*bzR(i)) - phi;

      #if CT_ENERGY_CORRECTION > 1
      emf->Bx3ey(i,j,k)    = aL*bzL(i) + aR*bzR(i);
      emf->Ex2e_phi(i,j,k) = -phi;
      #endif
    }
  }}
  #endif  // INCLUDE_KDIR 
#endif // endif INCLUDE_IDIR 

  // ------------------------------------------------------
  // X2a. Reconstruct along the x2-direction:
  //
  //      o Bz(k-1/2) and -vy(k-1/2)  --> (i, j-1/2, k-1/2)
  //      o Initialize Ex1 at (i,j-1/2, k-1/2)
  // ------------------------------------------------------

#if INCLUDE_JDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg - ns*INCLUDE_IDIR, iend = box->iend + INCLUDE_IDIR*ns;
  jbeg = box->jbeg,                   jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg,                   kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ARRAY_STRIDES (grid, JDIR, ni, nj);

  #if INCLUDE_KDIR
  ArrayReconstruct<JDIR> (data, Vs[BX3s],           KDIR, recB, grid, &box_rec, 0); 
  ArrayReconstruct<JDIR> (data, emf->e[IDIR][KDIR], KDIR, recV, grid, &box_rec, 1); 

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, emf, grid)
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
    long int offset0 = ni*(i + nj*k);
    long int offset1 = offset0 + size_tot;

    Ary1D bzL(&data->sweep->vL[offset0],ni);
    Ary1D bzR(&data->sweep->vR[offset0],ni);
    Ary1D vL (&data->sweep->vL[offset1],ni);
    Ary1D vR (&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (j = jbeg; j <= jend; j++){

      #if BACKGROUND_FIELD == YES
      {
      double Bgf[3];
      State state;
      STATE_SET(state, NULL, NULL, Bgf, x1[i], x2L[j], x3L[k]);
      BackgroundField (state, run_config);
      bzL(j) += state.Bgf[KDIR];
      bzR(j) += state.Bgf[KDIR];
      }
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      double aL = 0.5;
      double aR = 0.5;
      #else
      double aL = 0.5*(emf->aL[JDIR](i,j,k-1) + emf->aL[JDIR](i,j,k));
      double aR = 0.5*(emf->aR[JDIR](i,j,k-1) + emf->aR[JDIR](i,j,k));
      #endif
      double dL = 0.5*(emf->dL[JDIR](i,j,k-1) + emf->dL[JDIR](i,j,k));
      double dR = 0.5*(emf->dR[JDIR](i,j,k-1) + emf->dR[JDIR](i,j,k));

      double phi = dR*bzR(j) - dL*bzL(j);
      Ex1e(i,j,k) = (aL*vL(j)*bzL(j) + aR*vR(j)*bzR(j)) + phi;

      #if CT_ENERGY_CORRECTION > 1
      emf->Bx3ex(i,j,k) = aL*bzL(j) + aR*bzR(j);
      emf->Ex1e_phi(i,j,k) = phi;
      #endif
    }
  }}
  #endif // INCLUDE_KDIR 

  // ------------------------------------------------------
  // X2b. Reconstruct along the x2-direction:
  // 
  //      o Bx(i-1/2) and vy(i-1/2) --> (i-1/2, j-1/2, k)
  //      o Complete Ex3 at (i-1/2,j-1/2, k)
  // ------------------------------------------------------

  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                   iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg,                   jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg - ns*INCLUDE_KDIR, kend = box->kend + INCLUDE_KDIR*ns;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  #if INCLUDE_IDIR
  ArrayReconstruct<JDIR> (data, Vs[BX1s],           IDIR, recB, grid, &box_rec, 0); 
  ArrayReconstruct<JDIR> (data, emf->e[KDIR][IDIR], IDIR, recV, grid, &box_rec, 1); 

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, emf, grid)
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
    long int offset0 = ni*(i + nj*k);
    long int offset1 = offset0 + size_tot;

    Ary1D bxL(&data->sweep->vL[offset0],ni);
    Ary1D bxR(&data->sweep->vR[offset0],ni);
    Ary1D vL(&data->sweep->vL[offset1],ni);
    Ary1D vR(&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (j = jbeg; j <= jend; j++){

      #if BACKGROUND_FIELD == YES
      {
      double Bgf[3];
      State state;
      STATE_SET(state, NULL, NULL, Bgf, x1L[i], x2L[j], x3[k]);
      BackgroundField (state, run_config);
      bxL(j) += state.Bgf[IDIR];
      bxR(j) += state.Bgf[IDIR];
      }
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      double aL = 0.5;
      double aR = 0.5;
      #else
      double aL = 0.5*(emf->aL[JDIR](i-1,j,k) + emf->aL[JDIR](i,j,k));
      double aR = 0.5*(emf->aR[JDIR](i-1,j,k) + emf->aR[JDIR](i,j,k));
      #endif
      double dL = 0.5*(emf->dL[JDIR](i-1,j,k) + emf->dL[JDIR](i,j,k));
      double dR = 0.5*(emf->dR[JDIR](i-1,j,k) + emf->dR[JDIR](i,j,k));

      double phi = dR*bxR(j) - dL*bxL(j);
      Ex3e(i,j,k) += (aL*vL(j)*bxL(j) + aR*vR(j)*bxR(j)) - phi;

      #if CT_ENERGY_CORRECTION > 1
      emf->Bx1ez(i,j,k) = aL*bxL(j) + aR*bxR(j);
      emf->Ex3e_phi(i,j,k) -= phi;
emf->Ex3e_phiy(i,j,k) = phi;
      #endif
    }
  }}
  #endif // INCLUDE_IDIR 
#endif // endif INCLUDE_JDIR 

  // ------------------------------------------------------
  // X3a. Reconstruct along the x3-direction:
  //
  //      o Bx(i-1/2) and -vz(i-1/2)  -->  (i-1/2, j, k-1/2)
  //      o Complete Ex2 at (i-1/2,j,k-1/2)
  // ------------------------------------------------------

#if INCLUDE_KDIR

  #ifdef HIGH_ORDER
  ibeg = box->ibeg,                   iend = box->iend + INCLUDE_IDIR;
  jbeg = box->jbeg - ns*INCLUDE_JDIR, jend = box->jend + INCLUDE_JDIR*ns;
  kbeg = box->kbeg,                   kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  ARRAY_STRIDES (grid, KDIR, ni, nj);

  #if INCLUDE_IDIR
  ArrayReconstruct<KDIR> (data, Vs[BX1s],           IDIR, recB, grid, &box_rec, 0); 
  ArrayReconstruct<KDIR> (data, emf->e[JDIR][IDIR], IDIR, recV, grid, &box_rec, 1); 

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, emf, grid) 
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    long int offset0 = ni*(i + nj*j);
    long int offset1 = offset0 + size_tot;

    Ary1D bxL(&data->sweep->vL[offset0],ni);
    Ary1D bxR(&data->sweep->vR[offset0],ni);
    Ary1D vL (&data->sweep->vL[offset1],ni);
    Ary1D vR (&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (k = kbeg; k <= kend; k++){

      #if BACKGROUND_FIELD == YES
      {
      double Bgf[3];
      State state;
      STATE_SET(state, NULL, NULL, Bgf, x1L[i], x2[j], x3L[k]);
      BackgroundField (state, run_config);
      bxL(k) += state.Bgf[IDIR];
      bxR(k) += state.Bgf[IDIR];
      }
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      double aL = 0.5;
      double aR = 0.5;
      #else
      double aL = 0.5*(emf->aL[KDIR](i-1,j,k) + emf->aL[KDIR](i,j,k));
      double aR = 0.5*(emf->aR[KDIR](i-1,j,k) + emf->aR[KDIR](i,j,k));
      #endif
      double dL = 0.5*(emf->dL[KDIR](i-1,j,k) + emf->dL[KDIR](i,j,k));
      double dR = 0.5*(emf->dR[KDIR](i-1,j,k) + emf->dR[KDIR](i,j,k));

      double phi = dR*bxR(k) - dL*bxL(k);
      Ex2e(i,j,k) += (aL*vL(k)*bxL(k) + aR*vR(k)*bxR(k)) + phi;

      #if CT_ENERGY_CORRECTION > 1
      emf->Bx1ey(i,j,k) = aL*bxL(k) + aR*bxR(k);
      emf->Ex2e_phi(i,j,k) += phi;
      #endif
    }
  }}
  #endif // INCLUDE_IDIR 

  // ------------------------------------------------------
  // X3b. Reconstruct along the x3-direction:
  //
  //      o By(j-1/2) and vz(j-1/2)  --> (i, j-1/2, k-1/2)
  //      o Complete Ex1 at (i, j-1/2, k-1/2)
  // ------------------------------------------------------

  #ifdef HIGH_ORDER
  ibeg = box->ibeg - ns*INCLUDE_IDIR, iend = box->iend + INCLUDE_IDIR*ns;
  jbeg = box->jbeg,                   jend = box->jend + INCLUDE_JDIR;
  kbeg = box->kbeg,                   kend = box->kend + INCLUDE_KDIR;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  #endif

  #if INCLUDE_JDIR
  ArrayReconstruct<KDIR> (data, Vs[BX2s],           JDIR, recB, grid, &box_rec, 0); 
  ArrayReconstruct<KDIR> (data, emf->e[IDIR][JDIR], JDIR, recV, grid, &box_rec, 1); 

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, emf, grid) 
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    long int offset0 = ni*(i + nj*j);
    long int offset1 = offset0 + size_tot;

    Ary1D byL(&data->sweep->vL[offset0],ni);
    Ary1D byR(&data->sweep->vR[offset0],ni);
    Ary1D vL (&data->sweep->vL[offset1],ni);
    Ary1D vR (&data->sweep->vR[offset1],ni);

    #pragma acc loop vector
    for (k = kbeg; k <= kend; k++){

      #if BACKGROUND_FIELD == YES
      {
      double Bgf[3];
      State state;
      STATE_SET(state, NULL, NULL, Bgf, x1[i], x2L[j], x3L[k]);
      BackgroundField (state, run_config);
      byL(k) += state.Bgf[JDIR];
      byR(k) += state.Bgf[JDIR];
      }
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      double aL = 0.5;
      double aR = 0.5;
      #else
      double aL = 0.5*(emf->aL[KDIR](i,j-1,k) + emf->aL[KDIR](i,j,k));
      double aR = 0.5*(emf->aR[KDIR](i,j-1,k) + emf->aR[KDIR](i,j,k));
      #endif
      double dL = 0.5*(emf->dL[KDIR](i,j-1,k) + emf->dL[KDIR](i,j,k));
      double dR = 0.5*(emf->dR[KDIR](i,j-1,k) + emf->dR[KDIR](i,j,k));

      double phi = dR*byR(k) - dL*byL(k);
      Ex1e(i,j,k) += (aL*vL(k)*byL(k) + aR*vR(k)*byR(k)) - phi;

      #if CT_ENERGY_CORRECTION > 1
      emf->Bx2ex(i,j,k) = aL*byL(k) + aR*byR(k);
      emf->Ex1e_phi(i,j,k) -= phi;
      #endif
    }
  }}
  #endif // INCLUDE_JDIR 
#endif // endif INCLUDE_KDIR 

  mynvtxstop_();

#endif
}
