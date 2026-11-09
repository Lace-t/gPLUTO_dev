////////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Store or retrieve the Electromotive Force (EMF).

  This file provides a database functionality for storing or
  retrieving EMF components and related information at different
  points and times in the code.

  The CT_Allocate() allocates memory for EMF structure members.
  The CT_StoreEMF() function is called immediately after a 1D Riemann
  solver during the hydro sweeps in order to save Fluxes at cell
  interfaces and characteristic signal velocities into the emf
  structure for later reuse.
  The fluxes coming from different sweeps are the different components
  of the advective part (-v X B) part of the electric field.

  The function CT_GetEMF() is used to obtain the edge-centered electric
  field by properly averaging the EMF components previously stored
  at the zone faces during the 1D sweeps.

  \author  A. Mignone  (andrea.mignone@unito.it)
           V. Berta    (vittoria.berta@unito.it)
           M. Rossazza (marco.rossazza@unito.it)

  \date    Sep 26, 2026
*/
////////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

void CT_IntegrateEMF (EMF *emf, Grid *grid, RBox *box);

////////////////////////////////////////////////////////////////////////////
void CT_AllocateEMF (EMF *emf0, Grid *grid)
/*!
  Allocate memory for EMF structure and
  check for incompatible combinations of algorithms.
  
  CT AVERAGE  |E(sweep)|E(pnt)|E(dff)|sgn(vs)|S(L/R)|  vt  |a(L/R)|d(L/R)|
  ------------|--------|------|------|-------|------|------|------|------|
  ARITHMETIC  |  x x   |      |      |       |      |      |      |      |
  CT_CONTACT  |  x x   |      |      |   x   |      |      |      |      |
  CT_FLUX     |        | x  x | x  x |       |      |      |      |      |
  UCT_HLL     |        |      |      |       | x  x | x  x |      |      |
  UCT_HLLD    |        |      |      |       |      | x  x | x  x | x  x |
  UCT_GFORCE  |        |      |      |       |      | x  x |      | x  x |
  
  The "x" represent how many 3D arrays per interface are required.
 */
////////////////////////////////////////////////////////////////////////////
{
  int *ntot = grid->ntot;

#if PARTICLES
  #if (PARTICLES == PARTICLES_CR) && (PARTICLES_CR_FEEDBACK == YES)
    #if (CT_EMF_AVERAGE != ARITHMETIC) && (CT_EMF_AVERAGE != CT_CONTACT) && \
        (CT_EMF_AVERAGE != CT_FLUX)
        #error ! CT_Allocate(): EMF average not compatible with particles. Use ARITHMETIC / CT_CONTACT
    #endif
  #endif
#endif
 
  emf0->dL[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->dR[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->aL[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->aR[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->dL[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->dR[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->aL[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->aR[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->dL[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->dR[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->aL[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->aR[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  // emf0->dxL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->dxR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->axL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->axR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  // emf0->dyL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->dyR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->ayL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->ayR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  // emf0->dzL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->dzR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->azL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  // emf0->azR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  // ------------------------------------------------------
  // 1. Allocate memory for edge-centered ("e") electric
  //    and magnetic (only for ResRMHD) fields
  // ------------------------------------------------------

  DIM_EXPAND(                                                     ;  ,
           emf0->Ee[KDIR].create(ntot[IDIR]+1,ntot[JDIR]+1,ntot[KDIR]);  ,
           emf0->Ee[IDIR].create(ntot[IDIR],  ntot[JDIR]+1,ntot[KDIR]+1);
           emf0->Ee[JDIR].create(ntot[IDIR]+1,ntot[JDIR],  ntot[KDIR]+1);)

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  DIM_EXPAND(                                                     ;  ,
           emf0->Be[KDIR].create(ntot[IDIR]+1,ntot[JDIR]+1,ntot[KDIR]);  ,
           emf0->Be[IDIR].create(ntot[IDIR],  ntot[JDIR]+1,ntot[KDIR]+1);
           emf0->Be[JDIR].create(ntot[IDIR]+1,ntot[JDIR],  ntot[KDIR]+1);)
  #endif

  #if CT_EMF_AVERAGE == CT_CONTACT
  DIM_EXPAND(emf0->svx.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);  ,
             emf0->svy.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);  ,
             emf0->svz.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);)
  emf0->sv[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->sv[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->sv[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif

  // ------------------------------------------------------
  // 2. Allocate memory for face-centered electric and
  //    magnetic fields
  // ------------------------------------------------------

  #if    (CT_EMF_AVERAGE == UCT_HLL) \
      || (CT_EMF_AVERAGE == UCT_HLLD) || (CT_EMF_AVERAGE == UCT_GFORCE)
    DIM_EXPAND(emf0->SxL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
              emf0->SxR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);  ,

              emf0->SyL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
              emf0->SyR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);  ,

              emf0->SzL.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
              emf0->SzR.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);)
  #endif

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  emf0->Frho[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Frho[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Frho[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif

  emf0->e[KDIR][IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e[JDIR][IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e[KDIR][JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e[IDIR][JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e[IDIR][KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e[JDIR][KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  
  emf0->e_dff[KDIR][IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e_dff[JDIR][IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e_dff[KDIR][JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e_dff[IDIR][JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e_dff[IDIR][KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->e_dff[JDIR][KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  #if (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  emf0->B[KDIR][IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->B[JDIR][IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->B[KDIR][JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->B[IDIR][JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->B[IDIR][KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->B[JDIR][KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif

  // ------------------------------------------------------
  // 3. Additional arrays for conservative energy correction
  // ------------------------------------------------------

  #if CT_ENERGY_CORRECTION > 1
  emf0->Bx1ey.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Bx1ez.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->Bx2ex.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Bx2ez.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->Bx3ex.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Bx3ey.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->Ex1e_phi.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Ex2e_phi.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Ex3e_phi.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->Ex3e_phix.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Ex3e_phiy.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  emf0->Pfx.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Pfy.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  emf0->Pfz.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif
}

////////////////////////////////////////////////////////////////////////////
void CT_ComputeEMF (const Data *data, Grid *grid, RBox *box)
/*!
  Compute the total electromotive force (EMF) using the
  hyperbolic fluxes stored during the most recent upwind step.
  
  \param [in] d
  \param [in] grid
 */
////////////////////////////////////////////////////////////////////////////
{
  EMF  *emf = data->emf;

  // 0. Check if cell-centered emf is needed.
  //    This may depend on the scheme.

  #if    CT_EMF_AVERAGE == CT_CONTACT \
      || CT_EMF_AVERAGE == UCT0 \
      || ((PARTICLES == PARTICLES_CR) && (PHYSICS != ResRMHD))
  CT_ComputeCenterEMF (data, grid);
  #endif

  // 1. Select EMF-averaging scheme

  #if CT_EMF_AVERAGE == ARITHMETIC

  CT_EMF_ArithmeticAverage (emf, 0.25, box);

  #elif CT_EMF_AVERAGE == CT_CONTACT
  int i,j,k;
  int ibeg = box->ibeg, iend = box->iend+INCLUDE_IDIR;
  int jbeg = box->jbeg, jend = box->jend+INCLUDE_JDIR;
  int kbeg = box->kbeg, kend = box->kend+INCLUDE_KDIR;
  
  #if (PHYSICS == ResRMHD)
    #error CT_CONTACT not defined for ResRMHD
  #endif

  CT_EMF_ArithmeticAverage (emf, 1.0, box);
  CT_EMF_IntegrateToCorner (data, emf, box);

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(emf)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    #if INCLUDE_JDIR && INCLUDE_KDIR
    emf->Ee[IDIR](i,j,k) *= 0.25;
    #endif
    #if INCLUDE_IDIR && INCLUDE_KDIR
    emf->Ee[JDIR](i,j,k) *= 0.25;
    #endif
    #if INCLUDE_IDIR && INCLUDE_JDIR
    emf->Ee[KDIR](i,j,k) *= 0.25;
    #endif
  }}}

  #elif CT_EMF_AVERAGE == CT_FLUX

  CT_EMF_Flux(data, emf, grid, box);
 
  #elif CT_EMF_AVERAGE == CT_MAXWELL

  CT_MaxwellSolver (data, emf, grid, box);

  #elif    (CT_EMF_AVERAGE == UCT_HLLD) || (CT_EMF_AVERAGE == UCT_GFORCE) \
        || (CT_EMF_AVERAGE == UCT_HLL)

  #if (PHYSICS == ResRMHD)
    #error UCT AVERAGEs not defined for ResRMHD
  #endif

  CT_EMF_Riemann2D(data, emf, grid, box);

  #else
    #error Invalid CT_EMF_AVERAGE
  #endif

  // 2. Deconvolution of electric field (only for HO 3D)

  #if (defined HIGH_ORDER) && (DIMENSIONS == 3)
  CT_IntegrateEMF (emf, grid, box);
  #endif

}

#if (PHYSICS == MHD) && (RESISTIVITY != NO)
////////////////////////////////////////////////////////////////////////////
void CT_ResistiveEMF (const Data *data, int op, Grid *grid, RBox *box)
/*!
  Compute/Add resistive terms to EMF (used during parabolic update
  or STS operator splitting)
  
  \param [in] d
  \param [in] op     operation: op = 0 means initialize,
  op = 1 means add.
  \param [in] grid
  
  \return a pointer to an edge-centered EMF.
 */
////////////////////////////////////////////////////////////////////////////
{
  int    i, j, k;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  EMF *emf = data->emf;

  if (op == 0){
    BOX_STAG_LOOP(*box, i, j, k){
      #if INCLUDE_JDIR && INCLUDE_KDIR
      emf->Ee[IDIR](i,j,k) = data->eta[IDIR](i,j,k)*data->J[IDIR](i,j,k);
      #endif
      #if INCLUDE_IDIR && INCLUDE_KDIR
      emf->Ee[JDIR](i,j,k) = data->eta[JDIR](i,j,k)*data->J[JDIR](i,j,k);
      #endif
      #if INCLUDE_IDIR && INCLUDE_JDIR
      emf->Ee[KDIR](i,j,k) = data->eta[KDIR](i,j,k)*data->J[KDIR](i,j,k);
      #endif
    }
  }else if (op == 1){
    BOX_STAG_LOOP(*box, i, j, k){
      #if INCLUDE_JDIR && INCLUDE_KDIR
      emf->Ee[IDIR](i,j,k) += data->eta[IDIR](i,j,k)*data->J[IDIR](i,j,k);
      #endif
      #if INCLUDE_IDIR && INCLUDE_KDIR
      emf->Ee[JDIR](i,j,k) += data->eta[JDIR](i,j,k)*data->J[JDIR](i,j,k);
      #endif
      #if INCLUDE_IDIR && INCLUDE_JDIR
      emf->Ee[KDIR](i,j,k) += data->eta[KDIR](i,j,k)*data->J[KDIR](i,j,k);
      #endif
    }
  }else{
    QuitPLUTO("! CT_ResistiveEFM(): invalid op\n",1);
  }
}
#endif // RESISTIVITY != NO 

// #undef eps_CT_CONTACT

////////////////////////////////////////////////////////////////////////////
void CT_ComputeCenterEMF(const Data *data, const Grid *grid)
/*!
  Compute cell-center inductive electric field.
 */
////////////////////////////////////////////////////////////////////////////
{
  int i,j,k;
  double vx1, vx2, vx3;
  double Bx1, Bx2, Bx3;

  vx1 = vx2 = vx3 = 0.0;
  Bx1 = Bx2 = Bx3 = 0.0;
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = 0; k < grid->ntot[KDIR]; k++){
  for (j = 0; j < grid->ntot[JDIR]; j++){
  for (i = 0; i < grid->ntot[IDIR]; i++){

    vx1 = data->Vc(i,j,k,VX1);
    vx2 = data->Vc(i,j,k,VX2);
    vx3 = data->Vc(i,j,k,VX3);

    Bx1 = data->Vc(i,j,k,BX1);
    Bx2 = data->Vc(i,j,k,BX2);
    Bx3 = data->Vc(i,j,k,BX3);

  // -- Compute inductive electric field -- 

    data->Ex1(i,j,k) = (vx3*Bx2 - vx2*Bx3);
    data->Ex2(i,j,k) = (vx1*Bx3 - vx3*Bx1);
    data->Ex3(i,j,k) = (vx2*Bx1 - vx1*Bx2);

  }}}
}

#ifdef HIGH_ORDER
////////////////////////////////////////////////////////////////////////////
void CT_IntegrateEMF (EMF *emf, Grid *grid, RBox *box)
/*!
  Deconvolve the edge-centered electric field into 
  edge-averaged electric field
  
  \param [in] d
  \param [in] emf
  \param [in] grid
  \param [in] RBox
 */
////////////////////////////////////////////////////////////////////////////
{
  int i, j, k;
  #if HO_AVERAGE < 3
  int ns = 0;
  #else
  int ns = 1;
  #endif
  int ibeg = (box->ibeg - ns), iend = (box->iend + ns);
  int jbeg = (box->jbeg - ns), jend = (box->jend + ns);
  int kbeg = (box->kbeg - ns), kend = (box->kend + ns);
  CT_ASSIGN(static double *d2e;,
            static double *d2b;)
  int nmax = grid->nmax;

  // 0. Allocate memory for edge-centered auxiliary 
  //    electric fields

  if (d2e == NULL) {
    CT_ASSIGN(d2e = ARRAY_1D(nmax, double);,
              d2b = ARRAY_1D(nmax, double);)
  }

  // ------------------------------------------------------
  // 1. Integration of the edge-centered E/B field:
  //    Use edge-centered point value to obtain 
  //    line-integral.
  // ------------------------------------------------------

  // 1a. Integration of Ex1s, Bx1s along x 
  // #pragma omp target should be added but is not simple
  #pragma acc parallel loop collapse(2) private(d2e[:nmax])
  for (k = kbeg; k <= kend+1; k++){
  for (j = jbeg; j <= jend+1; j++){
    #pragma acc loop
    for (i = ibeg; i <= iend; i++) {
      d2e[i] = HO_Average3D<IDIR>(emf->Ee[IDIR],i,j,k);
    }
    #pragma acc loop
    for (i = ibeg; i <= iend; i++) {
      emf->Ee[IDIR](i,j,k) += d2e[i]/24.0;
    }

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #pragma acc loop
    for (i = ibeg; i <= iend; i++) {
      d2b[i] = HO_Average3D<IDIR>(emf->Bx1e,i,j,k);
    }
    #pragma acc loop
    for (i = ibeg; i <= iend; i++) {
      emf->Bx1e(i,j,k) += d2b[i]/24.0;
    }
    #endif // (PHYSICS == ResRMHD) 
  }} // End for (i,j,k)

  // 1b. Integration of Ex2s, Bx2s along y

  #pragma acc parallel loop collapse(2) private(d2e[:nmax])
  for (k = kbeg; k <= kend+1; k++){
  for (i = ibeg; i <= iend+1; i++){
    #pragma acc loop
    for (j = jbeg; j <= jend; j++) {
      d2e[j] = HO_Average3D<JDIR>(emf->Ee[JDIR],i,j,k);
    }
    #pragma acc loop
    for (j = jbeg; j <= jend; j++) {
      emf->Ee[JDIR](i,j,k) += d2e[j]/24.0;
    }

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #pragma acc loop
    for (j = jbeg; j <= jend; j++) {
      d2b[j] = HO_Average3D<JDIR>(emf->Bx2e,i,j,k);
    }
    #pragma acc loop
    for (j = jbeg; j <= jend; j++) {
      emf->Bx2e(i,j,k) += d2b[j]/24.0;
    }
    #endif // (PHYSICS == ResRMHD) 
  }} // End for (i,j,k)

  // 1c. Integration of Ex3s, Bx3s along z

  #pragma acc parallel loop collapse(2)   private(d2e[:nmax])
  for (i = ibeg; i <= iend+1; i++){
  for (j = jbeg; j <= jend+1; j++){
    #pragma acc loop
    for (k = kbeg; k <= kend; k++) {
      d2e[k] = HO_Average3D<KDIR>(emf->Ee[KDIR],i,j,k);
    }
    #pragma acc loop
    for (k = kbeg; k <= kend; k++) {
      emf->Ee[KDIR](i,j,k) += d2e[k]/24.0;
    }

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #pragma acc loop
    for (k = kbeg; k <= kend; k++) {
      d2b[k] = HO_Average3D<KDIR>(emf->Bx3e,i,j,k);
    }
    #pragma acc loop
    for (k = kbeg; k <= kend; k++) {
      emf->Bx3e(i,j,k) += d2b[k]/24.0;
    }
    #endif // (PHYSICS == ResRMHD) 
  }} // End for (i,j,k)

}
#endif
