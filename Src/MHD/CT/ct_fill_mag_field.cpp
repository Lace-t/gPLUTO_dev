////////////////////////////////////////////////////////////////////////
/*! \file
  \brief Assign the normal component of the staggered fields
  in the ghost zone-faces.
  
  Using div.B = 0 in staggered MHD or div.E = q (in ResRMHD), this
  function computes the staggered component of magnetic/electric field 
  lying on the zone-faces parallel to the boundary specified by "side".
  This is preformed by solving the div.B = 0 condition for one
  variable only which in 2-D requires the knowledge of the other
  3 components while in 3-D required the knowledge of the other 5
  staggered components.
  
  Note that this operation is performed in the outermost ghost zones
  only since the face values at grid->lbeg[IDIR] or grid->lend[IDIR]+1 
  are supposed to be part of the solution and are not changed 
  by this function.
  
  \verbatim
  +-----+-----+-----+-----+-----+--
  |     |     |     |     |     |
  X     X     X     |     |     |
  |     |     |     |     |     |
  +-----+-----+-----+-----+-----+--
  |
  <-----------------> BEG
  Physical boundary
  
  X = components assigned in this function.
  \endverbatim
  
  \author A. Mignone (andrea.mignone@unito.it)
  \date   Feb 2, 2026
  
   ///////////////////////////////////////////////////////////////////// */
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void CT_FillStaggeredFields (const Ary3D (&Vs)[6], const Ary3D &q,
                             int side, Grid *grid)
/*!
  
  \param [in]    Vs   pointer to staggered Ary3D
  \param [in]     q   pointer to charge Ary
  \param [in]  side   side of computational domain where ghost zones are filled
  \param [in]  grid   pointer to grid structure
*/
////////////////////////////////////////////////////////////////////////
{ 
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  int  i,j,k;
  DIM_EXPAND(double Bxp;, double Byp;, double Bzp;)
  DIM_EXPAND(double Bxm;, double Bym;, double Bzm;)
  double dBx = 0.0, dBy = 0.0, dBz = 0.0;
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  double qdV;
  double Exp, Eyp, Ezp;
  double Exm, Eym, Ezm;
  double dEx = 0.0, dEy = 0.0, dEz = 0.0;
  #endif

  DIM_EXPAND(double Axp;, double Ayp;, double Azp;)
  DIM_EXPAND(double Axm;, double Aym;, double Azm;)

  // ------------------------------------------------------
  // 0. Set initial and final index for each side
  // ------------------------------------------------------

  ibeg = 0; iend = grid->ntot[IDIR]-1;
  jbeg = 0; jend = grid->ntot[JDIR]-1;
  kbeg = 0; kend = grid->ntot[KDIR]-1;

  // ------------------------------------------------------
  // Loop over cells in the boundary
  // Divergence is written (for magnetic field) as
  //
  //       (Axp*Bxp - Axm*Bxm)
  //    + (Ayp*Byp - Aym*Bym)
  //    + (Azp*Bzp - Azm*Bzm) = 0
  //
  //    Divergence is written (for electric field) as
  //
  //      (Axp*Exp - Axm*Exm)
  //    + (Ayp*Eyp - Aym*Eym)
  //    + (Azp*Ezp - Azm*Ezm) = q*dV
  //
  //    and it is used to recover one components knowing 
  //    the other components and the charge.
  // ------------------------------------------------------

  // ------------------------------------------------------
  // 1. X1_BEG & X1_END Boundaries
  // ------------------------------------------------------
  
  #if INCLUDE_IDIR
  if (side == X1_BEG) {
    iend = grid->lbeg[IDIR]-1;

    #pragma omp target teams distribute parallel for collapse(2)
    #pragma acc parallel loop collapse(2) present(grid, Vs)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    #pragma acc loop seq
    for (i = iend; i >= ibeg; i--){

      DIM_EXPAND(Bxp = Vs[BX1s](i+1,j,k); Bxm = Vs[BX1s](i,j,k);  ,
                 Byp = Vs[BX2s](i,j+1,k); Bym = Vs[BX2s](i,j,k);  ,
                 Bzp = Vs[BX3s](i,j,k+1); Bzm = Vs[BX3s](i,j,k);  )

      DIM_EXPAND(Axp = grid->A[IDIR](i+1,j,k); Axm = grid->A[IDIR](i,j,k);  ,
                 Ayp = grid->A[JDIR](i,j+1,k); Aym = grid->A[JDIR](i,j,k);  ,
                 Azp = grid->A[KDIR](i,j,k+1); Azm = grid->A[KDIR](i,j,k);)

      DIM_EXPAND(dBx = (Axp*Bxp - Axm*Bxm);  ,
                 dBy = (Ayp*Byp - Aym*Bym);  ,
                 dBz = (Azp*Bzp - Azm*Bzm); )

      Vs[BX1s](i,j,k) = (Axp*Bxp + dBy + dBz)/Axm;

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      DIM_EXPAND(Exp = Vs[EX1s](i+1,j,k); Exm = Vs[EX1s](i,j,k);  ,
                 Eyp = Vs[EX2s](i,j+1,k); Eym = Vs[EX2s](i,j,k);  ,
                 Ezp = Vs[EX3s](i,j,k+1); Ezm = Vs[EX3s](i,j,k);  )

      DIM_EXPAND(dEx = (Axp*Exp - Axm*Exm);  ,
                 dEy = (Ayp*Eyp - Aym*Eym);  ,
                 dEz = (Azp*Ezp - Azm*Ezm); )

      qdV = q(i,j,k)*grid->dV(i,j,k);
      Vs[EX1s](i,j,k) = (Axp*Exp + dEy + dEz - qdV)/Axm;
      #endif

    }}}
  }

  if (side == X1_END) {
    ibeg = grid->lend[IDIR]+1;

    #pragma omp target teams distribute parallel for collapse(2)
    #pragma acc parallel loop collapse(2) present(grid, Vs)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    #pragma acc loop seq
    for (i = ibeg; i <= iend; i++){

      DIM_EXPAND(Bxp = Vs[BX1s](i+1,j,k); Bxm = Vs[BX1s](i,j,k);  ,
                 Byp = Vs[BX2s](i,j+1,k); Bym = Vs[BX2s](i,j,k);  ,
                 Bzp = Vs[BX3s](i,j,k+1); Bzm = Vs[BX3s](i,j,k);  )

      DIM_EXPAND(Axp = grid->A[IDIR](i+1,j,k); Axm = grid->A[IDIR](i,j,k);  ,
                 Ayp = grid->A[JDIR](i,j+1,k); Aym = grid->A[JDIR](i,j,k);  ,
                 Azp = grid->A[KDIR](i,j,k+1); Azm = grid->A[KDIR](i,j,k);)

      DIM_EXPAND(dBx = (Axp*Bxp - Axm*Bxm);  ,
                 dBy = (Ayp*Byp - Aym*Bym);  ,
                 dBz = (Azp*Bzp - Azm*Bzm); )

      Vs[BX1s](i+1,j,k) = (Axm*Bxm - (dBy + dBz))/Axp;

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      DIM_EXPAND(Exp = Vs[EX1s](i+1,j,k); Exm = Vs[EX1s](i,j,k);  ,
                 Eyp = Vs[EX2s](i,j+1,k); Eym = Vs[EX2s](i,j,k);  ,
                 Ezp = Vs[EX3s](i,j,k+1); Ezm = Vs[EX3s](i,j,k);  )

      DIM_EXPAND(dEx = (Axp*Exp - Axm*Exm);  ,
                 dEy = (Ayp*Eyp - Aym*Eym);  ,
                 dEz = (Azp*Ezp - Azm*Ezm); )

      qdV = q(i,j,k)*grid->dV(i,j,k);
      Vs[EX1s](i+1,j,k) = (Axm*Exm - (dEy + dEz) + qdV)/Axp;
      #endif

    }}}
  }
  #endif // INCLUDE_IDIR 

  // ------------------------------------------------------
  // 2. X2_BEG & X2_END Boundaries
  // ------------------------------------------------------

  #if INCLUDE_JDIR
  if (side == X2_BEG) {
    jend = grid->lbeg[JDIR]-1;

    #pragma omp target teams distribute parallel for collapse(2)
    #pragma acc parallel loop collapse(2) present(grid,Vs)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    #pragma acc loop seq
    for (j = jend; j >= jbeg; j--){

      DIM_EXPAND(Bxp = Vs[BX1s](i+1,j,k); Bxm = Vs[BX1s](i,j,k);  ,
                 Byp = Vs[BX2s](i,j+1,k); Bym = Vs[BX2s](i,j,k);  ,
                 Bzp = Vs[BX3s](i,j,k+1); Bzm = Vs[BX3s](i,j,k);  )

      DIM_EXPAND(Axp = grid->A[IDIR](i+1,j,k); Axm = grid->A[IDIR](i,j,k);  ,
                 Ayp = grid->A[JDIR](i,j+1,k); Aym = grid->A[JDIR](i,j,k);  ,
                 Azp = grid->A[KDIR](i,j,k+1); Azm = grid->A[KDIR](i,j,k);)

      DIM_EXPAND(dBx = (Axp*Bxp - Axm*Bxm);  ,
                 dBy = (Ayp*Byp - Aym*Bym);  ,
                 dBz = (Azp*Bzp - Azm*Bzm); )

      Vs[BX2s](i,j,k) = (Ayp*Byp + dBx + dBz)/Aym;

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      DIM_EXPAND(Exp = Vs[EX1s](i+1,j,k); Exm = Vs[EX1s](i,j,k);  ,
                 Eyp = Vs[EX2s](i,j+1,k); Eym = Vs[EX2s](i,j,k);  ,
                 Ezp = Vs[EX3s](i,j,k+1); Ezm = Vs[EX3s](i,j,k);  )

      DIM_EXPAND(dEx = (Axp*Exp - Axm*Exm);  ,
                 dEy = (Ayp*Eyp - Aym*Eym);  ,
                 dEz = (Azp*Ezp - Azm*Ezm); )

      qdV = q(i,j,k)*grid->dV(i,j,k);
      Vs[EX2s](i,j,k) = (Ayp*Eyp + dEx + dEz - qdV)/Aym;
      #endif

    }}}
  }

  if (side == X2_END) {
    jbeg = grid->lend[JDIR]+1;
    #pragma omp target teams distribute parallel for collapse(2)
    #pragma acc parallel loop collapse(2) present(grid, Vs)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    #pragma acc loop seq
    for (j = jbeg; j <= jend; j++){

      DIM_EXPAND(Bxp = Vs[BX1s](i+1,j,k); Bxm = Vs[BX1s](i,j,k);  ,
                 Byp = Vs[BX2s](i,j+1,k); Bym = Vs[BX2s](i,j,k);  ,
                 Bzp = Vs[BX3s](i,j,k+1); Bzm = Vs[BX3s](i,j,k);  )

      DIM_EXPAND(Axp = grid->A[IDIR](i+1,j,k); Axm = grid->A[IDIR](i,j,k);  ,
                 Ayp = grid->A[JDIR](i,j+1,k); Aym = grid->A[JDIR](i,j,k);  ,
                 Azp = grid->A[KDIR](i,j,k+1); Azm = grid->A[KDIR](i,j,k);)

      DIM_EXPAND(dBx = (Axp*Bxp - Axm*Bxm);  ,
                 dBy = (Ayp*Byp - Aym*Bym);  ,
                 dBz = (Azp*Bzp - Azm*Bzm); )
      Vs[BX2s](i,j+1,k) = (Aym*Bym - (dBx + dBz))/Ayp;

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      DIM_EXPAND(Exp = Vs[EX1s](i+1,j,k); Exm = Vs[EX1s](i,j,k);  ,
                 Eyp = Vs[EX2s](i,j+1,k); Eym = Vs[EX2s](i,j,k);  ,
                 Ezp = Vs[EX3s](i,j,k+1); Ezm = Vs[EX3s](i,j,k);  )

      DIM_EXPAND(dEx = (Axp*Exp - Axm*Exm);  ,
                 dEy = (Ayp*Eyp - Aym*Eym);  ,
                 dEz = (Azp*Ezp - Azm*Ezm); )

      qdV = q(i,j,k)*grid->dV(i,j,k);
      Vs[EX2s](i,j+1,k) = (Aym*Eym - (dEx + dEz) + qdV )/Ayp;
      #endif

    }}}
  }
  #endif  // INCLUDE_JDIR 

  // ------------------------------------------------------
  // 3. X3_BEG & X3_END Boundaries
  // ------------------------------------------------------

  #if INCLUDE_KDIR
  if (side == X3_BEG) {
    kend = grid->lbeg[KDIR]-1;
    
    #pragma omp target teams distribute parallel for collapse(2)
    #pragma acc parallel loop collapse(2) present(grid, Vs)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    #pragma acc loop seq
    for (k = kend; k >= kbeg; k--){

      DIM_EXPAND(Bxp = Vs[BX1s](i+1,j,k); Bxm = Vs[BX1s](i,j,k);  ,
                 Byp = Vs[BX2s](i,j+1,k); Bym = Vs[BX2s](i,j,k);  ,
                 Bzp = Vs[BX3s](i,j,k+1); Bzm = Vs[BX3s](i,j,k);  )


      DIM_EXPAND(Axp = grid->A[IDIR](i+1,j,k); Axm = grid->A[IDIR](i,j,k);  ,
                 Ayp = grid->A[JDIR](i,j+1,k); Aym = grid->A[JDIR](i,j,k);  ,
                 Azp = grid->A[KDIR](i,j,k+1); Azm = grid->A[KDIR](i,j,k);)

      DIM_EXPAND(dBx = (Axp*Bxp - Axm*Bxm);  ,
                 dBy = (Ayp*Byp - Aym*Bym);  ,
                 dBz = (Azp*Bzp - Azm*Bzm); )

      Vs[BX3s](i,j,k)= (Azp*Bzp + dBx + dBy)/Azm;

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      DIM_EXPAND(Exp = Vs[EX1s](i+1,j,k); Exm = Vs[EX1s](i,j,k);  ,
                 Eyp = Vs[EX2s](i,j+1,k); Eym = Vs[EX2s](i,j,k);  ,
                 Ezp = Vs[EX3s](i,j,k+1); Ezm = Vs[EX3s](i,j,k);  )

      DIM_EXPAND(dEx = (Axp*Exp - Axm*Exm);  ,
                 dEy = (Ayp*Eyp - Aym*Eym);  ,
                 dEz = (Azp*Ezp - Azm*Ezm); )

      qdV = q(i,j,k)*grid->dV(i,j,k);
      Vs[EX3s](i,j,k) = (Azp*Ezp + dEx + dEy - qdV)/Azm;
      #endif

    }}}
  }

  if (side == X3_END) {
    kbeg = grid->lend[KDIR]+1;
    #pragma omp target teams distribute parallel for collapse(2)
    #pragma acc parallel loop collapse(2) present(grid, Vs)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    #pragma acc loop seq
    for (k = kbeg; k <= kend; k++){

      DIM_EXPAND(Bxp = Vs[BX1s](i+1,j,k); Bxm = Vs[BX1s](i,j,k);  ,
                 Byp = Vs[BX2s](i,j+1,k); Bym = Vs[BX2s](i,j,k);  ,
                 Bzp = Vs[BX3s](i,j,k+1); Bzm = Vs[BX3s](i,j,k);  )

      DIM_EXPAND(Axp = grid->A[IDIR](i+1,j,k); Axm = grid->A[IDIR](i,j,k);  ,
                 Ayp = grid->A[JDIR](i,j+1,k); Aym = grid->A[JDIR](i,j,k);  ,
                 Azp = grid->A[KDIR](i,j,k+1); Azm = grid->A[KDIR](i,j,k);)

      DIM_EXPAND(dBx = (Axp*Bxp - Axm*Bxm);  ,
                 dBy = (Ayp*Byp - Aym*Bym);  ,
                 dBz = (Azp*Bzp - Azm*Bzm); )

      Vs[BX3s](i,j,k+1) = (Azm*Bzm - (dBx + dBy))/Azp;

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      DIM_EXPAND(Exp = Vs[EX1s](i+1,j,k); Exm = Vs[EX1s](i,j,k);  ,
                 Eyp = Vs[EX2s](i,j+1,k); Eym = Vs[EX2s](i,j,k);  ,
                 Ezp = Vs[EX3s](i,j,k+1); Ezm = Vs[EX3s](i,j,k);  )

      DIM_EXPAND(dEx = (Axp*Exp - Axm*Exm);  ,
                 dEy = (Ayp*Eyp - Aym*Eym);  ,
                 dEz = (Azp*Ezp - Azm*Ezm); )

      qdV = q(i,j,k)*grid->dV(i,j,k);
      Vs[EX3s](i,j,k+1) = (Azm*Ezm - (dEx + dEy) + qdV)/Azp;
      #endif

    }}}
  }
  #endif // INCLUDE_KDIR 

  #if DIVE_CONTROL != CONSTRAINED_TRANSPORT
  #endif
}

