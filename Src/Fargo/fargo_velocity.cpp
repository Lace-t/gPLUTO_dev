////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Functions for computing/retrieving the mean aziumthal velocity.

  \authors A. Mignone (andrea.mignone@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it) \n
           G. Muscianisi (g.muscianisi@cineca.it)
  \date    Dec 06, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define USE_MPI_COMMUNICATOR  YES

/*! Defines a 2D array containing the azimuthally-averaged velocity.
    When the average is done along the X2 direction (Cartesian and polar
    geometries) the array should have size (NX3_TOT, NX1_TOT).
    When the average is performed along the X3 direction (spherical geometry)
    the array has dimensions (NX2_TOT, NX1_TOT)                           */
static double **wA;  
/*static double **wAx1, **wAx2, **wAx3;*/  // -- average orbital speed at
                                            // x1, x2 and x3 faces -- 
                                            
////////////////////////////////////////////////////////////////////////
/*!
  Initialize FARGO module.
*/
////////////////////////////////////////////////////////////////////////
void FARGO_Initialize(Grid *grid)
{
  int nx1_tot=grid->ntot[IDIR];
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  int nx3_tot=grid->ntot[KDIR];
  if (wA == NULL) {
    wA = ARRAY_2D(nx3_tot, nx1_tot, double);
    #pragma acc enter data copyin(wA[nx3_tot][:nx1_tot])
  }
  #else
  int nx2_tot=grid->ntot[JDIR];
  if (wA == NULL) {
    wA = ARRAY_2D(nx2_tot, nx1_tot, double);
    #pragma acc enter data copyin(wA[:nx2_tot][:nx1_tot])
  }
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Compute the background orbital velocity as a 2D array by averaging
  the total azimuthal velocity along the orbital coordinate y/phi/phi.
  Update the residual velocity as well.
*/
////////////////////////////////////////////////////////////////////////
void FARGO_AverageVelocity(Data *data, Grid *grid)
{
#if FARGO_NSTEP_AVERAGE > 0
mynvtxstart_("FARGO_AverageVelocity",BLUE);

  int i,j,k;
  double w;
  int nx1_tot = grid->ntot[IDIR];
  RBox box;

  RBoxDefine (grid->lbeg, grid->lend, &box);

  int condition = (g_stepNumber%FARGO_NSTEP_AVERAGE == 0) && g_stepNumber != 0;
  if (!condition) return;

  // ----------------------------------------------------
  // 0. Create new communicators for each row of procs
  //    lying in the same orbital direction.
  // ----------------------------------------------------

  #ifdef PARALLEL
  int    *nproc      = grid->nproc;
  int    *rank_coord = grid->rank_coord;
  double w_sum = 0.0;
  #if USE_MPI_COMMUNICATOR == NO
  double w_recv = 0.0;
  int count, coords[3], src, dst;
  MPI_Comm cartcomm;
  MPI_Status status;
  #endif
  MPI_Comm wComm;

  #if USE_MPI_COMMUNICATOR == YES
  #if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
  int color = rank_coord[KDIR]*nproc[IDIR] + rank_coord[IDIR];
  #elif GEOMETRY == SPHERICAL
  int color = rank_coord[JDIR]*nproc[IDIR] + rank_coord[IDIR];
  #endif
  MPI_Comm_split(MPI_COMM_WORLD, color, g_rank, &wComm);
  #endif
  #endif  // PARALLEL

  // ------------------------------------------------------
  // 1. Fill ghost zones
  // ------------------------------------------------------

  Boundary(data, 0, grid); 

  // ------------------------------------------------------
  // 2. get ranks of the upper and lower procs
  // ------------------------------------------------------

  // why not using #if USE_MPI_COMMUNICATOR == NO?
  #if !USE_MPI_COMMUNICATOR == YES
    #ifdef PARALLEL
  #if USE_AL == NO
  cartcomm = new_SZ.cart_comm;
  #else
    AL_Get_cart_comm(SZ, &cartcomm); 
  #endif
    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[SDIR] += 1;
    MPI_Cart_rank(cartcomm, coords, &dst);

    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[SDIR] -= 1;
    MPI_Cart_rank(cartcomm, coords, &src);
    #endif
  #endif

  // ------------------------------------------------------
  //  3. Compute average orbital speed in CARTESIAN 
  //     or CYLINDRICAL geometries
  // ------------------------------------------------------

  #if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
  int nx3_tot = grid->ntot[KDIR];
  int nx2     = grid->nint[JDIR];
  int jbeg    = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  #ifndef PARALLEL
  #pragma acc parallel loop collapse(2) present(data, wA)
  #endif
  for(k=0; k<nx3_tot; k++){
  for(i=0; i<nx1_tot; i++){
    
    // 3a. Average orbital speed

    w = 0.0;
    #ifndef PARALLEL
    #pragma acc loop vector reduction(+:w)
    #endif
    for(j=jbeg; j<=jend; j++){
      #if GEOMETRY == CARTESIAN
      double vphi = data->Vc(i,j,k,VX2);
      #else
      double vphi = data->Vc(i,j,k,iVPHI);
      #endif
      w += wA[k][i] + vphi;  //  Sum total velocity  
    }

    #if USE_MPI_COMMUNICATOR == YES
    if (grid->nproc[JDIR] > 1){
      #ifdef PARALLEL
      MPI_Allreduce (&w, &w_sum, 1, MPI_DOUBLE, MPI_SUM, wComm);
      w = w_sum/(double)(grid->nint_glob[JDIR]);
      #endif
    }else{
      w = w/(double)nx2;
    }
    #endif

    #if USE_MPI_COMMUNICATOR == NO
    if (grid->nproc[JDIR] > 1){
      #ifdef PARALLEL
      w_sum = w;
      for (count=1; count < grid->nproc[JDIR]; count++ ){
        MPI_Sendrecv(&w,      1, MPI_DOUBLE, dst, 0,
                      &w_recv, 1, MPI_DOUBLE, src, 0, cartcomm, &status);
        w      = w_recv;
        w_sum += w;
      }
      w = w_sum/(double)(grid->nint_glob[JDIR]);
      #endif
    }else{
      w = w/(double)nx2;
    }
    #endif

    // 3b. Update residual

    #ifndef PARALLEL
    #pragma acc loop vector
    #endif
    for(j=jbeg; j<=jend; j++){
      #if GEOMETRY == CARTESIAN
      data->Vc(i,j,k,VX2)   += wA[k][i] - w;
      #else
      data->Vc(i,j,k,iVPHI) += wA[k][i] - w;
      #endif
    }
    wA[k][i] = w;
  }}
  #endif

  // ------------------------------------------------------
  // 4. Compute average orbital speed in SPHERICAL 
  //    geometry.
  // ------------------------------------------------------

  #if GEOMETRY == SPHERICAL
  int nx2_tot = grid->ntot[JDIR];
  int nx3     = grid->nint[KDIR];
  int kbeg    = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  #ifndef PARALLEL
  #pragma acc parallel loop collapse(2) present(data, wA)
  #endif
  for(j=0; j<nx2_tot; j++){
  for(i=0; i<nx1_tot; i++){

    // 4a. Compute average orbital speed
      
    w = 0.0;
    #ifndef PARALLEL
    #pragma acc loop vector reduction(+:w)
    #endif
    for(k=kbeg; k<=kend; k++){
      w += wA[j][i] + data->Vc(i,j,k,iVPHI);
    }

    #if USE_MPI_COMMUNICATOR == YES
    if (grid->nproc[KDIR] > 1){
      #ifdef PARALLEL
      MPI_Allreduce (&w, &w_sum, 1, MPI_DOUBLE, MPI_SUM, wComm);
      w = w_sum/(double)(grid->nint_glob[KDIR]);
      #endif
    }else{
      w = w/(double)nx3;
    }
    #endif

    #if USE_MPI_COMMUNICATOR == NO
    if (grid->nproc[KDIR] > 1){
      #ifdef PARALLEL
      w_sum = w;
      for (count=1; count < grid->nproc[KDIR]; count++ ){
        MPI_Sendrecv(&w, 1, MPI_DOUBLE, dst, 0,  &w_recv, 1,
                            MPI_DOUBLE, src, 0, cartcomm, &status);
        w      = w_recv;
        w_sum += w;
      }
      w = w_sum/(double)(grid->nint_glob[KDIR]);
      #endif
    }else{
        w = w/(double)nx3;
    }
    #endif

    // 4b. Update residual
  
    #ifndef PARALLEL
    #pragma acc loop vector
    #endif
    for(k=kbeg; k<=kend; k++){
      #if GEOMETRY == CARTESIAN
      data->Vc(i,j,k,VX2) += wA[j][i] - w;
      #else             
      data->Vc(i,j,k,iVPHI) += wA[j][i] - w;
      #endif
    }

    wA[j][i] = w;
  }}
  #endif

  // ------------------------------------------------------
  //  5. Free communicator & convert to conservative
  //  -----------------------------------------------------

  #ifdef PARALLEL
  #if USE_MPI_COMMUNICATOR == YES
  MPI_Comm_free(&wComm);
  #endif
  #pragma acc update device(wA)
  #endif
  
mynvtxstop_();

  PrimToCons3D(data, 0, grid, &box);

#endif // FARGO_NSTEP_AVERAGE > 0 

}

////////////////////////////////////////////////////////////////////////
/*!
  Add the mean background contribution to the residual
  velocity in order to obtain the total velocity.
  
  \param [in]  d      pointer to PLUTO Data structure. It is assumed
  that data->Vc[VX2] contains the residual.
  \param [out] vtot   array where the total velocity is stored.
  (can be also data->Vc[VX2])
  \param [in]  grid   pointer to array of Grid structures
  
  \return This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
void FARGO_ComputeTotalVelocity(const Data *data, int nv, Grid *grid)
{
  int i,j,k;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  // ------------------------------------------------------
  // 1. Add background shear velocity to the fluid
  // ------------------------------------------------------

  for(k = 0; k < nx3_tot; k++){
  for(j = 0; j < nx2_tot; j++){
  for(i = 0; i < nx1_tot; i++){
    #if GEOMETRY == CARTESIAN 
    data->Vc(i,j,k,nv) = wA[k][i] + data->Vc(i,j,k,VX2);
    #elif GEOMETRY == CYLINDRICAL
    data->Vc(i,j,k,nv) = wA[k][i] + data->Vc(i,j,k,iVPHI);
    #elif GEOMETRY == SPHERICAL
    data->Vc(i,j,k,nv) = wA[j][i] + data->Vc(i,j,k,iVPHI);
    #else
    QuitPLUTO("! FARGO not supported in this geometry",1);
    #endif 
  }}}

  // ------------------------------------------------------
  // 2. Add background shear velocity to particles.
  //    Note: we do not call Particles_Interpolate since
  //    the array to be interpolated has only 2 dimensions
  //    and not three.
  // ------------------------------------------------------

#if PARTICLES
  int i1,j1,k1;
  particleNode *curNode;
  Particle *p;
  double v;
  static double ***wp;

  if (wp == NULL)  wp = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);

  PARTICLES_LOOP(curNode, d->PHead){
    p = &(curNode->p);

  // -- Interpolate fargo shift velocity at particle position -- 

    Particles_GetWeights(p, p->cell, wp, grid);  
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    v = 0.0;
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      v += wp[k1][j1][i1]*wA[k + k1][i + i1];
    }}}
//  v = -SB_Q*SB_OMEGA*p->coord[IDIR]; 
    p->speed[JDIR] += v;
    #elif GEOMETRY == SPHERICAL
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      v += wp[k1][j1][i1]*wA[j + j1][i + i1];
    }}}
    p->speed[KDIR] += v;
    #endif

  }
#endif  // PARTICLES 
}

////////////////////////////////////////////////////////////////////////
void FARGO_ComputeResidualVelocity(const Data *data, double ***vres, Grid *grid)
/*!
  Add the mean background contribution to the residual
  velocity in order to obtain the total velocity.
  
  \param [in]  d      pointer to PLUTO Data structure. It is assumed
  that data->Vc[VX2] contains the total vel.
  \param [out] vres   array where the residual velocity is stored.
  (can be also data->Vc[VX2])
  \param [in]  grid   pointer to array of Grid structures
  
  \return This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
{
//   int i,j,k;
// #if GEOMETRY == CARTESIAN
//   double ***vtot = data->Vc[VX2];
// #else
//   double ***vtot = data->Vc[iVPHI];
// #endif

// /* --------------------------------------------------------
//    2a. Subtract background shear velocity to the fluid
//    -------------------------------------------------------- */

//   TOT_LOOP(k,j,i){
//     #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
//     vres[k][j][i] = vtot[k][j][i] - wA[k][i];
//     #elif GEOMETRY == SPHERICAL
//     vres[k][j][i] = vtot[k][j][i] - wA[j][i];
//     #endif 
//   }

// /* --------------------------------------------------------
//    2a. Add background shear velocity to particles.
//        Note: we do not call Particles_Interpolate since
//        the array to be interpolated has only 2 dimensions
//        and not three.
//    -------------------------------------------------------- */

// #if PARTICLES
//   int i1,j1,k1;
//   particleNode *curNode;
//   Particle *p;
//   double v;
//   static double ***wp;

//   if (wp == NULL)  wp = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);

//   PARTICLES_LOOP(curNode, d->PHead){
//     p = &(curNode->p);

//   /* -- Interpolate fargo shift velocity at particle position -- */

//     Particles_GetWeights(p, p->cell, wp, grid);  
//     i = p->cell[IDIR];
//     j = p->cell[JDIR];
//     k = p->cell[KDIR];

//     v = 0.0;
//     #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
//     for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
//     for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
//     for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
//       v -= wp[k1][j1][i1]*wA[k + k1][i + i1];
//     }}}
// /*  v = -SB_Q*SB_OMEGA*p->coord[IDIR]; */
//     p->speed[JDIR] += v;
//     #elif GEOMETRY == SPHERICAL
//     for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
//     for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
//     for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
//       v -= wp[k1][j1][i1]*wA[j + j1][i + i1];
//     }}}
//     p->speed[KDIR] -= v;
//     #endif
//   }
// #endif  /* PARTICLES */
}


////////////////////////////////////////////////////////////////////////
double **FARGO_Velocity(void)
/*!
  Return a pointer to the background orbital velocity ::wA.
*/
////////////////////////////////////////////////////////////////////////
{
  return wA;
}
