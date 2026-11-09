////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Ring-average method to treat singularity at a polar axis

  Implement the ring-average technique by Zhang et al. (2019) JCP, to
  remove the severe time step restriction for multi-dimensional
  curvilinear finite-volume MHD solvers near a polar axis.

  The average is carried out on conservative variables and magnetic 
  fluxes inside a chunk:
  
       Original grid:            |---|---|---|---|---|---|---|---| 
                                   A1  A2  B1  B2  C1  C2  D1  D2
      
       Averaged cells (chunk):   |-------|-------|-------|-------| 
                                     A       B       C       D
  
       After average:            |---|---|---|---|---|---|---|---| 
                                   A   A   B   B   C   C   D   D

  The Ring Average technique is applied before each RK stage and 
  thus requires no changes to the existing data structure,
  grid definition or numerical methods in the original MHD solver.

  \b References:
    - "Conservative averaging-reconstruction techniques (Ring Average)
       for 3-D finite-volume MHD solvers with axis singularity",
       Zhang et al, JCP (2019) 376, 276-294
       
  \authors A. Mignone   (andrea.mignone@unito.it)\n
           M. Rossazza 
  \date    March 04, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if RING_AVERAGE > 1
////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void RingAverageCons(Data *data, Ary3D (&Rs)[6], Grid *grid)
{
mynvtxstart_("RingAverageCons",SILVER);

  int    nv, i;
  int    nchunks;

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  #ifdef STAGGERED_MHD
  auto &Ax1 = grid->A[IDIR];
  auto &Bx1 = data->Vs[IDIR];
  auto &Phix1 = data->Phi[IDIR];
  auto &PhiMx1 = data->PhiM[IDIR];
  #if GEOMETRY == CYLINDRICAL
  auto &Bx3 = data->Vs[KDIR];
  auto &PhiMx3 = data->PhiM[KDIR];
  auto &Phix3 = data->Phi[KDIR];
  #elif GEOMETRY == SPHERICAL
  auto &Bx2 = data->Vs[JDIR];
  auto &Phix2 = data->Phi[JDIR];
  auto &PhiMx2 = data->PhiM[JDIR];
  #endif
  #if RING_AVERAGE_FOURIER == YES
  auto &A0x1 = data->A0[IDIR];
  auto &A0x2 = data->A0[JDIR];
  auto &A0x3 = data->A0[KDIR];

  auto &A1x1 = data->A1[IDIR];
  auto &A1x2 = data->A1[JDIR];
  auto &A1x3 = data->A1[KDIR];

  auto &B1x1 = data->B1[IDIR];
  auto &B1x2 = data->B1[JDIR];
  auto &B1x3 = data->B1[KDIR];
  #endif
  #endif

  // 0a. Check if we're close the axis. Return if processor does not own
  //     its left boundary or if csize = 1

  #if GEOMETRY == CYLINDRICAL
  if (grid->lbound[IDIR] != SINGULAR_AXIS) return;
  #elif GEOMETRY == SPHERICAL
  if (   grid->lbound[JDIR] != SINGULAR_AXIS 
      && grid->rbound[JDIR] != SINGULAR_AXIS) return;
  #endif

  RBox box_ring;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_ring);
  #pragma omp target enter data map(to: box_ring)
  #pragma acc enter data copyin(box_ring)

  // ------------------------------------------------------
  // 1. CYLINDRICAL geometry: average conservative variables 
  //    and magnetic flux over all rings. 
  //    Continue until chunk size (= csize[]) > 1.
  // ------------------------------------------------------

  #if GEOMETRY == CYLINDRICAL
  int k;

  // 1a. Convert to primitive

  ConsToPrim3D (data, 0, grid, &box_ring);

  int nphi = grid->nint[JDIR], jl;

  // 1b. Average conservative variables

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid) 
  NVAR_LOOP(nv){
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= grid->end_ring; i++){

    int csize = grid->ring_av_csize[i];
 
    nchunks = nphi / csize;
 
    #pragma acc loop vector
    for (int m = 0; m < nchunks; m++) {
      
      int jl_beg = jbeg + m*csize;
      double uav  = 0.0;
      double dVav = 0.0;

      #pragma acc loop seq 
      for (jl = jl_beg; jl < jl_beg+csize; jl++) {
        dVav += grid->dV(i,jl,k);
        uav  += data->Uc(i,jl,k,nv)*grid->dV(i,jl,k);
      }

      // Replace all zones inside a chunk with average value
      #pragma acc loop seq
      for (jl = jl_beg; jl < jl_beg+csize; jl++) data->Uc(i,jl,k,nv) = uav/dVav;
    }
  }}}

  // 1c. Average magnetic flux (Compute phiM_r and phiM_z mean)

  #ifdef STAGGERED_MHD
  auto &Ax3 = grid->A[KDIR];

  // 1c1. Apply Fourier decomposition (Eq. [16] in Zhang et al.)

  #if RING_AVERAGE_FOURIER == YES
  #pragma omp target teams distribute parallel for loop collapse(2)
  #pragma acc parallel loop collapse(2) present(data)
  for (k = kbeg; k <= kend+INCLUDE_KDIR; k++){
  for (i = ibeg; i <= iend;              i++){
    A0x1(i,k) = A1x1(i,k) = B1x1(i,k) = 0.0;
    #if INCLUDE_KDIR == YES
    A0x3(i,k) = A1x3(i,k) = B1x3(i,k) = 0.0;
    #endif
  }}

  #pragma omp target teams distribute parallel for loop collapse(2)
  #pragma acc parallel loop collapse(2) present(data)
  for (i = ibeg; i <= grid->end_ring; i ++){
  for (k = kbeg; k <= kend; k++){
    #pragma acc loop seq
    for (j = jbeg; j <= jend; j++){
      A0x1(i,k) +=   Bx1(i+1,j,k)*Ax1(i+1,j,k);
      A1x1(i,k) +=   Bx1(i+1,j,k)*Ax1(i+1,j,k)*cos(2.0*M_PI*j/nphi);
      B1x1(i,k) += 2*Bx1(i+1,j,k)*Ax1(i+1,j,k)*sin(2.0*M_PI*j/nphi);
    }
    A0x1(i,k) /= nphi;
    A1x1(i,k) /= nphi;
    B1x1(i,k) /= nphi;
  }}
  #if INCLUDE_KDIR == YES
  #pragma omp target teams distribute parallel for loop collapse(2)
  #pragma acc parallel loop collapse(2) present(data)
  for (i = ibeg; i <= grid->end_ring; i ++){
  for (k = kbeg; k <= kend+1;         k++){
    #pragma acc loop seq
    for (j = jbeg; j <= jend; j++){
      A0x3(i,k) +=   Bx3(i,j,k)*Ax3(i,j,k);
      A1x3(i,k) +=   Bx3(i,j,k)*Ax3(i,j,k)*cos(2.0*M_PI*j/nphi);
      B1x3(i,k) += 2*Bx3(i,j,k)*Ax3(i,j,k)*sin(2.0*M_PI*j/nphi);
    }
    A0x3(i,k) /= nphi;
    A1x3(i,k) /= nphi;
    B1x3(i,k) /= nphi;
  }}
  #endif
  #endif
  
  // 1c2. Compute magnetic fluxes

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, grid)
  for (i = ibeg; i <= grid->end_ring; i ++){
  for (k = kbeg; k <= kend + INCLUDE_KDIR; k++){ 
    
    int csize = grid->ring_av_csize[i];
    nchunks = nphi / csize;

    #pragma acc loop vector
    for (int m = 0; m < nchunks; m++) {
      
      int jl_beg = jbeg + m*csize;
      DIM_EXPAND(double uav_r = 0.0;,    
                                    ,
                 double uav_z = 0.0;)

      #pragma acc loop seq 
      for (jl = jl_beg; jl < jl_beg+csize; jl++) {

        DIM_EXPAND(Phix1(i+1,jl,k) =  Bx1(i+1,jl,k)*Ax1(i+1,jl,k);,
                                                                  ,
                   Phix3(i,jl,k)   =  Bx3(i,jl,k)*Ax3(i,jl,k);     )

        #if RING_AVERAGE_FOURIER == YES
        double cj = cos(2.0*M_PI*jl/nphi);
        double sj = sin(2.0*M_PI*jl/nphi);

        DIM_EXPAND(Phix1(i+1,jl,k) -=   A0x1(i,k) 
                                      + A1x1(i,k)*cj 
                                      + B1x1(i,k)*sj;,
                                                     , 
                   Phix3(i,  jl,k) -=   A0x3(i,k) 
                                      + A1x3(i,k)*cj 
                                      + B1x3(i,k)*sj;)
        #endif

        DIM_EXPAND(uav_r += Phix1(i+1,jl,k);,    
                                                   ,
                   uav_z += Phix3(i,  jl,k);)
      }

      #pragma acc loop seq
      for (jl = jl_beg; jl < jl_beg+csize; jl++) {
        DIM_EXPAND(PhiMx1(i+1,jl,k) = uav_r/csize;,
                                                  ,  
                   PhiMx3(i,  jl,k) = uav_z/csize;)
      }
    }
  }}
  #endif // STAGGERED_MHD

  #endif // GEOMETRY == CYLINDRICAL

  // ------------------------------------------------------
  // 2. Repeat for SPHERICAL geometry
  // ------------------------------------------------------

  #if GEOMETRY == SPHERICAL
  int j;

  // 2a. Convert to primitive

  ConsToPrim3D (data, 0, grid, &box_ring);

  int nphi = grid->nint[KDIR], kl;

  // 2b. Average conservative variables (Loop and average 
  //     from North to South pole)

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  NVAR_LOOP(nv){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    int csize = grid->ring_av_csize[j];
    if (csize == 1) continue; // Move to next ring
    nchunks = nphi / csize;

    #pragma acc loop vector
    for (int m = 0; m < nchunks; m++) {
      int kl_beg = kbeg + m*csize;
      double dVav = 0.0;
      double uav = 0.0;

      #pragma acc loop seq 
      for (kl = kl_beg; kl < kl_beg+csize; kl++) {
        dVav += grid->dV(i,j,kl);
        uav  += data->Uc(i,j,kl,nv)*grid->dV(i,j,kl);
      }

      #pragma acc loop seq
      for (kl = kl_beg; kl < kl_beg+csize; kl++) {
        data->Uc(i,j,kl,nv) = uav/dVav;
      }
    }
  }}}

  // 2c. Average magnetic flux

  #ifdef STAGGERED_MHD
  auto &Ax2 = grid->A[JDIR];
  double *th = grid->xC[JDIR];

  // 2c1. Apply Fourier decomposition (Eq. [16] in Zhang et al.)

  #if RING_AVERAGE_FOURIER == YES
  #pragma omp target teams distribute parallel for loop collapse(2)
  #pragma acc parallel loop collapse(2) present(data)
  for (j = jbeg; j <= jend;              j++){
  for (i = ibeg; i <= iend+INCLUDE_IDIR; i++){
    A0x2(i,j) = A1x2(i,j) = B1x2(i,j) = 0.0;
    A0x1(i,j) = A1x1(i,j) = B1x1(i,j) = 0.0;
  }}

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, grid)
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    int csize = grid->ring_av_csize[j];
    if (csize == 1) continue; // Move to next ring

    int dnp = (th[j] < 0.5*M_PI ? 1:0);  // Extra index shift at north pole

    #pragma acc loop seq
    for (k = kbeg; k <= kend; k++){
      A0x2(i,j) +=   Bx2(i,j+dnp,k)*Ax2(i,j+dnp,k);
      A1x2(i,j) +=   Bx2(i,j+dnp,k)*Ax2(i,j+dnp,k)*cos(2.0*M_PI*k/nphi);
      B1x2(i,j) += 2*Bx2(i,j+dnp,k)*Ax2(i,j+dnp,k)*sin(2.0*M_PI*k/nphi);
    }
    A0x2(i,j) /= nphi;
    A1x2(i,j) /= nphi;
    B1x2(i,j) /= nphi;
  }}

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, grid)
  for (j = jbeg; j <= jend;              j++){
  for (i = ibeg; i <= iend+INCLUDE_IDIR; i++){
    #pragma acc loop seq
    for (k = kbeg; k <= kend; k++){
      A0x1(i,j) +=   Bx1(i,j,k)*Ax1(i,j,k);
      A1x1(i,j) +=   Bx1(i,j,k)*Ax1(i,j,k)*cos(2.0*M_PI*k/nphi);
      B1x1(i,j) += 2*Bx1(i,j,k)*Ax1(i,j,k)*sin(2.0*M_PI*k/nphi);
    }
    A0x1(i,j) /= nphi;
    A1x1(i,j) /= nphi;
    B1x1(i,j) /= nphi;
  }}
  #endif
    
  // 2c2. Compute magnetic fluxes

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, grid)
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend+INCLUDE_IDIR; i++){

    int csize = grid->ring_av_csize[j];
    if (csize == 1) continue; // Move to next ring

    int dnp = (th[j] < 0.5*M_PI ? 1:0);  // Extra index shift at north pole
    
    nchunks = nphi / csize;

    #pragma acc loop vector
    for (int m = 0; m < nchunks; m++) {
      int kl_beg = kbeg + m*csize;
      double uav_r  = 0.0;
      double uav_th = 0.0;

      #pragma acc loop seq 
      for (kl = kl_beg; kl < kl_beg+csize; kl++) {

        Phix2(i,j+dnp,kl) = Bx2(i,j+dnp,kl)*Ax2(i,j+dnp,kl); 
        Phix1(i,j,    kl) = Bx1(i,j,    kl)*Ax1(i,j,kl);

        #if RING_AVERAGE_FOURIER == YES
        double ck = cos(2.0*M_PI*kl/nphi);
        double sk = sin(2.0*M_PI*kl/nphi);

        Phix2(i,j+dnp,kl) -=   A0x2(i,j) 
                             + A1x2(i,j)*ck 
                             + B1x2(i,j)*sk;  
                                                                 
        Phix1(i,j,    kl) -=   A0x1(i,j) 
                             + A1x1(i,j)*ck 
                             + B1x1(i,j)*sk;
        #endif

        uav_r  += Phix1(i,j,    kl); 
        uav_th += Phix2(i,j+dnp,kl);
      }

      #pragma acc loop seq
      for (kl = kl_beg; kl < kl_beg+csize; kl++) {
        PhiMx2(i,j+dnp,kl) = uav_th/csize; 
        PhiMx1(i,j,    kl) = uav_r/csize;
      }
    }
  }}
  #endif

  #endif  // GEOMETRY == SPHERICAL

  // ------------------------------------------------------
  // 3. Reconstruct profile in each chunk
  // ------------------------------------------------------

  // Convert to primitive

  ConsToPrim3D (data, 0, grid, &box_ring);

  #if GEOMETRY == CYLINDRICAL
  RingAverageReconstruct<JDIR>(data, Rs, grid, &box_ring);
  #elif GEOMETRY == SPHERICAL
  RingAverageReconstruct<KDIR>(data, Rs, grid, &box_ring);
  #endif

  // Convert to primitive

  ConsToPrim3D (data, 0, grid, &box_ring);
  
  #pragma omp target exit data map(delete:box_ring)
  #pragma acc exit data delete(box_ring)

mynvtxstop_();

}
////////////////////////////////////////////////////////////////////////
/*!
  - Check requirements for using RING_AVERAGE.
  - Determine the size of chunks on each ring, as a function of
  R (in polar) or \theta (spherical).
*/
////////////////////////////////////////////////////////////////////////
void RingAverageSize(Grid *grid)
{
  int csize;

  #if GEOMETRY == CYLINDRICAL

  // ------------------------------------------------------
  // 1a. Check requirements for CYLINDRICAL GEOMETRY:
  // 
  //     - NX2 must be a power of 2
  //     - NX2 must be divisible by RING_AVERAGE
  //     - phi direction cannot be parallel
  // ------------------------------------------------------ 

  int i;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int nx2 = grid->nint[JDIR];
  int j = nx2;
  while (j != 1) {
    if (j%2 != 0) {
      QuitPLUTO("! RingAverageCons(): nx2 is not a power of 2",1);
    }
    j /= 2;  
  }

  if (nx2%RING_AVERAGE != 0){
  printf ("nx2 = %d, RING_AVERAGE = %d\n", nx2, RING_AVERAGE);
    QuitPLUTO ("from RingAverageCons(): nx2 not divisible by RING_AVERAGE",1);
  }

  if (grid->nproc[JDIR] > 1){
    QuitPLUTO("from RingAverageSize(): RING_AVERAGE requires nprocs(phi) = 1",1);
  }

  // ------------------------------------------------------
  // 1b. Determine chunk size
  // ------------------------------------------------------

  if (grid->lbound[IDIR] == SINGULAR_AXIS){
    csize = RING_AVERAGE;    // Inner ring start with RING_AVERAGE

    for (i = ibeg; i <= iend; i++){
      grid->ring_av_csize[i] = csize;
      if (csize > 1) { 
        csize >>= 1;   // Same as csize /= 2 
        grid->end_ring = i;  // Last ring with csize > 1
      } 
    }
  }else{
    for (i = ibeg; i <= iend; i++) grid->ring_av_csize[i] = 1; 
  }

#elif GEOMETRY == SPHERICAL

  // ------------------------------------------------------
  //  2a. Check requirements for SPHERICAL GEOMETRY:
  // 
  //      - NX3 must be a power of 2
  //      - NX3 must be divisible by RING_AVERAGE
  //      - phi direction cannot be parallel
  // ------------------------------------------------------

  int j;
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int nx3 = grid->nint[KDIR];
  int k = nx3;
  while (k != 1) {
   if (k%2 != 0) {
     QuitPLUTO("! RingAverageCons(): nx3 is not a power of 2",1);
   }
   k /= 2;  
  }

  if (nx3%RING_AVERAGE != 0){
    QuitPLUTO ("from RingAverageCons(): nx3 not divisible by RING_AVERAGE",1);
  }

  if (grid->nproc[KDIR] > 1){
    QuitPLUTO("! RING_AVERAGE requires nprocs(phi) = 1",1);
  }

  // ------------------------------------------------------
  // 2b. Determine chunk size
  // ------------------------------------------------------

  if (grid->lbound[JDIR] == SINGULAR_AXIS){  // Close to north pole ? 

    csize = RING_AVERAGE;
    for (j = jbeg; j <= jend; j++){    
      grid->ring_av_csize[j] = csize;
      if (csize > 1) csize >>= 1;   //  Same as csize /= 2 
    }
  }else{
    for (j = jbeg; j <= jend; j++) grid->ring_av_csize[j] = 1;
  }

  if (grid->rbound[JDIR] == SINGULAR_AXIS){  // Close to south pole ? 
    csize = RING_AVERAGE;
    for (j = jend; j >= jbeg; j--){    
      grid->ring_av_csize[j] = MAX(grid->ring_av_csize[j], csize);
      if (csize > 1) csize >>= 1;   // Same as csize /= 2 
    }
  }else{
    for (j = jend; j >= jbeg; j--) {
      grid->ring_av_csize[j] = MAX(grid->ring_av_csize[j], 1);
    }
  }

#endif

}
#endif // RING_AVERAGE > 1 
