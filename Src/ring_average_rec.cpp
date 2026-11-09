////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Ring-average reconstruction methods 

  Reconstruct values inside each chunk using standard limited 
  reconstruction. This yields Um, Up and UM.  
  Use these value for a parabolic fit Up(x) and redefine cell averages 
  in the original grid:
 
       Initial Chunk:   |----U----|----U----|----U----|----U----|
       Reconstruct:     |Um-----------------U-----------------Up|
 
       Using parabolic fit Up(),to reassign values:
 
                        |--Up(1)--|--Up(2)--|--Up(3)--|--Up(4)--|

  See ring_average.cpp.

  \b Reference:
    - "Conservative averaging-reconstruction techniques (Ring Average)
       for 3-D finite-volume MHD solvers with axis singularity",
       Zhang et al, JCP (2019) 376, 276-294
       
  \authors A. Mignone (andrea.mignone@unito.it)\n
           M. Rossazza 
  \date    March 06, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "States/states.hpp"

#if RING_AVERAGE > 1

// Explicit template instantiations
template void RingAverageReconstruct<IDIR> (Data *, Ary3D (&Rs)[6], Grid *, RBox *);
template void RingAverageReconstruct<JDIR> (Data *, Ary3D (&Rs)[6], Grid *, RBox *);
template void RingAverageReconstruct<KDIR> (Data *, Ary3D (&Rs)[6], Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void RingAverageReconstruct (Data *data, Ary3D (&Rs)[6], Grid *grid, RBox *box)
{
mynvtxstart_("PLMRingAverageReconstruct",GREEN);
  
  int i, j, k, nv, m;
  int ti, tj, tk;
  int csize;        // Size of a single chunk 
  int mbeg, mend;   // Initial and final active indices on chunked grid
  int nchunks;      // Number of chunks (= grid size on reduced grid) 
  int ngh  = GetNghost();
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  RBox box_loc;

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_loc);

  #if GEOMETRY == CYLINDRICAL
  int nphi = grid->nint[JDIR];

  #elif GEOMETRY == SPHERICAL
  int nphi = grid->nint[KDIR];
  #endif
  int ni, nj;

  #if RING_AVERAGE_REC_PRIM == YES
  RunConfig &run_config = *(data->run_config);
  #endif

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(&box_loc, ibeg, iend, jbeg, jend, kbeg, kend);
  ARRAY_STRIDES (grid, dir, ni, nj);

  #if RING_AVERAGE_REC == LINEAR
  double dvap, dvam, dva;
  #endif

  // #if RING_AVERAGE_REC == 1
  // return;
  // #endif

  // ------------------------------------------------------
  // 1. Assign chunk averages to 1D arrays
  // 
  //    ua(m) = U(j)
  //    va(m) = V(j)
  // 
  //    m is a contiguous index running on chunks (not cells), 
  //    starting at the same ibeg as the global grid in order to 
  //    have the same number of ghost zones.
  //    For csize == 4 (for instance):
  //
  //    i = 3, m = 3
  //    i = 7, m = 4 
  //    ...
  //    i = 127, m = 34
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (nv = 0; nv < NVAR; nv++){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){

    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j]; // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k]; // t_j == k
    #endif
    if (csize == 1) continue; // No need to reconstruct 
    
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D va(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D ua(&data->sweep->uC[offset1],ni,NVAR);

    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
    mbeg = ibeg; //dbeg;
    mend = ibeg + nchunks - 1;

    // m is the macrocell index, a contiguous index 
    // running on chunks (not cells)

    #pragma acc loop vector
    for (m = mbeg; m <= mend; m++) {
      i = ibeg + (m-mbeg)*csize; 
      if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
      if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

      ua(m,nv) = data->Uc(ti,tj,tk,nv);   
      va(m,nv) = data->Vc(ti,tj,tk,nv);
    }
    
    // 1c. Set periodic b.c. on the reduced grid

    #pragma acc loop seq
    for (m = 0; m < mbeg; m++) {
      ua(m,nv) = ua(m+nchunks,nv);
      va(m,nv) = va(m+nchunks,nv);
    }

    #pragma acc loop seq
    for (m = mend+1; m <= mend + ngh; m++){
      ua(m,nv) = ua(m-nchunks,nv);
      va(m,nv) = va(m-nchunks,nv);
    }
  }}}

  // ------------------------------------------------------
  // 2. Reconstruct left and right interface values 
  //    vam() and vap()  on the reduced grid (chunk):
  //
  //    Original:     |-A1-|-A2-|-B1-|-B2-|-C1-|-C2-|
  //    Chunk:        |----A----|----B----|----C----|
  //    Rec. values:            Bm       Bp
  // ------------------------------------------------------
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (nv = 0; nv < NVAR; nv++){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    
    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j]; // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k]; // t_j == k
    #endif
    if (csize == 1) continue; // No need to reconstruct further

    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D va(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D ua(&data->sweep->uC[offset1],ni,NVAR);
    Ary2D vap(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D vam(&data->sweep->fR[offset1],ni,NVAR);

    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
    mbeg = ibeg; //dbeg;
    mend = ibeg + nchunks - 1;

    #if RING_AVERAGE_REC == FLAT
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      #if RING_AVERAGE_REC_PRIM == NO
      vap(m,nv) = ua(m,nv);
      vam(m,nv) = ua(m,nv);
      #else
      vap(m,nv) = va(m,nv);
      vam(m,nv) = va(m,nv);
      #endif
    }
    #elif RING_AVERAGE_REC == LINEAR
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      #if RING_AVERAGE_REC_PRIM == NO
      dvap = ua(m+1,nv) - ua(m,nv);
      dvam = ua(m,nv) - ua(m-1,nv);
      dva  = VANLEER_LIMITER(dvap, dvam);
      vap(m,nv) = ua(m,nv) + 0.5*dva;
      vam(m,nv) = ua(m,nv) - 0.5*dva;
      #else
      dvap = va(m+1,nv) - va(m,nv);
      dvam = va(m,nv) - va(m-1,nv);
      dva  = VANLEER_LIMITER(dvap, dvam);
      vap(m,nv) = va(m,nv) + 0.5*dva;
      vam(m,nv) = va(m,nv) - 0.5*dva;
      #endif
    }
    #elif RING_AVERAGE_REC == PARABOLIC
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      double w[5];
      #if RING_AVERAGE_REC_PRIM == NO
      w[4] = ua(m+2,nv);
      w[3] = ua(m+1,nv);
      w[2] = ua(m,  nv);
      w[1] = ua(m-1,nv);
      w[0] = ua(m-2,nv);
      #else
      w[4] = va(m+2,nv);
      w[3] = va(m+1,nv);
      w[2] = va(m,  nv);
      w[1] = va(m-1,nv);
      w[0] = va(m-2,nv);
      #endif

      double vp, vm; 
      PPM_States(w+2, vm, vp);
      vap(m,nv) = vp;
      vam(m,nv) = vm;
    }
    #elif RING_AVERAGE_REC == MP5
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      double w[5];
      #if RING_AVERAGE_REC_PRIM == NO
      w[4] = ua(m+2,nv);
      w[3] = ua(m+1,nv);
      w[2] = ua(m,  nv);
      w[1] = ua(m-1,nv);
      w[0] = ua(m-2,nv);
      #else
      w[4] = va(m+2,nv);
      w[3] = va(m+1,nv);
      w[2] = va(m,  nv);
      w[1] = va(m-1,nv);
      w[0] = va(m-2,nv);
      #endif

      vap(m,nv) = MP5_States (&w[2], +1);
      vam(m,nv) = MP5_States (&w[2], -1);
    }
    #else
    #error Invalid RING_AVERAGE_REC value
    #endif

  }}}

  #if RING_AVERAGE_REC_PRIM == YES
  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){

    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j]; // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k]; // t_j == k
    #endif
    if (csize == 1) continue; // No need to reconstruct further

    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D vap(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D vam(&data->sweep->fR[offset1],ni,NVAR);

    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
    mbeg = ibeg; //dbeg;
    mend = ibeg + nchunks - 1;

    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      Metric metric;
      double vloc[NVAR], uloc[NVAR];

      // Note: tj (CYLINDRICAL) and tk (SPHERICAL) should never be employed
      // since we're averaging along the azimuthal direction

      if constexpr (dir == JDIR) {ti = j; tj = -1; tk = k;}
      if constexpr (dir == KDIR) {ti = j; tj = k; tk = -1;}

      NVAR_LOOP(nv) vloc[nv] = vap(m,nv); 
      PrimToCons(vloc, uloc, run_config);
      MetricSet(metric, grid->rC[ti], grid->sC[tj]);
      OrthoToCons (uloc, metric);
      NVAR_LOOP(nv) vap(m,nv) = uloc[nv];

      NVAR_LOOP(nv) vloc[nv] = vam(m,nv); 
      PrimToCons(vloc, uloc, run_config);
      MetricSet(metric, grid->rC[ti], grid->sC[tj]);
      OrthoToCons (uloc, metric);
      NVAR_LOOP(nv) vam(m,nv) = uloc[nv];
    }  
  }}
  #endif

  // ------------------------------------------------------
  // 3. Use va, vap and vam to fit a unique parabola through
  //    the chunk and replace center for cells of the 
  //    original grid. This step has to be done on 
  //    conservative variables.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (nv = 0; nv < NVAR; nv++){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    long int offset1 = NVAR*ni*(j + nj*k);
    Ary2D va(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D ua(&data->sweep->uC[offset1],ni,NVAR);
    Ary2D vap(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D vam(&data->sweep->fR[offset1],ni,NVAR);

    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j]; // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k]; // t_j == k
    #endif
    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
    mbeg = ibeg; 
    mend = ibeg + nchunks - 1;
    if (csize == 1) continue; // No need to reconstruct further
  
    // 3a. Loop on original grid

    double urec;
    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){
      
      m = (i - ibeg)/csize + ibeg; // Macrocell index

      // Local grid index inside chunk,  1 <= ck <= csize 
      int ck  = (i - ibeg)%csize + 1;  

      // Replace zone-center 

      double fm = vam(m,nv), fp = vap(m,nv), fM = ua(m,nv);
      // double xi = double(ck-0.5)/double(csize);
      urec =    fm 
             + (3.0*fM - 2.0*fm - fp)*(2.0*ck - 1.0)/csize 
             + (-2.0*fM + fm + fp)*(3.0*ck*ck - 3.0*ck + 1.0)/(csize*csize);

      if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
      if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

      data->Uc(ti,tj,tk,nv) = urec;
    }

  }}}

  ConsToPrim3D (data, 0, grid, box);

  // Now take care of face-centered fields

  #ifdef STAGGERED_MHD
  #if GEOMETRY == CYLINDRICAL
  CT_RingAverageReconstruct<JDIR>(data, Rs, grid, box);
  #elif GEOMETRY == SPHERICAL
  CT_RingAverageReconstruct<KDIR>(data, Rs, grid, box);
  #endif
  #else
  #endif

mynvtxstop_();

}

////////////////////////////////////////////////////////////////////////
/*!
   Update magnetic field with the emf previously defined (it is zero 
   everywhere except around the chunks). 
   This will consistently replace the magnetic field with its average
   or reconstructed value
*/
////////////////////////////////////////////////////////////////////////
#ifdef STAGGERED_MHD
// Explicit template instantiations
template void CT_RingAverageReconstruct<IDIR> (Data *, Ary3D (&Rs)[6], Grid *, RBox *);
template void CT_RingAverageReconstruct<JDIR> (Data *, Ary3D (&Rs)[6], Grid *, RBox *);
template void CT_RingAverageReconstruct<KDIR> (Data *, Ary3D (&Rs)[6], Grid *, RBox *);
template<int dir>

void CT_RingAverageReconstruct (Data *data, Ary3D (&Rs)[6], Grid *grid, RBox *box)
{
mynvtxstart_("CT_RingAverageReconstruct",GREEN);
  
  int i, j, k, nv, il, m;
  int ti, tj, tk;
  int csize;        // Size of a single chunk  
  int mbeg, mend;   // Initial and final active indices on chunked grid    
  int nchunks;      // Number of chunks (= grid size on reduced grid) 
  int ngh  = GetNghost();

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  int *ntot = grid->ntot;

  auto &Ex1e = data->emf->Ee[IDIR];
  auto &Ex2e = data->emf->Ee[JDIR];
  auto &Ex3e = data->emf->Ee[KDIR];

  auto &Phix1 = data->Phi[IDIR];
  auto &PhiMx1 = data->PhiM[IDIR];
  #if GEOMETRY == CYLINDRICAL
  auto &Phix3 = data->Phi[KDIR];
  auto &PhiMx3 = data->PhiM[KDIR];
  #elif GEOMETRY == SPHERICAL
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
 
  double *dx1 = grid->dx[IDIR];
  #if GEOMETRY == CYLINDRICAL
  double *dx3 = grid->dx[KDIR];
  #elif GEOMETRY == SPHERICAL
  double *dx2 = grid->dx[JDIR];
  #endif

  RBox box_loc;

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_loc);

  #if GEOMETRY == CYLINDRICAL
  int nphi = grid->nint[JDIR];
  #elif GEOMETRY == SPHERICAL
  int nphi   = grid->nint[KDIR];
  double *th = grid->xC[JDIR];
  #endif
  int ni, nj;

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(&box_loc, ibeg, iend, jbeg, jend, kbeg, kend);
  ARRAY_STRIDES (grid, dir, ni, nj);

  #if RING_AVERAGE_CT_REC == LINEAR
  double dvap, dvam, dva;
  #endif

  #pragma acc parallel loop collapse(3) present(data, grid)
  for (k = 0; k < ntot[KDIR]; k++){
  for (j = 0; j < ntot[JDIR]; j++){
  for (i = 0; i < ntot[IDIR]; i++){
    DIM_EXPAND(                     ,
               Ex3e(i,j,k) = 0.0;   ,
               Ex1e(i,j,k) = 0.0;
               Ex2e(i,j,k) = 0.0;)
  }}}

  // ------------------------------------------------------
  // 1. Assign chunk averages to 1D arrays
  // 
  //    ua(m) = U(j)
  //    va(m) = V(j)
  // 
  //    m is a contiguous index running on chunks (not cells), 
  //    starting at the same ibeg as the global grid in order to 
  //    have the same number of ghost zones.
  // ------------------------------------------------------

  int nfields = DIMENSIONS - 1;

  int dk = INCLUDE_KDIR*(GEOMETRY == CYLINDRICAL);  // Extra point for edge-
  int dj = INCLUDE_IDIR*(GEOMETRY == SPHERICAL);    // centered emf

  #pragma omp target teams distribute parallel for collapse(2)
  #pragma acc parallel loop collapse(2) present(data, grid)
  for (k = kbeg; k <= kend + dk; k++){
  for (j = jbeg; j <= jend + dj; j++){
    
    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j];  // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k];  // t_j == k
    int dnp = (th[k] < 0.5*M_PI ? 1:0);  // Extra index shift at north pole
    // dnp = 0; // for south pole
    // if (k < ntot[JDIR]/2) dnp = 1;  // for north pole
    #endif
    if (csize == 1) continue; // No need to reconstruct further
    
    long int offset1 = nfields*ni*(j + nj*k);
    Ary2D va(&data->sweep->uC[offset1],ni,nfields);

    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
    mbeg = ibeg; //dbeg;
    mend = ibeg + nchunks - 1;

    // m is the macrocell index, a contiguous index 
    // running on chunks (not cells)

    #pragma acc loop vector
    for (m = mbeg; m <= mend; m++) {
      i = ibeg + (m-mbeg)*csize; 

      #if GEOMETRY == CYLINDRICAL 
      {
        ti = j; tj = i; tk = k;
        DIM_EXPAND(                             ,
                   va(m,0) = PhiMx1(ti+1,tj,tk);,
                   va(m,1) = PhiMx3(ti,tj,tk);   )
      }
      #endif

      #if GEOMETRY == SPHERICAL
      {
        ti = j; tj = k; tk = i;
        DIM_EXPAND(                             ,
                   va(m,0) = PhiMx2(ti,tj+dnp,tk);,
                   va(m,1) = PhiMx1(ti,tj,tk);   )
      }
      #endif

    }

    // 1c. Set periodic b.c. on the reduced grid

    #pragma acc loop seq
    for (m = 0; m < mbeg; m++) {
      DIM_EXPAND(                          ,
                 va(m,0) = va(m+nchunks,0);,
                 va(m,1) = va(m+nchunks,1); )
    }

    #pragma acc loop seq
    for (m = mend+1; m <= mend + ngh; m++) {
      DIM_EXPAND(                          ,
                 va(m,0) = va(m-nchunks,0);,
                 va(m,1) = va(m-nchunks,1); )
    }

  }}

  // ------------------------------------------------------
  // 2. Reconstruct left and right interface values 
  //    vam() and vap()  on the reduced grid (chunk):
  //
  //    Original:     |-A1-|-A2-|-B1-|-B2-|-C1-|-C2-|
  //    Chunk:        |----A----|----B----|----C----|
  //    Rec. values:            Bm       Bp
  // ------------------------------------------------------
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (nv = 0; nv < nfields; nv++){
  for (k = kbeg; k <= kend + dk; k++){
  for (j = jbeg; j <= jend + dj; j++){
    long int offset1 = nfields*ni*(j + nj*k);
    Ary2D va(&data->sweep->uC[offset1],ni,nfields);
    Ary2D vap(&data->sweep->fL[offset1],ni,nfields);
    Ary2D vam(&data->sweep->fR[offset1],ni,nfields);

    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j]; // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k]; // t_j == k
    #endif
    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
    mbeg = ibeg; //dbeg;
    mend = ibeg + nchunks - 1;
    if (csize == 1) continue; // No need to reconstruct further

    #if RING_AVERAGE_CT_REC == FLAT
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      vap(m,nv) = va(m,nv);
      vam(m,nv) = va(m,nv);
    }
    #elif RING_AVERAGE_CT_REC == LINEAR
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      dvap = va(m+1,nv) - va(m,nv);
      dvam = va(m,nv) - va(m-1,nv);
      dva  = VANLEER_LIMITER(dvap, dvam);
      vap(m,nv) = va(m,nv) + 0.5*dva;
      vam(m,nv) = va(m,nv) - 0.5*dva;
    }
    #elif RING_AVERAGE_CT_REC == PARABOLIC
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend+1; m++){
      double w[5];
      w[4] = va(m+2,nv);
      w[3] = va(m+1,nv);
      w[2] = va(m,  nv);
      w[1] = va(m-1,nv);
      w[0] = va(m-2,nv);

      double vp, vm; 
      PPM_States(w+2, vm, vp);
      vap(m,nv) = vp;
      vam(m,nv) = vm;
    }
    #elif RING_AVERAGE_CT_REC == MP5
    #pragma acc loop vector
    for (m = mbeg-1; m <= mend; m++){
      double w[5];
      w[4] = va(m+2,nv);
      w[3] = va(m+1,nv);
      w[2] = va(m,  nv);
      w[1] = va(m-1,nv);
      w[0] = va(m-2,nv);

      vap(m,nv) = MP5_States (&w[2], +1);
      vam(m,nv) = MP5_States (&w[2], -1);
    }
    #else
    #error Invalid RING_AVERAGE_CT_REC value
    #endif

  }}}

  // ------------------------------------------------------
  // 3. Use va, vap and vam to fit a unique parabola through
  //    the chunk and obtain the central value for each cell
  //    of the original grid.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (nv = 0; nv < nfields; nv++){
  for (k = kbeg; k <= kend + dk; k++){
  for (j = jbeg; j <= jend + dj; j++){
    
    #if GEOMETRY == CYLINDRICAL
    csize = grid->ring_av_csize[j]; // t_i == j
    #elif GEOMETRY == SPHERICAL
    csize = grid->ring_av_csize[k]; // t_j == k
    int dnp = (th[k] < 0.5*M_PI ? 1:0);  // Extra index shift at north pole
    #endif
    if (csize == 1) continue; // No need to reconstruct 
    
    long int offset1 = nfields*ni*(j + nj*k);
    Ary2D va(&data->sweep->uC[offset1],ni,nfields);
    Ary2D vap(&data->sweep->fL[offset1],ni,nfields);
    Ary2D vam(&data->sweep->fR[offset1],ni,nfields);

    nchunks = nphi/csize; //  >> (NPHI_AVERAGE - (g_i-ibeg));
  
    // 3a. Loop on original grid

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){
      #if GEOMETRY == CYLINDRICAL 
      ti = j; tj = i; tk = k;
      #elif GEOMETRY == SPHERICAL 
      ti = j; tj = k; tk = i;
      #endif
      #if RING_AVERAGE_FOURIER == YES
      double ci = cos(2.0*M_PI*i/nphi);
      double si = sin(2.0*M_PI*i/nphi);
      #endif
      
      m = (i - ibeg)/csize + ibeg; // Macrocell index

      // Local grid index inside chunk,  1 <= ck <= csize 
      int ck  = (i - ibeg)%csize + 1;  

      // double xp = ck/(double)csize;
      // double xm = (ck-1.0)/(double)csize;

      // reconstruct central value

      double fm = vam(m,nv), fp = vap(m,nv), fM = va(m,nv);
      // double xi = double(ck-0.5)/double(csize);

      // Map nv to actual array component (nv can only be 0 or 1 and it 
      // refers to the staggered field component moving away from the axis
      // (i.e., Br in CYLINDRICAL coords, Bth in SPHERICAL coords)

      double phi_mean =     fm 
                        + (3.0*fM - 2.0*fm - fp)*(2.0*ck - 1.0)/csize 
                        + (-2.0*fM + fm + fp)*(3.0*ck*ck - 3.0*ck + 1.0)
                           /(csize*csize);

      #if GEOMETRY == CYLINDRICAL 
      if (nv == 0) PhiMx1(ti+1,tj,tk) = phi_mean;
      else         PhiMx3(ti,tj,  tk) = phi_mean;
      #if RING_AVERAGE_FOURIER == YES
      if (nv == 0) {
        PhiMx1(ti+1,tj,tk) += A0x1(ti,tk) + A1x1(ti,tk)*ci + B1x1(ti,tk)*si;
        Phix1 (ti+1,tj,tk) += A0x1(ti,tk) + A1x1(ti,tk)*ci + B1x1(ti,tk)*si;
      } else {
        PhiMx3(ti,  tj,tk) += A0x3(ti,tk) + A1x3(ti,tk)*ci + B1x3(ti,tk)*si;
        Phix3 (ti,  tj,tk) += A0x3(ti,tk) + A1x3(ti,tk)*ci + B1x3(ti,tk)*si;
      }
      #endif
      #elif GEOMETRY == SPHERICAL
      if (nv == 0) PhiMx2(ti,tj+dnp,tk) = phi_mean;
      else         PhiMx1(ti,tj,    tk) = phi_mean;
      #if RING_AVERAGE_FOURIER == YES
      if (nv == 0) {
        PhiMx2(ti,tj+dnp,tk) += A0x2(ti,tj) + A1x2(ti,tj)*ci + B1x2(ti,tj)*si;
        Phix2 (ti,tj+dnp,tk) += A0x2(ti,tj) + A1x2(ti,tj)*ci + B1x2(ti,tj)*si;
      } else {
        PhiMx1(ti,tj,    tk) += A0x1(ti,tj) + A1x1(ti,tj)*ci + B1x1(ti,tj)*si;
        Phix1 (ti,tj,    tk) += A0x1(ti,tj) + A1x1(ti,tj)*ci + B1x1(ti,tj)*si;
      }
      #endif
      #endif
    }

    // 3c. Compute Electric fields using magnetic flux

    #pragma acc loop vector
    for (m = 0; m < nchunks; m++) {
      int il_beg = ibeg + m*csize;
      
      double dPhi0 = 0.0;
      double dPhi1 = 0.0;

      #pragma acc loop seq 
      for (il = il_beg; il < il_beg+csize-1; il++) {

        double invalpha;

        #if GEOMETRY == CYLINDRICAL 
        ti = j; tj = il; tk = k;
        double dr   = dx1[ti];
        double dz   = (INCLUDE_KDIR ? dx3[tk]:1.0);
        #elif GEOMETRY == SPHERICAL
        ti = j; tj = k; tk = il;
        double rdth = grid->xL[IDIR][ti]*dx2[tj];
        double dr   = dx1[ti];
        #endif

        if (nv == 0) {  
          #if GEOMETRY == CYLINDRICAL
          dPhi0   += Phix1(ti+1,tj,tk) - PhiMx1(ti+1,tj,tk);
          invalpha = 1.0/dz;                                                  
          Ex3e(ti+1,tj+1,tk) = invalpha*dPhi0;
          #elif GEOMETRY == SPHERICAL       
          dPhi0   += Phix2(ti,tj+dnp,tk) - PhiMx2(ti,tj+dnp,tk);
          invalpha = 1.0/dr;
          Ex1e(ti,tj+dnp,tk+1) = invalpha*dPhi0;
          #endif
        } else if (nv == 1) {  
          #if GEOMETRY == CYLINDRICAL
          dPhi1   += Phix3(ti,tj,tk) - PhiMx3(ti,tj,tk);
          invalpha = 1.0/dr;
          Ex1e(ti,tj+1,tk)   = -invalpha*dPhi1;
          #elif GEOMETRY == SPHERICAL
          dPhi1   += Phix1(ti,tj,tk) - PhiMx1(ti,tj,tk);
          invalpha = 1.0/rdth;
          Ex2e(ti,tj,tk+1) = -invalpha*dPhi1;
          #endif
        }
      }
    }
  }}}

  // ------------------------------------------------------
  // 4. Upate magnetic field with the emf previously defined 
  //    (it is zero everywhere except around the chunks).
  //    This will consistently replace the magnetic field 
  //    with its average value
  // ------------------------------------------------------
  
  // return to matching sweep ordering
  RBoxPermuteIndices<IDIR>(&box_loc, ibeg, iend, jbeg, jend, kbeg, kend);
  ARRAY_STRIDES (grid, dir, ni, nj);
  #pragma acc enter data copyin(box_loc)
  #pragma omp target enter data map(to: box_loc)

  #ifdef STAGGERED_MHD
  CT_Update(data, Rs, 1.0, grid, &box_loc); 
  #pragma acc parallel loop collapse(3) present(data, Rs) 
  BOX_STAG_LOOP(box_loc,i,j,k) {
    CT_EXPAND(nv, data->Vs[nv](i,j,k) += Rs[nv](i,j,k);)
  }

  CT_AverageStaggeredFields (data, 1, &box_loc, grid);
  #endif

  // Convert to primitive

  ConsToPrim3D (data, 0, grid, &box_loc);
  #pragma acc exit data delete(box_loc)
  #pragma omp target exit data map(delete: box_loc)

mynvtxstop_();

}
#endif // STAGGERED_MHD

#endif // RING_AVERAGE > 1 
