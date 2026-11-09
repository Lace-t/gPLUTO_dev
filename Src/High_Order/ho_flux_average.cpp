////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Apply high-order surface averaging to interface fluxes
          in the sweep direction \p dir.
  
   This function corrects the 1D interface fluxes computed during a
   directional sweep by adding transverse high-order contributions.
   The 1D fluxes are first copied into a temporary 3D array, then a
   quadrature rule based on transverse second derivatives is applied.
   Limiting is enforced when HO_LAP_LIMITER flags are present.
  
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// Explicit template instantiations
template void HO_FluxAverage<IDIR>(Data *, Grid *, RBox *, RBox *);
template void HO_FluxAverage<JDIR>(Data *, Grid *, RBox *, RBox *);
template void HO_FluxAverage<KDIR>(Data *, Grid *, RBox *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
  \tparam dir  Sweep direction (IDIR, JDIR, or KDIR).
  
  \param data       Simulation data (reads flags, reads/writes sweep fluxes).
  \param grid       Grid structure (dimensions, indexing, strides).
  \param box        Active region where corrected fluxes are required.
  \param sweep_box  Region containing the 1D fluxes from the sweep,
  typically including ghost zones.
  
  \note Uses OpenACC/OpenMP target pragmas and allocates a static
  temporary 3D flux array on first call.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HO_FluxAverage(Data *data, Grid *grid, RBox *box, RBox *sweep_box)
{ 
  int i,j,k,nv;
  int ni, nj;

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  int bibeg = sweep_box->ibeg, biend = sweep_box->iend;
  int bjbeg = sweep_box->jbeg, bjend = sweep_box->jend;
  int bkbeg = sweep_box->kbeg, bkend = sweep_box->kend;

  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  double d2fx, d2fy, d2fz, df;

  static Ary4D Flux;
  
  // ------------------------------------------------------
  // 0. Memory allocation
  // ------------------------------------------------------
  
  if (Flux.isAllocated() == 0) {
    Flux.create(nx1_tot,nx2_tot,nx3_tot,NVAR);
    #pragma acc enter data copyin(Flux)
    #pragma omp target enter data map(to:Flux)
  }

  ARRAY_STRIDES (grid, dir, ni, nj);

  // ------------------------------------------------------
  // 1. Surface quadrature for X1-Flux
  // ------------------------------------------------------

  if (dir == IDIR){

    // 1a. Copy 1D flux to 3D array

    #pragma acc parallel loop collapse(4) present(data, Flux)
    #pragma omp target teams distribute parallel for collapse(4) map(to: data[:1],Flux)
    for (nv = 0; nv < NVAR; nv++){
    for (k = bkbeg; k <= bkend; k++){
    for (j = bjbeg; j <= bjend; j++){
    for (i = ibeg;  i <= iend+1; i++){ 
      long int offset1 = NVAR*ni*(j + nj*k);
      Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

      Flux(i,j,k,nv) = flux(i,nv);
    }}}}

    // 1b. Apply high-order quadrature rule to x1-flux

    d2fx = d2fy = d2fz = df = 0.0;
    #pragma acc parallel loop collapse(4) present(data, Flux)
    #pragma omp target teams distribute parallel for collapse(4) map(to: data[:1],Flux)
    for (nv = 0; nv < NVAR; nv++){
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){
      long int offset1 = NVAR*ni*(j + nj*k);
      Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

      DIM_EXPAND(                                           ,
                 d2fy = HO_Average4D<JDIR>(Flux,i,j,k,nv);  ,
                 d2fz = HO_Average4D<KDIR>(Flux,i,j,k,nv);)

      df = (d2fy + d2fz)/24.0;

      #if HO_LAP_LIMITER != NO
      if (   (data->flag(i-1,j,k) & FLAG_HO_LAP_LIMITER) 
          || (data->flag(i,j,k)   & FLAG_HO_LAP_LIMITER) ) df = 0.0;
      #endif

      flux(i,nv) = Flux(i,j,k,nv) + df;
    }}}}
  }

  // ------------------------------------------------------
  // 2. Surface quadrature for X2-Flux
  // ------------------------------------------------------

  if (dir == JDIR){

    // 2a. Copy 1D flux to 3D array

    #pragma acc parallel loop collapse(4) present(data, Flux)
    #pragma omp target teams distribute collapse(4) map(to: data[:1], Flux)
    for (nv = 0; nv < NVAR; nv++){
    for (k = bkbeg;   k <= bkend; k++){
    for (i = bibeg;   i <= biend; i++){
    for (j = bjbeg; j <= bjend+1; j++){

      long int offset1 = NVAR*ni*(i + nj*k);
      Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
      Flux(i,j,k,nv) = flux(j,nv);
    }}}}

    // 2b. Apply high-order quadrature rule to x2-flux

    d2fx = d2fy = d2fz = df = 0.0;
    #pragma acc parallel loop collapse(4) present(data, Flux)
    #pragma omp target teams distribute collapse(4) map(to: data[:1], Flux)
    for (nv = 0; nv < NVAR; nv++){
    for (k = kbeg;   k <= kend; k++){
    for (i = ibeg;   i <= iend; i++){
    for (j = jbeg; j <= jend+1; j++){

      long int offset1 = NVAR*ni*(i + nj*k);
      Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

      DIM_EXPAND(d2fx = HO_Average4D<IDIR>(Flux,i,j,k,nv);  ,
                                                              ,
                  d2fz = HO_Average4D<KDIR>(Flux,i,j,k,nv);   
      )

      df = (d2fx + d2fz)/24.0;

      #if HO_LAP_LIMITER != NO
      if (   (data->flag(i,j-1,k) & FLAG_HO_LAP_LIMITER) 
          || (data->flag(i,j,k)   & FLAG_HO_LAP_LIMITER) ) df = 0.0;
      #endif

      flux(j,nv) = Flux(i,j,k,nv) + df;
    }}}}
  }

  // ------------------------------------------------------
  // 3. Surface quadrature for X3-Flux
  // ------------------------------------------------------

  if (dir == KDIR){

    // 3a. Copy 1D flux to 3D array

    #pragma acc parallel loop collapse(4) present(data, Flux)
    #pragma omp target teams distribute collapse(4) //map(present: data[:1], Flux)
    for (nv = 0; nv < NVAR; nv++){
    for (j = bjbeg;   j <= bjend; j++){
    for (i = bibeg;   i <= biend; i++){
    for (k = bkbeg; k <= bkend+1; k++){
      long int offset1 = NVAR*ni*(i + nj*j);
      Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
      
      //#pragma acc loop vector
      Flux(i,j,k,nv) = flux(k,nv);
    }}}}

    // 3b. Apply high-order quadrature rule to x3-flux

    d2fx = d2fy = d2fz = df = 0.0;
    #pragma acc parallel loop collapse(4) present(data, Flux)
    #pragma omp target teams distribute collapse(4) //map(present: data[:1], Flux)
    for (nv = 0; nv < NVAR; nv++){
    for (j = jbeg;   j <= jend; j++){
    for (i = ibeg;   i <= iend; i++){
    for (k = kbeg; k <= kend+1; k++){
      long int offset1 = NVAR*ni*(i + nj*j);
      Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

      DIM_EXPAND(d2fx = HO_Average4D<IDIR>(Flux,i,j,k,nv);  ,
                  d2fy = HO_Average4D<JDIR>(Flux,i,j,k,nv);  , 
      )
      df = (d2fx + d2fy)/24.0;

      #if HO_LAP_LIMITER != NO
      if (   (data->flag(i,j,k-1) & FLAG_HO_LAP_LIMITER) 
          || (data->flag(i,j,k)   & FLAG_HO_LAP_LIMITER) ) df = 0.0;
      #endif

      flux(k,nv) = Flux(i,j,k,nv) + df;
    }}}}
  }
}