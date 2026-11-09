////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Copy 3D data into a 1D sweep buffers.

  This function extracts 1D slice ("sweep") of the 3D simulation data
  along the direction specified by the template parameter \p dir
  (IDIR, JDIR, or KDIR).  

  It copies cell variables, flags, and optional fields (gravity, staggered
  MHD components) into contiguous arrays stored inside \c data->sweep,
  preparing them for 1D reconstruction, Riemann solvers, or directional
  sweeps.
 
  The region to extract is defined by @p box. Strides and offsets are
  computed according to the grid geometry and the chosen direction.
 

  \authors A. Mignone  (andrea.mignone@unito.it)
  \date    Nov 15, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// Explicit template instantiations
template void DataToSweep<IDIR>(Data *, Grid *, RBox *);
template void DataToSweep<JDIR>(Data *, Grid *, RBox *);
template void DataToSweep<KDIR>(Data *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
  \tparam dir  Sweep direction (IDIR, JDIR, or KDIR).
 
  \param[in] data   Pointer to the simulation data structure containing
                    conservative/primitive variables and sweep buffers.
  \param[in] grid   Pointer to the grid structure containing sizes, strides,
                    and topology.
  \param[in] box    Pointer to the region of interest to be copied
                    (index bounds in i, j, k).
 
  \note This routine is heavily optimized for GPU offloading (OpenACC/OpenMP)
        and relies on pre-allocated sweep buffers.
 
  \warning The sweep arrays in @c data->sweep must be allocated with
           sufficient size before calling this function.
 
*/ ///////////////////////////////////////////////////////////////////// */
template<int dir>
void DataToSweep (Data *data, Grid *grid, RBox *box)
{
  int i,j,k,nv;
  int ni, nj;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  int nx1_tot = grid->ntot[IDIR];

  #ifdef HIGH_ORDER
    const Ary4D& Vc      = data->Vpc;   
    const Ary3D (&Vs)[6] = data->Vps;
    #pragma omp target update to(Vc,Vs)
  #else
    const Ary4D& Vc      = data->Vc;
    const Ary3D (&Vs)[6] = data->Vs;
    #pragma omp target update to(Vc,Vs)
  #endif

  ARRAY_STRIDES (grid, dir, ni, nj);

  // ------------------------------------------------------
  // 1. Copy point values to local arrays (IDIR)
  // ------------------------------------------------------

  if (dir == IDIR){
    #pragma acc parallel loop collapse(3) present(data, grid)
    #pragma omp target teams distribute parallel for collapse(3) \
                map(to: data[:1], grid[:1]) 
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = 0;    i < nx1_tot; i++){
      
      long int offset = ni*(j + nj*k);
      long int offset1 = NVAR*offset;

      // Offset vc to the correct sweep index. 
      // Note: for every pair of (j,k) indices vc() 
      // points at different memory area.

      Ary2D vc(&data->sweep->vC[offset1],ni,NVAR);
      Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

      #if SELF_GRAV
      Ary1D g_c(&data->sweep->g_c[offset],ni);
      #endif

      #ifdef STAGGERED_MHD
      CT_ASSIGN(Ary1D Bn(&data->sweep->Bn[offset],ni);  ,
                Ary1D En(&data->sweep->En[offset],ni);)
      #endif	    

      flag(i) = data->flag(i,j,k);

      #pragma acc loop seq
      NVAR_LOOP(nv) vc(i,nv) = Vc(i,j,k,nv);

      #if SELF_GRAV
      g_c(i) = data->acc_grav(IDIR,i,j,k);
      #endif

      #ifdef STAGGERED_MHD
      CT_ASSIGN(Bn(i) = Vs[BX1s](i,j,k);  ,
                En(i) = Vs[EX1s](i,j,k);)
      #endif

    }}}
  }

  // ------------------------------------------------------
  // 2. Copy point values to local arrays (JDIR)
  // ------------------------------------------------------

  #if INCLUDE_JDIR
  int nx2_tot = grid->ntot[JDIR];
  if (dir == JDIR){
    #pragma acc parallel loop collapse(3) present(data, grid)
    #pragma omp target teams distribute parallel for collapse(3) \
                map(to: data[:1], grid[:1])
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = 0;    j < nx2_tot; j++){
      long int offset  = ni*(i + nj*k);
      long int offset1 = NVAR*offset;

      // Offset vc to the correct sweep index. 
      // Note: for every pair of (i,k) indices vc() 
      // points at different memory area.

      Ary2D vc(&data->sweep->vC[offset1],ni,NVAR);
      Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

      #if SELF_GRAV
      Ary1D g_c(&data->sweep->g_c[offset],ni);
      #endif
      
      #ifdef STAGGERED_MHD
      CT_ASSIGN(Ary1D Bn(&data->sweep->Bn[offset],ni);  ,
                Ary1D En(&data->sweep->En[offset],ni);)
      #endif

      flag(j) = data->flag(i,j,k);

      #pragma acc loop  seq
      NVAR_LOOP(nv) vc(j,nv) = Vc(i,j,k,nv);

      #if SELF_GRAV
      g_c(j) = data->acc_grav(JDIR,i,j,k);
      #endif

      #ifdef STAGGERED_MHD
      CT_ASSIGN(Bn(j) = Vs[BX2s](i,j,k);  ,
                En(j) = Vs[EX2s](i,j,k);)
      #endif

    }}}
  }
  #endif

  // ------------------------------------------------------
  // 3. Copy point values to local arrays (KDIR)
  // ------------------------------------------------------

  #if INCLUDE_KDIR
  int nx3_tot = grid->ntot[KDIR];
  if (dir == KDIR){
    #pragma acc parallel loop collapse(3) present(data, grid)
    #pragma omp target teams distribute parallel for collapse(3) 
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = 0;    k < nx3_tot; k++){

      long int offset  = ni*(i + nj*j);
      long int offset1 = NVAR*offset;

      // Offset vc to the correct sweep index. 
      // Note: for every pair of (i,j) indices vc() 
      // points at different memory area.

      Ary2D vc(&data->sweep->vC[offset1],ni,NVAR);
      Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
      
      #if SELF_GRAV
      Ary1D g_c(&data->sweep->g_c[offset],ni);
      #endif
      
      #ifdef STAGGERED_MHD
      CT_ASSIGN(Ary1D Bn(&data->sweep->Bn[offset],ni);  ,
                Ary1D En(&data->sweep->En[offset],ni);)
      #endif

      flag(k) = data->flag(i,j,k);

      #pragma acc loop  seq
      NVAR_LOOP(nv) vc(k,nv) = Vc(i,j,k,nv);

      #if SELF_GRAV
      g_c(k) = data->acc_grav(KDIR,i,j,k);
      #endif

      #ifdef STAGGERED_MHD
      CT_ASSIGN(Bn(k) = Vs[BX3s](i,j,k);  ,
                En(k) = Vs[EX3s](i,j,k);)
      #endif

    }}}
  }
  #endif

  #ifndef STAGGERED_MHD
  #endif
  #if !INCLUDE_JDIR || !INCLUDE_KDIR
  #endif
}

