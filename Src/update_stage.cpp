////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Single stage integration for RK time stepping.

  Compute the solution increment in conservative form,
  \f[
     \begin{array}{l}
      \delta U_c \equiv \Delta t R_c(V) \\[0.5em]
      \delta V_s \equiv \Delta t R_s(V)
     \end{array}
  \f]
  where \f$ \delta U_c \f$ and \f$ \delta V_s \f$ are 3D arrays
  containing, respectively, the right hand sides (increments) for 
  conservative and staggered variables.
  Here \c R_c(V) and \c R_s(V) contain contributions from all directions.

  When the integrator stage is the first one (predictor), this function
  also computes the maximum of inverse time steps for hyperbolic and
  parabolic terms (if the latters are included explicitly).

  \authors A. Mignone       (andrea.mignone@unito.it)\n
           V. Berta         (vittoria.berta@unito.it)\n
           M. Rossazza      (marco.rossazza@unito.it)\n
           C. Zanni         (zanni@oato.inaf.it) \n
           D. Melon Fuksman (fuksman@mpia.de)
  \date    May 7, 2026

*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

double getMax(Data *, const RBox *);
double getMaxPar(Data *, const RBox *);

////////////////////////////////////////////////////////////////////////
/*!
  Compute the right hand side of the system of conservation laws.
  
  \param [in]     data    pointer to PLUTO data structure
  \param [out]    dUc     cell-centered conservative variables increment
  \param [out]    dUs     staggered variables increment
  \param [in]     grid    pointer to Grid structure
  \param [in,out] Dts     time step structure
  \param [in]     box     active domain integration box
  \param [in]     dt      time step 
*/
////////////////////////////////////////////////////////////////////////
// Preserve the original entry point used by other physics/time integrators.
// Only the PCP driver can prove that a stage has no correction faces.
void UpdateStage(Data *data, Ary4D &dUc, Ary3D (&dUs)[6], Grid *grid,
                 timeStep *Dts, RBox *box, double dt)
{
  UpdateStage(data, dUc, dUs, grid, Dts, box, dt, true);
}

void UpdateStage (Data *data, Ary4D &dUc, Ary3D (&dUs)[6], Grid *grid, 
                  timeStep *Dts, RBox *box, double dt, bool pcp_flux_possible)
{
  mynvtxstart_("US",RED);
  int nv, i, j, k;
  int ni, nj;
  int ntr = 0;  // Default number of transverse zones 

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  RBox sweep_box;

  #pragma acc data present(data)
  //#pragma omp target enter data map(to: data[:1])

  // 0b. Set the stencil in the transverse direction

#ifdef STAGGERED_MHD
  #if CT_EMF_AVERAGE == ARITHMETIC || CT_EMF_AVERAGE == CT_CONTACT   
  ntr = 1;
  #else
  ntr = grid->nghost[IDIR];
  #endif
#else
  #ifdef HIGH_ORDER     
  ntr = 1;
  #endif
#endif

  // ------------------------------------------------------
  // 1. X1-Update
  // ------------------------------------------------------
 
  ARRAY_STRIDES (grid, IDIR, ni, nj);
  RBoxDefine(ibeg, iend,
             jbeg - ntr*INCLUDE_JDIR, jend + ntr*INCLUDE_JDIR,
             kbeg - ntr*INCLUDE_KDIR, kend + ntr*INCLUDE_KDIR, &sweep_box);

  #pragma acc enter data copyin(sweep_box)
  #pragma omp target enter data map(to:sweep_box)

  // 1a. Copy 3D arrays into 1D buffers

  DataToSweep<IDIR>(data, grid, &sweep_box);

  // 1b. Reconstruct and compute x1-flux      

  if (data->use_flat) ReconstructFlat<IDIR> (data, grid, &sweep_box);
  else                Reconstruct<IDIR>     (data, grid, &sweep_box);

  // 1c. Obtain flux (hyperbolic + parabolic) at i+1/2 interfaces

  mynvtxstart_("RiemannSolver#1",SILVER);
  #if PHYSICS == RMHD && PCP_NR == YES
  RMHD_PCP_LimitFaces(data, grid, &sweep_box, IDIR);
  #endif
  data->fluidRiemannSolver(data, Dts, grid, &sweep_box, IDIR);
  #if PHYSICS == RMHD && PCP_FOFC == YES
  if (pcp_flux_possible) {
    RMHD_PCP_Flux(data, Dts, grid, &sweep_box, IDIR);
  }
  #endif
  #if RESISTIVITY == EXPLICIT
  ParabolicFlux_T (data, Dts, grid, &sweep_box, IDIR);
  #endif 
  mynvtxstop_();

  #if RADIATION_RELATIVISTIC
  mynvtxstart_("RadRiemannSolver#1",BLUE);
  data->radiationRiemannSolver(data, Dts, grid, &sweep_box, IDIR);
  mynvtxstop_();
  #endif

  // 1d. Obtain surface-averaged flux 

  #ifdef HIGH_ORDER
  HO_FluxAverage<IDIR>(data, grid, box, &sweep_box);
  #endif 

  // 1e. Compute right hand side for x1-dir.

  RightHandSide (data, Dts, grid, dt, box, IDIR);

  #if (PHYSICS == GRMHD && GR_METRIC != CARTESIAN) 
  RightHandSideSource(data, grid, dt, box, IDIR);
  #endif

  #pragma acc parallel loop collapse(4) present(data, dUc)
  #pragma omp target teams distribute parallel for collapse(4)  map(to: data[:1],dUc)
  NFLUID_LOOP(nv){
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary2D rhs(&data->sweep->rhs[offset1],ni,NVAR);

    dUc(i,j,k,nv) = rhs(i,nv);
  }}}}
  #pragma acc exit data delete(sweep_box)
  #pragma omp target exit data map(delete:sweep_box)

  // ------------------------------------------------------
  // 2. X2-Update
  // ------------------------------------------------------

#if INCLUDE_JDIR

  ARRAY_STRIDES (grid, JDIR, ni, nj);
  RBoxDefine(ibeg - ntr*INCLUDE_IDIR, iend + ntr*INCLUDE_IDIR,
             jbeg, jend,
             kbeg - ntr*INCLUDE_KDIR, kend + ntr*INCLUDE_KDIR, &sweep_box);

  #pragma acc enter data copyin(sweep_box)
  #pragma omp target enter data map(to:sweep_box)

  // 2a. Copy 3D arrays into 1D buffers

  DataToSweep<JDIR>(data, grid, &sweep_box);

  // 2b. Reconstruct and compute x2-flux

   if (data->use_flat) ReconstructFlat<JDIR> (data, grid, &sweep_box);
  else                Reconstruct<JDIR>     (data, grid, &sweep_box);

  // 2c. Obtain flux (hyperbolic + parabolic) at j+1/2 interfaces

  mynvtxstart_("RiemannSolver#2",SILVER);
  #if PHYSICS == RMHD && PCP_NR == YES
  RMHD_PCP_LimitFaces(data, grid, &sweep_box, JDIR);
  #endif
  data->fluidRiemannSolver(data, Dts, grid, &sweep_box, JDIR);
  #if PHYSICS == RMHD && PCP_FOFC == YES
  if (pcp_flux_possible) {
    RMHD_PCP_Flux(data, Dts, grid, &sweep_box, JDIR);
  }
  #endif
  #if RESISTIVITY == EXPLICIT
  ParabolicFlux_T (data, Dts, grid, &sweep_box, JDIR);
  #endif 
  mynvtxstop_();

  #if RADIATION_RELATIVISTIC
  mynvtxstart_("RadRiemannSolver#2",BLUE);
  data->radiationRiemannSolver(data, Dts, grid, &sweep_box, JDIR);
  mynvtxstop_();
  #endif

  // 2d. Obtain surface-averaged flux 

  #ifdef HIGH_ORDER
  HO_FluxAverage<JDIR>(data, grid, box, &sweep_box);
  #endif

  // 2e. Compute right hand side for x2-dir

  RightHandSide (data, Dts, grid, dt, box, JDIR);

  #if (PHYSICS == GRMHD && GR_METRIC != CARTESIAN) 
  RightHandSideSource(data, grid, dt, box, JDIR);
  #endif
  
  #pragma acc parallel loop collapse(4) present(data, dUc)
  #pragma omp target teams distribute parallel for collapse(4) map(to: data[:1],dUc)
  NFLUID_LOOP(nv){
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
  for (j = jbeg; j <= jend; j++){
    long int offset = NVAR*ni*(i + nj*k);

    Ary1D cs2(&data->sweep->cs2[offset],ni);
    Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);

    dUc(i,j,k,nv) += rhs(j,nv);
  }}}}
  #pragma acc exit data delete(sweep_box)
  #pragma omp target exit data map(delete:sweep_box)

#endif // INCLUDE_JDIR

  // ------------------------------------------------------
  // 3. X3-Update
  // ------------------------------------------------------

#if INCLUDE_KDIR

  ARRAY_STRIDES (grid, KDIR, ni, nj);
  RBoxDefine(ibeg - ntr*INCLUDE_IDIR, iend + ntr*INCLUDE_IDIR,
             jbeg - ntr*INCLUDE_JDIR, jend + ntr*INCLUDE_JDIR,
             kbeg, kend,                    &sweep_box);

  #pragma acc enter data copyin(sweep_box)
  #pragma omp target enter data map(to:sweep_box)

  // 3a. Copy 3D array into 1D buffers

  DataToSweep<KDIR>(data, grid, &sweep_box);

  // 3b. Reconstruct and compute x3-flux

  if (data->use_flat) ReconstructFlat<KDIR> (data, grid, &sweep_box);
  else                Reconstruct<KDIR>     (data, grid, &sweep_box);

  // 3c. Obtain flux (hyperbolic + parabolic) at k+1/2 interfaces

  mynvtxstart_("RiemannSolver#3",SILVER);
  #if PHYSICS == RMHD && PCP_NR == YES
  RMHD_PCP_LimitFaces(data, grid, &sweep_box, KDIR);
  #endif
  data->fluidRiemannSolver(data, Dts, grid, &sweep_box, KDIR);
  #if PHYSICS == RMHD && PCP_FOFC == YES
  if (pcp_flux_possible) {
    RMHD_PCP_Flux(data, Dts, grid, &sweep_box, KDIR);
  }
  #endif
  #if RESISTIVITY == EXPLICIT
  ParabolicFlux_T (data, Dts, grid, &sweep_box, KDIR);
  #endif   
  mynvtxstop_();

  #if RADIATION_RELATIVISTIC
  mynvtxstart_("RadRiemannSolver#3",BLUE);
  data->radiationRiemannSolver(data, Dts, grid, &sweep_box, KDIR);
  mynvtxstop_();
  #endif

  // 3d. Obtain surface-averaged flux 

  #ifdef HIGH_ORDER
  HO_FluxAverage<KDIR>(data, grid, box, &sweep_box);
  #endif 

  // 3e. Compute right hand side for x3-dir

  RightHandSide (data, Dts, grid, dt, box, KDIR);
  #pragma acc parallel loop collapse(4) present(data, dUc)
  #pragma omp target teams distribute parallel for collapse(4) map(to: data[:1],dUc)
  NFLUID_LOOP(nv){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  for (k = kbeg; k <= kend; k++){
    long int offset = NVAR*ni*(i + nj*j);
    Ary1D cs2(&data->sweep->cs2[offset],ni);
    Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
   
    dUc(i,j,k,nv) += rhs(k,nv);
  }}}}
  #pragma acc exit data delete(sweep_box)
  #pragma omp target exit data map(delete:sweep_box)
#endif // INCLUDE_KDIR

  // 4. Add source terms

  #ifdef HIGH_ORDER
  Sources (data, dUc, dt, grid, box); //questo
  #endif

  // 5. Internal Boundary Reset
  
  #if INTERNAL_BOUNDARY == YES
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, dUc)
  BOX_LOOP(*box,i,j,k){
    if (data->flag(i,j,k) & FLAG_INTERNAL_BOUNDARY){
      NVAR_LOOP(nv) dUc(i,j,k,nv) = 0.0;
      #if INTERNAL_BOUNDARY_CFL == NO
      data->C_dt(i,j,k) = 0.0;
      #endif
    }
  }
  #endif

  // 6. CT Update

  #ifdef STAGGERED_MHD
  CT_ComputeEMF(data, grid, box); 
  #if RESISTIVITY == EXPLICIT
  CT_ResistiveEMF (data, 1, grid, box);
  #endif   
  CT_Update(data, dUs, dt, grid, box); 
  #if CT_ENERGY_CORRECTION > 1
  CT_EnergyCorrection (data, dUc, dt, grid, box);
  #endif
  #endif

  // 7. Compute inverse dt

  if (g_intStage == 1){
    double dim_fact = static_cast<double>(INCLUDE_IDIR + INCLUDE_JDIR + INCLUDE_KDIR);
    double inv_dt_ = getMax(data, box);
    Dts->invDt_hyp = inv_dt_ / dim_fact;

    //add parabolic timestep

    #if RESISTIVITY
    double inv_dt_par_ = getMaxPar(data, box);
    Dts->invDt_par = inv_dt_par_ / dim_fact;
    #endif
    
    // #pragma acc update device(Dts)
  }

  mynvtxstop_();
}

////////////////////////////////////////////////////////////////////////
/*! 
  Note: This function was modified on 18 Mar 2025 since after 
  some changes, the ResRMHD/KH conf #4 and #5 gave some errors.
*/
////////////////////////////////////////////////////////////////////////
double getMax(Data *data, const RBox *box)
{
  int i,j,k;
  double local = -1.0;

  #pragma acc data copy(local)
  {
    #pragma acc parallel loop collapse(3) reduction(max:local) present(local, data, box)
    #pragma omp target teams distribute parallel for collapse(3) reduction(max:local) map(to: data[:1],box[:1])    
    BOX_LOOP(*box, i,j,k) {
      local = MAX(local, data->C_dt(i,j,k));
    }
  }

  return local;
}

double getMaxPar(Data *data, const RBox *box)
{
  int i,j,k;
  double local = -1.0;

  #pragma acc data copy(local)
  {
    #pragma acc parallel loop collapse(3) reduction(max:local) present(local, data, box)
    #pragma omp target teams distribute parallel for collapse(3) reduction(max:local) map(to: data[:1],box[:1])      
    BOX_LOOP(*box, i,j,k) {
      local = MAX(local, data->C_dtpar[IDIR](i,j,k));
      local = MAX(local, data->C_dtpar[JDIR](i,j,k));
      local = MAX(local, data->C_dtpar[KDIR](i,j,k));
    }
  }
  return local;
}
