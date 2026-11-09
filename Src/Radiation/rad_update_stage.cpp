////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Substepping integration routines.

  Integrate the inhomogeneous radiation transport equations
  by carrying IMEX integrations with dt=dt_rad until t+dt_hydro
  is reached.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 8, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if RADIATION_NR

double getMax(Data *, const RBox *);

////////////////////////////////////////////////////////////////////////
/*!
  \param[in,out] data     pointer to Data structure
  \param[out]    dUc      cell-centered conservative variables increment
                          for radiation variables
  \param[in,out] Dts_rad  pointer to time step structure
  \param[in]     grid     pointer to Grid structure
  \param[in]     box      active domain integration box
  \param[in]     dt       time step for the current update step
*/
////////////////////////////////////////////////////////////////////////
void UpdateRadStage (Data *data, Ary4D &dUc, timeStep *Dts_rad, Grid *grid,
                     RBox *box, double dt)
{
  mynvtxstart_("URS", BLUE);

  int i, j, k;
  int ni, nj;
  int ntr = 0;
  const int comps = RADIATION_NEQS - 1;

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  RBox sweep_box;

  #ifdef HIGH_ORDER
  ntr = 1;
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

  // 1c. Obtain radiation flux at i+1/2 interfaces

  mynvtxstart_("RadRiemannSolver#1", BLUE);
  data->radiationRiemannSolver(data, Dts_rad, grid, &sweep_box, IDIR);
  mynvtxstop_();

  // 1d. Obtain surface-averaged flux

  #ifdef HIGH_ORDER
  HO_FluxAverage<IDIR>(data, grid, box, &sweep_box);
  #endif

  // 1e. Compute right hand side for x1-dir.

  RadRightHandSide(data, Dts_rad, grid, dt, box, IDIR);

  #pragma acc parallel loop collapse(3) present(data, dUc)
  #pragma omp target teams distribute parallel for collapse(3) map(to: data[:1], dUc)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    long int offset1 = ni*(j + nj*k);
    long int offset = NVAR*offset1;

    Ary1D cmax(&Dts_rad->cmax[offset1], ni);
    Ary2D rhs(&data->sweep->rhs[offset], ni, NVAR);

    dUc(i,j,k,ENR) = rhs(i,ENR);
    dUc(i,j,k,ENG) = rhs(i,ENG);
    for (int n = 0; n < comps; n++){
      dUc(i,j,k,FR1+n) = rhs(i,FR1+n);
      dUc(i,j,k,MX1+n) = rhs(i,MX1+n);
    }

    if (g_intStage == 1){
      data->C_dt(i,j,k) = 0.5*(cmax(i) + cmax(i+1))/grid->dx[IDIR][i];
    }
  }}}

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

  // 2c. Obtain radiation flux at j+1/2 interfaces

  mynvtxstart_("RadRiemannSolver#2", BLUE);
  data->radiationRiemannSolver(data, Dts_rad, grid, &sweep_box, JDIR);
  mynvtxstop_();

  // 2d. Obtain surface-averaged flux

  #ifdef HIGH_ORDER
  HO_FluxAverage<JDIR>(data, grid, box, &sweep_box);
  #endif

  // 2e. Compute right hand side for x2-dir

  RadRightHandSide(data, Dts_rad, grid, dt, box, JDIR);

  #pragma acc parallel loop collapse(3) present(data, dUc)
  #pragma omp target teams distribute parallel for collapse(3) map(to: data[:1], dUc)
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
  for (j = jbeg; j <= jend; j++){
    long int offset1 = ni*(i + nj*k);
    long int offset = NVAR*offset1;

    Ary1D cmax(&Dts_rad->cmax[offset1], ni);
    Ary2D rhs(&data->sweep->rhs[offset], ni, NVAR);

    dUc(i,j,k,ENR) += rhs(j,ENR);
    dUc(i,j,k,ENG) += rhs(j,ENG);
    for (int n = 0; n < comps; n++){
      dUc(i,j,k,FR1+n) += rhs(j,FR1+n);
      dUc(i,j,k,MX1+n) += rhs(j,MX1+n);
    }

    if (g_intStage == 1){
      double dl = grid->dx[JDIR][j];
      #if GEOMETRY == CYLINDRICAL || GEOMETRY == SPHERICAL
      dl *= fabs(grid->rC[i]);
      #endif
      data->C_dt(i,j,k) += 0.5*(cmax(j) + cmax(j+1))/dl;
    }
  }}}

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

  // 3c. Obtain radiation flux at k+1/2 interfaces

  mynvtxstart_("RadRiemannSolver#3", BLUE);
  data->radiationRiemannSolver(data, Dts_rad, grid, &sweep_box, KDIR);
  mynvtxstop_();

  // 3d. Obtain surface-averaged flux

  #ifdef HIGH_ORDER
  HO_FluxAverage<KDIR>(data, grid, box, &sweep_box);
  #endif

  // 3e. Compute right hand side for x3-dir

  RadRightHandSide(data, Dts_rad, grid, dt, box, KDIR);

  #pragma acc parallel loop collapse(3) present(data, dUc)
  #pragma omp target teams distribute parallel for collapse(3) map(to: data[:1], dUc)
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  for (k = kbeg; k <= kend; k++){
    long int offset1 = ni*(i + nj*j);
    long int offset = NVAR*offset1;

    Ary1D cmax(&Dts_rad->cmax[offset1], ni);
    Ary2D rhs(&data->sweep->rhs[offset], ni, NVAR);

    dUc(i,j,k,ENR) += rhs(k,ENR);
    dUc(i,j,k,ENG) += rhs(k,ENG);
    for (int n = 0; n < comps; n++){
      dUc(i,j,k,FR1+n) += rhs(k,FR1+n);
      dUc(i,j,k,MX1+n) += rhs(k,MX1+n);
    }

    if (g_intStage == 1){
      double dl = grid->dx[KDIR][k];
      #if GEOMETRY == SPHERICAL
      dl *= fabs(grid->rC[i]*grid->sC[j]);
      #endif
      data->C_dt(i,j,k) += 0.5*(cmax(k) + cmax(k+1))/dl;
    }
  }}}

  #pragma acc exit data delete(sweep_box)
  #pragma omp target exit data map(delete:sweep_box)

#endif // INCLUDE_KDIR

  // 4. Internal Boundary Reset

  #if INTERNAL_BOUNDARY == YES
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, dUc, box)
  BOX_LOOP(*box,i,j,k){
    if (data->flag(i,j,k) & FLAG_INTERNAL_BOUNDARY){
      dUc(i,j,k,ENR) = 0.0;
      dUc(i,j,k,ENG) = 0.0;

      for (int n = 0; n < comps; n++){
        dUc(i,j,k,FR1+n) = 0.0;
        dUc(i,j,k,MX1+n) = 0.0;
      }

      #if INTERNAL_BOUNDARY_CFL == NO
      data->C_dt(i,j,k) = 0.0;
      #endif
    }
  }
  #endif

  // 5. Compute inverse dt

  if (g_intStage == 1){
    double dim_fact = static_cast<double>(INCLUDE_IDIR + INCLUDE_JDIR + INCLUDE_KDIR);
    double inv_dt_ = getMax(data, box);
    Dts_rad->invDt_hyp = inv_dt_ / dim_fact;

    // #pragma acc update device(Dts_rad)
  }

  mynvtxstop_();
}

#endif // RADIATION_NR
