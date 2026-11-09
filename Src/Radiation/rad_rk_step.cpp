////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Substepping integration routines.

  Integrate the inhomogeneous radiation transport equations
  by carrying IMEX integrations with dt=dt_rad until t+dt_hydro
  is reached.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 9, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

double getMax(Data *, const RBox *);

// ***********************************************************************
#if RADIATION_IMEX_SSP2

////////////////////////////////////////////////////////////////////////
/*!
  Apply a single IMEX source function: U += dt*S.
  Used by the IMEX-SSP2(2,2,2) scheme.

  \param[in,out] data  pointer to Data structure
  \param[in]     S     source-term array
  \param[in]     box   active integration box
  \param[in]     dt    time increment
*/
////////////////////////////////////////////////////////////////////////
void AddRadSource1 (Data *data, Grid *grid, Ary4D &S, RBox *box, double dt)
{
  int i, j, k, nv;
  const int comps = RADIATION_NEQS - 1;
  double src[NVAR];

  #if RADIATION_NR
  RunConfig &run_config = *(data->run_config);

  #pragma omp target teams distribute parallel for collapse(3) private(src)
  #pragma acc parallel loop collapse(3) present(data, S, run_config, grid, box) private(src[:NVAR])
  BOX_LOOP(*box, i, j, k){
    #pragma acc loop seq
    NVAR_LOOP(nv) src[nv] = 0.0;
    src[ENG] = run_config.rad_c * dt * S(i,j,k,0);
    src[ENR] = -run_config.rad_reduced_c * dt * S(i,j,k,0);
    for (int n = 0; n < comps; n++){
      src[MX1+n] = dt * S(i,j,k,n+1);
      src[FR1+n] = -run_config.rad_reduced_c * dt * S(i,j,k,n+1);
    }
    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    OrthoToCons(src, metric);
    #else
    (void)grid;
    #endif
    #pragma acc loop seq
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) += src[nv];
  }
  #else

  #pragma omp target teams distribute parallel for collapse(3) private(src)
  #pragma acc parallel loop collapse(3) present(data, S, grid, box) private(src[:NVAR])
  BOX_LOOP(*box, i, j, k){
    #pragma acc loop seq
    NVAR_LOOP(nv) src[nv] = 0.0;
    src[ENG] = dt * S(i,j,k,0);
    src[ENR] = -dt * S(i,j,k,0);
    for (int n = 0; n < comps; n++){
      src[MX1+n] = dt * S(i,j,k,n+1);
      src[FR1+n] = -dt * S(i,j,k,n+1);
    }
    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    OrthoToCons(src, metric);
    #else
    (void)grid;
    #endif
    #pragma acc loop seq
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) += src[nv];
  }
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Apply two IMEX source functions: U += dt1*S1 + dt2*S2.
  Used by the IMEX-SSP2(2,2,2) scheme.

  \param[in,out] data  pointer to Data structure
  \param[in]     S1    first source-term array
  \param[in]     S2    second source-term array
  \param[in]     box   active integration box
  \param[in]     dt1   time step for S1 integration
  \param[in]     dt2   time step for S2 integration
*/
////////////////////////////////////////////////////////////////////////
void AddRadSource2 (Data *data, Grid *grid, Ary4D &S1, Ary4D &S2, RBox *box,
                    double dt1, double dt2)
{
  int i, j, k, nv;
  const int comps = RADIATION_NEQS - 1;
  double src[NVAR];

  #if RADIATION_NR
  RunConfig &run_config = *(data->run_config);

  #pragma omp target teams distribute parallel for collapse(3) private(src)
  #pragma acc parallel loop collapse(3) present(data, S1, S2, run_config, grid, box) private(src[:NVAR])
  BOX_LOOP(*box, i, j, k){
    #pragma acc loop seq
    NVAR_LOOP(nv) src[nv] = 0.0;
    double se = dt1*S1(i,j,k,0) + dt2*S2(i,j,k,0);
    src[ENG] = run_config.rad_c * se;
    src[ENR] = -run_config.rad_reduced_c * se;
    for (int n = 0; n < comps; n++){
      double sm = dt1*S1(i,j,k,n+1) + dt2*S2(i,j,k,n+1);
      src[MX1+n] = sm;
      src[FR1+n] = -run_config.rad_reduced_c * sm;
    }
    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    OrthoToCons(src, metric);
    #else
    (void)grid;
    #endif
    #pragma acc loop seq
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) += src[nv];
  }
  #else

  #pragma omp target teams distribute parallel for collapse(3) private(src)
  #pragma acc parallel loop collapse(3) present(data, S1, S2, grid, box) private(src[:NVAR])
  BOX_LOOP(*box, i, j, k){
    #pragma acc loop seq
    NVAR_LOOP(nv) src[nv] = 0.0;
    double se = dt1*S1(i,j,k,0) + dt2*S2(i,j,k,0);
    src[ENG] = se;
    src[ENR] = -se;
    for (int n = 0; n < comps; n++){
      double sm = dt1*S1(i,j,k,n+1) + dt2*S2(i,j,k,n+1);
      src[MX1+n] = sm;
      src[FR1+n] = -sm;
    }
    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    OrthoToCons(src, metric);
    #else
    (void)grid;
    #endif
    #pragma acc loop seq
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) += src[nv];
  }
  #endif
}

#endif // RADIATION_IMEX_SSP2

// ***********************************************************************
#if RADIATION_NR

// Weight factors for 2nd stage of RK integrators.

#if TIME_STEPPING == RK2
 #define w0  0.5
 #define wc  0.5
#elif TIME_STEPPING == RK3
 #define w0 0.75
 #define wc 0.25
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Advance radiation transport equations by a single time step using
  IMEX-Runge-Kutta integrators (Euler/RK2/RK3).
  When RADIATION_IMEX_SSP2 is enabled, the IMEX-SSP2(2,2,2) scheme of
  Pareschi & Russo (2005) is used (requires TIME_STEPPING == RK2).

  \param[in,out] data     pointer to Data structure
  \param[in]     dt       time step for this substep
  \param[in,out] Dts_rad  pointer to the radiation-substep time structure
  \param[in]     grid     pointer to Grid structure

  \return 0 on success.
*/
////////////////////////////////////////////////////////////////////////
int AdvanceRadStep (Data *data, double dt, timeStep *Dts_rad, Grid *grid)
{
  int i, j, k, nv;
  const double one_third = 1.0/3.0;
  const int comps = RADIATION_NEQS - 1;
  static bool first_call = true;

  static Ary4D Rc;

  #if RADIATION_IMEX_SSP2
   #if TIME_STEPPING != RK2
   QuitPLUTO("! AdvanceRadStep(): RADIATION_IMEX_SSP2 requires TIME_STEPPING == RK2\n",1);
   #endif
  static Ary4D Srad1, Srad2;
  const double rk_gamma = 0.29289321881;
  double dt_rk1 = dt * rk_gamma;
  double dt_rk2 = dt * (1.0 - 2.0*rk_gamma);
  double dt_rk3 = dt * 0.5*(1.0 - rk_gamma);
  #endif

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  static RBox box;
  if(first_call){
    RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box);
    #pragma acc enter data copyin(box)
    #pragma omp target enter data map(to: box)
    first_call = false;
  }

  const Ary4D &U0 = data->Uc0;

  // -------------------------------------------------------
  // 0. Allocate IMEX-SSP2 source arrays on first call
  // -------------------------------------------------------

  if (!Rc.isAllocated()){
    int nx1 = grid->ntot[IDIR], nx2 = grid->ntot[JDIR], nx3 = grid->ntot[KDIR];
    Rc.create(nx1, nx2, nx3, NVAR);
    #pragma acc enter data copyin(Rc)
    #pragma omp target enter data map(to: Rc)
  }

  #if RADIATION_IMEX_SSP2
  if (!Srad1.isAllocated()){
    int nx1 = grid->ntot[IDIR], nx2 = grid->ntot[JDIR], nx3 = grid->ntot[KDIR];
    Srad1.create(nx1, nx2, nx3, NVAR);
    Srad2.create(nx1, nx2, nx3, NVAR);
    #pragma acc enter data copyin(Srad1, Srad2)
    #pragma omp target enter data map(to: Srad1, Srad2)
  }
  #endif

  // -------------------------------------------------------
  // 1. Predictor step (EULER, RK2, RK3)
  // -------------------------------------------------------

  // 1a. Set boundary conditions.

  g_intStage = 1;
  Boundary(data, 0, grid);

  // 1b. Convert primitive to conservative, save initial stage.

  PrimToCons3D(data, 0, grid, &box);
  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc0(i,j,k,nv) = data->Uc(i,j,k,nv);
  }

  // 1c. First implicit step (IMEX-SSP2 only).

  #if RADIATION_IMEX_SSP2
  ConsToPrim3D(data, 0, grid, &box);
  RadStep3D(data, Srad1, grid, dt_rk1);
  Boundary(data, 0, grid);
  #endif

  // 1d. Explicit radiation update.

  UpdateRadStage(data, Rc, Dts_rad, grid, &box, dt);

  #if RADIATION_IMEX_SSP2
  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = U0(i,j,k,nv);
  }

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, Rc)
  BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,ENR) += Rc(i,j,k,ENR);
    data->Uc(i,j,k,ENG) += Rc(i,j,k,ENG);
    for (int n = 0; n < comps; n++){
      data->Uc(i,j,k,FR1+n) += Rc(i,j,k,FR1+n);
      data->Uc(i,j,k,MX1+n) += Rc(i,j,k,MX1+n);
    }
  }
  #else
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, Rc)
  BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,ENR) += Rc(i,j,k,ENR);
    data->Uc(i,j,k,ENG) += Rc(i,j,k,ENG);
    for (int n = 0; n < comps; n++){
      data->Uc(i,j,k,FR1+n) += Rc(i,j,k,FR1+n);
      data->Uc(i,j,k,MX1+n) += Rc(i,j,k,MX1+n);
    }
  }
  #endif

  // 1e. Implicit step (second for IMEX-SSP2, only for standard IMEX1).

  #if RADIATION_IMEX_SSP2
  AddRadSource1(data, grid, Srad1, &box, dt_rk2);
  ConsToPrim3D(data, 0, grid, &box);
  RadStep3D(data, Srad2, grid, dt_rk1);
  #else
  ConsToPrim3D(data, 0, grid, &box);
  RadStep3D(data, grid, dt);
  #endif

  // -------------------------------------------------------
  // 2. Corrector step (RK2, RK3)
  // -------------------------------------------------------

#if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)

  // 2a. Set boundary conditions.

  g_intStage = 2;
  Boundary(data, 0, grid);

  // 2b. Explicit radiation update.

  UpdateRadStage(data, Rc, Dts_rad, grid, &box, dt);

  #if RADIATION_IMEX_SSP2
  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = 0.5*U0(i,j,k,nv) + 0.5*data->Uc(i,j,k,nv);
  }

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, Rc)
  BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,ENR) += 0.5*Rc(i,j,k,ENR);
    data->Uc(i,j,k,ENG) += 0.5*Rc(i,j,k,ENG);
    for (int n = 0; n < comps; n++){
      data->Uc(i,j,k,FR1+n) += 0.5*Rc(i,j,k,FR1+n);
      data->Uc(i,j,k,MX1+n) += 0.5*Rc(i,j,k,MX1+n);
    }
  }
  #else
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, Rc)
  BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,ENR) += Rc(i,j,k,ENR);
    data->Uc(i,j,k,ENG) += Rc(i,j,k,ENG);
    for (int n = 0; n < comps; n++){
      data->Uc(i,j,k,FR1+n) += Rc(i,j,k,FR1+n);
      data->Uc(i,j,k,MX1+n) += Rc(i,j,k,MX1+n);
    }
  }
  #endif

  // 2c. Implicit step (IMEX1 only; IMEX-SSP2 skips this).

  #if !RADIATION_IMEX_SSP2
  ConsToPrim3D(data, 0, grid, &box);
  RadStep3D(data, grid, dt);
  #endif

  // 2d. RK averaging: Uc = w0*Uc0 + wc*Uc.

  #if !RADIATION_IMEX_SSP2
  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = w0*U0(i,j,k,nv) + wc*data->Uc(i,j,k,nv);
  }
  #endif

  // 2e. IMEX-SSP2: apply combined source.

  #if RADIATION_IMEX_SSP2
  AddRadSource2(data, grid, Srad1, Srad2, &box, dt_rk1, dt_rk3);
  #endif

  // 2f. Convert to primitive variables.

  ConsToPrim3D(data, 0, grid, &box);

#endif // TIME_STEPPING == RK2 || RK3

  // -------------------------------------------------------
  // 3. Last corrector step (RK3 only; no IMEX-SSP2 here)
  // -------------------------------------------------------

#if TIME_STEPPING == RK3

  // 3a. Set boundary conditions.

  g_intStage = 3;
  Boundary(data, 0, grid);

  // 3b. Explicit radiation update.

  UpdateRadStage(data, Rc, Dts_rad, grid, &box, dt);

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, Rc)
  BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,ENR) += Rc(i,j,k,ENR);
    data->Uc(i,j,k,ENG) += Rc(i,j,k,ENG);
    for (int n = 0; n < comps; n++){
      data->Uc(i,j,k,FR1+n) += Rc(i,j,k,FR1+n);
      data->Uc(i,j,k,MX1+n) += Rc(i,j,k,MX1+n);
    }
  }

  // 3c. Implicit step.

  ConsToPrim3D(data, 0, grid, &box);
  RadStep3D(data, grid, dt);

  // 3d. RK averaging: Uc = 1/3*Uc0 + 2/3*Uc.

  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(data)
  NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = one_third*(U0(i,j,k,nv) + 2.0*data->Uc(i,j,k,nv));
  }

  // 3e. Convert to primitive variables.

  ConsToPrim3D(data, 0, grid, &box);

#endif // TIME_STEPPING == RK3

  return 0;
}
#endif // RADIATION_NR
