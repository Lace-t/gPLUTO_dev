////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Compute the right hand side of the nonrelativistic radiation
         transport equations.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 8, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if RADIATION_NR

////////////////////////////////////////////////////////////////////////
/*!
  \param[in,out] data  pointer to Data structure
  \param[in]     Dts   pointer to time step structure
  \param[in]     grid  pointer to Grid structure
  \param[in]     dt    time increment
  \param[in]     box   pointer to RBox structure
  \param[in]     dir   sweep direction

  \return This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
void RadRightHandSide (Data *data, timeStep *Dts, const Grid *grid,
                       double dt, RBox *box, int dir)
{
  int i, j, k, nv;
  int ni, nj;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  double *x1  = grid->xC[IDIR];
  double *dx1 = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];
  double *x1m = grid->xL[IDIR];
  #if GEOMETRY != CARTESIAN
  const Ary3D &Ax1 = grid->A[IDIR];
  const Ary3D &Ax2 = grid->A[JDIR];
  const Ary3D &Ax3 = grid->A[KDIR];
  #endif
  #if GEOMETRY == SPHERICAL
  double *sC = grid->sC, *sL = grid->sL;
  #endif
  [[maybe_unused]] RunConfig &run_config = *(data->run_config);
  static double *FA;

  // 0. Allocate memory

  #if GEOMETRY != CARTESIAN
  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];
  if (FA == NULL) {
    FA = ARRAY_1D(size_tot*NVAR, double);
  }
  #endif

  ARRAY_STRIDES (grid, dir, ni, nj);

  // ------------------------------------------------------
  // I. Compute radiation rhs in the x1-direction
  // ------------------------------------------------------

  if (dir == IDIR){

    #if GEOMETRY != CARTESIAN

    // I1. Build area-weighted radiation flux

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){

      long int offset = NVAR*ni*(j + nj*k);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      #pragma acc loop seq
      NRAD_LOOP(nv) fA(i,nv) = flux(i,nv)*Ax1(i,j,k);

      #if GEOMETRY == CYLINDRICAL
      fA(i,iFRPHI) *= fabs(x1m[i]);
      #elif GEOMETRY == SPHERICAL
      IF (COVARIANT_MX2, fA(i,FR2) *= fabs(x1m[i]));
      fA(i,FR3) *= fabs(x1m[i]*sC[j]);
      #endif

    }}}
    #endif // GEOMETRY != CARTESIAN

    // I2. Compute radiation right hand side

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){

      long int offset1 = ni*(j + nj*k);
      long int offset  = NVAR*offset1;

      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary2D fA(&FA[offset],ni,NVAR);
      #endif

      double v[NVAR];
      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(i,nv);

      // I2a. Radiation flux divergence

      #pragma acc loop seq
      NVAR_LOOP(nv) rhs(i,nv) = 0.0;

      #if GEOMETRY == CARTESIAN
      {
        double dtdx1 = dt/dx1[i];
        #pragma acc loop seq
        NRAD_LOOP(nv) rhs(i,nv) = -dtdx1*(flux(i+1,nv) - flux(i,nv));
      }
      #else
      {
        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NRAD_LOOP(nv) rhs(i,nv) = -dtdV*(fA(i+1,nv) - fA(i,nv));
      }
      #endif

      // I2b. v/c source terms

      #if RADIATION_IMPLICIT_NR
      {
        Ary2D vl(&data->sweep->vL[offset],ni,NVAR);
        Ary2D vr(&data->sweep->vR[offset],ni,NVAR);
        double vmid[NVAR];
        #pragma acc loop seq
        NVAR_LOOP(nv) vmid[nv] = 0.5*(vl(i+1,nv) + vr(i,nv));
        double src[NVAR] = {0.0};
        AddRelTerms(vmid, src, dt, run_config);
        #if GEOMETRY != CARTESIAN
        {
          Metric metric;
          #if GEOMETRY == SPHERICAL
          MetricSet(metric, fabs(x1[i]), sC[j]);
          #else
          MetricSet(metric, fabs(x1[i]), 0.0);
          #endif
          OrthoToCons(src, metric);
        }
        #endif
        #pragma acc loop seq
        NVAR_LOOP(nv) rhs(i,nv) += src[nv];
      }
      #endif

      // I2b'. Irradiation
      #if IRRADIATION
      rhs(i, FIR) = 0.0;
      #if (!RADIATION_IMPLICIT_NR) || RADIATION_IMEX_SSP2
      IF_ENERGY(rhs(i,ENG) -= dt*v[FIR];)
      #endif
      #endif

      // I2c. Geometric source terms

      #if GEOMETRY == CYLINDRICAL
      {
        double Sr = EddTensor(v, iFRPHI, iFRPHI) * v[ENR];
        rhs(i,iFRR) += dt*Sr/fabs(x1[i])*run_config.rad_reduced_c;
      }
      #elif GEOMETRY == SPHERICAL
      {
        double Pthth = EddTensor(v, iFRTH, iFRTH) * v[ENR];
        double Pphph = EddTensor(v, iFRPHI, iFRPHI) * v[ENR];

        rhs(i,iFRR) += dt*(Pthth + Pphph)/fabs(x1[i])*run_config.rad_reduced_c;
      }
      #endif

    }}}

  // ------------------------------------------------------
  // J. Compute radiation rhs in the x2-direction
  // ------------------------------------------------------

  }else if (dir == JDIR){

    #if GEOMETRY != CARTESIAN

    // J1. Build area-weighted radiation flux

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend+1; j++){

      long int offset = NVAR*ni*(i + nj*k);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      #pragma acc loop seq
      NRAD_LOOP(nv) fA(j,nv) = flux(j,nv)*Ax2(i,j,k);

      #if GEOMETRY == CYLINDRICAL
      fA(j,iFRPHI) *= fabs(x1[i]);
      #elif GEOMETRY == SPHERICAL
      IF (COVARIANT_MX2, fA(j,FR2) *= fabs(x1[i]));
      fA(j,FR3) *= fabs(x1[i]*sL[j]);
      #endif

    }}}
    #endif // GEOMETRY != CARTESIAN

    // J2. Compute radiation right hand side

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){

      long int offset1 = ni*(i + nj*k);
      long int offset  = NVAR*offset1;

      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary2D fA(&FA[offset],ni,NVAR);
      #endif

      // J2a. Radiation flux divergence

      #pragma acc loop seq
      NVAR_LOOP(nv) rhs(j,nv) = 0.0;

      #if GEOMETRY == CARTESIAN
      {
        double dtdx2 = dt/dx2[j];
        #pragma acc loop seq
        NRAD_LOOP(nv) rhs(j,nv) = -dtdx2*(flux(j+1,nv) - flux(j,nv));
      }
      #else
      {
        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NRAD_LOOP(nv) rhs(j,nv) = -dtdV*(fA(j+1,nv) - fA(j,nv));
      }
      #endif

      // J2a'. Irradiation
      #if IRRADIATION
      rhs(j, FIR) = 0.0;
      #endif

      // J2b. Geometric source terms (spherical only)

      #if GEOMETRY == SPHERICAL
      {
        Ary2D vc(&data->sweep->vC[offset],ni,NVAR);
        double v[NVAR];
        #pragma acc loop seq
        NVAR_LOOP(nv) v[nv] = vc(j,nv);

        double ct = grid->cot[j];
        double Pphph = EddTensor(v, iFRPHI, iFRPHI) * v[ENR];
        double Prth = EddTensor(v, iFRR, iFRTH) * v[ENR];
        double Sth = (ct*Pphph - Prth);
        #if COVARIANT_MX2 == YES
        double src_th = dt*Sth*run_config.rad_reduced_c;
        #else
        double src_th = dt*Sth*run_config.rad_reduced_c/x1[i];
        #endif
        rhs(j,iFRTH) += src_th;
      }
      #endif

    }}}

  // ------------------------------------------------------
  // K. Compute radiation rhs in the x3-direction
  // ------------------------------------------------------

  }else if (dir == KDIR){

    #if GEOMETRY != CARTESIAN

    // K1. Build area-weighted radiation flux

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend+1; k++){

      long int offset = NVAR*ni*(i + nj*j);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      #pragma acc loop seq
      NRAD_LOOP(nv) fA(k,nv) = flux(k,nv)*Ax3(i,j,k);

      #if GEOMETRY == CYLINDRICAL
      fA(k,iFRPHI) *= fabs(x1[i]);
      #elif GEOMETRY == SPHERICAL
      IF (COVARIANT_MX2, fA(k,FR2) *= fabs(x1[i]));
      fA(k,FR3) *= fabs(x1[i]*sC[j]);
      #endif

    }}}
    #endif // GEOMETRY != CARTESIAN

    // K2. Compute radiation right hand side

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend; k++){

      long int offset1 = ni*(i + nj*j);
      long int offset  = NVAR*offset1;

      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary2D fA(&FA[offset],ni,NVAR);
      #endif

      // K2a. Radiation flux divergence

      #pragma acc loop seq
      NVAR_LOOP(nv) rhs(k,nv) = 0.0;

      #if GEOMETRY == CARTESIAN
      {
        double dtdx3 = dt/dx3[k];
        #pragma acc loop seq
        NRAD_LOOP(nv) rhs(k,nv) = -dtdx3*(flux(k+1,nv) - flux(k,nv));
      }
      #else
      {
        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NRAD_LOOP(nv) rhs(k,nv) = -dtdV*(fA(k+1,nv) - fA(k,nv));
      }
      #endif

      // K2a'. Irradiation
      #if IRRADIATION
      rhs(k, FIR) = 0.0;
      #endif

    }}}

  }

}

#endif // RADIATION_NR