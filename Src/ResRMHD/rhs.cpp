////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Compute the right hand side of the conservative 
         resistive RMHD equations

  See also the \ref RHS_page.

  Here we have the following possible combinations to control the 
  magnetic and electric field divergences:

  <CENTER>
  div(B)|div(E)|NVAR|         q              | Included 
  ------|------|----|------------------------|------------------
   NO   |  NO  | 11 |    (En+ - En-)/dx      | RightHandSide()
   CT   |  NO  | 11 |    (En+ - En-)/dx      | RightHandSide()
   CT   |  CT  | 11 |(En[stg]+ - En[stg]-)/dx| CT_ComputeCharge(), CT_Update()
   GLM  |  NO  | 12 |    (En+ - En-)/dx      | RightHandSide()
   GLM  |  GLM | 13 | (F(phi)+ - F(phi)-)/dx | RightHandSide()
  </CENTER>

  \b References
   
  \author A. Mignone  (andrea.mignone@unito.it)
          G. Mattia   (mattia@mpia.de)
  \date   Oct 23, 2025

*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Compute components of tensor sources
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void TensorSource(State &state, RunConfig& run_config, double T[][3])
{
  double *v   = state.v;
  double wg, vmag2, inv_lor2, inv_lor;
  double Emag2, Bmag2, ptot;

  wg       = v[RHO]*Enthalpy(v,run_config);
  vmag2    = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  inv_lor2 = 1.0 - vmag2;
  inv_lor  = sqrt(inv_lor2);

  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  Emag2 = v[EX1]*v[EX1] + v[EX2]*v[EX2] + v[EX3]*v[EX3];

  ptot = v[PRS] + 0.5*(Bmag2 + Emag2);

  #if GEOMETRY == CYLINDRICAL
  T[JDIR][JDIR] = ptot + wg/inv_lor2*v[VX2]*v[VX2] - v[BX2]*v[BX2] - v[EX2]*v[EX2]; 
  #elif GEOMETRY == SPHERICAL
  T[JDIR][IDIR] =        wg/inv_lor2*v[VX2]*v[VX1] - v[BX2]*v[BX1] - v[EX2]*v[EX1];
  T[JDIR][JDIR] = ptot + wg/inv_lor2*v[VX2]*v[VX2] - v[BX2]*v[BX2] - v[EX2]*v[EX2];
  T[KDIR][KDIR] = ptot + wg/inv_lor2*v[VX3]*v[VX3] - v[BX3]*v[BX3] - v[EX3]*v[EX3];
  #else
  #endif
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
void RightHandSide (Data *data, timeStep *Dts, const Grid * const grid,
                    double dt, RBox *box, const int dir)
/*
  Solve the Riemann problem between L/R states using a
  Rusanov-Lax Friedrichs flux.
*/
////////////////////////////////////////////////////////////////////////
 {
 int  i, j, k, nv;
  int  ni, nj; 
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  double *x1  = grid->xC[IDIR], *x2  = grid->xC[JDIR], *x3  = grid->xC[KDIR];
  double *dx1 = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];
  double *x1m = grid->xL[IDIR];

  #if GEOMETRY == SPHERICAL
  double *sC = grid->sC, *sL  = grid->sL;
  double *dmu = grid->dmu;
  #endif

  int BXn = BX1 + dir;
  int EXn = EX1 + dir;

  #if GEOMETRY != CARTESIAN
  const Ary3D &Ax1 = grid->A[IDIR];
  const Ary3D &Ax2 = grid->A[JDIR];
  const Ary3D &Ax3 = grid->A[KDIR];
  #endif

  RunConfig &run_config = *(data->run_config);

  // 0. Allocate memory

  static double *FA;
  #if GEOMETRY != CARTESIAN
  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];
  if (FA == NULL) FA = ARRAY_1D(size_tot*NVAR, double);
  #endif

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("RHS",GREEN);

  if (dir == IDIR){

    // ----------------------------------------------------
    // I. Compute rhs in the x1-direction
    // ----------------------------------------------------
     
    #if GEOMETRY != CARTESIAN

    // ----------------------------------------------------
    // I1. Modify flux to account for geometry or gravity
    // ----------------------------------------------------
    #pragma omp target teams distribute parallel for collapse(3) map (to:data[:1], run_config[:1], grid[:1])
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){
      long int offset = NVAR*ni*(j + nj*k);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      
      // I1a. Multiply flux x area. 
      //      Avoid singularities in fluxes of controvariant quantities
      //      (like Bth, Bphi) by simplifying metric and area factors.
      //      See RMHD/rhs.cpp for comments

      #pragma acc loop seq
      NVAR_LOOP(nv) fA(i,nv) = flux(i,nv)*Ax1(i,j,k);

      #if GEOMETRY == CYLINDRICAL
      fA(i,MX2) *= x1m[i];
      double A = DIM_EXPAND(1.0, *dx2[j], *dx3[k]); 
      IF (COVARIANT_BX2, fA(i,BX2)  = flux(i,BX2)*A);
      IF (COVARIANT_EX2, fA(i,EX2)  = flux(i,EX2)*A);
      #endif
      #if GEOMETRY == SPHERICAL
      IF (COVARIANT_MX2, fA(i,MX2) *= x1m[i]);
      fA(i,MX3) *= x1m[i]*sC[j];

      double A = DIM_EXPAND(x1m[i], *dmu[j], *dx3[k]);
      IF (COVARIANT_BX2, fA(i,BX2) = flux(i,BX2)*A); 
      IF (COVARIANT_EX2, fA(i,EX2) = flux(i,EX2)*A);
      fA(i,BX3)  = flux(i,BX3)*A/sC[j];  
      fA(i,EX3)  = flux(i,EX3)*A/sC[j];
      #endif

      }
    }}
    #endif

    // ----------------------------------------------------
    // I2. Build rhs in the x1-dir
    // ----------------------------------------------------
    
    #pragma omp target teams distribute parallel for collapse(3) map (to:data[:1], run_config[:1], grid[:1])
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){

      long int offset1 = ni*(j + nj*k);
      long int offset  = NVAR*offset1;

      Ary1D Bn   (&data->sweep->Bn[offset1], ni);
      Ary1D cmax (&Dts->cmax[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);

      Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset],ni,NVAR);
      Ary1D En(&data->sweep->En[offset1],ni);

      double dtdx1 = dt/dx1[i];
      double v[NVAR];
      double src[NVAR] = {0.0};  // Aggregate initialization
      
      State state;
      STATE_SET (state, v, NULL, NULL, x1[i], x2[j], x3[k]);

      #if GEOMETRY != CARTESIAN
      double T[3][3];
      Metric metric;
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(i,nv);

      // I2a. Conservative flux difference

      #if GEOMETRY == CARTESIAN

      #pragma acc loop seq
      NVAR_LOOP(nv)  rhs(i,nv) = -dtdx1*(flux(i+1,nv) - flux(i,nv));

      #else

      double dtdV = dt/grid->dV(i,j,k);

      #pragma acc loop seq
      NVAR_LOOP(nv) rhs(i,nv) = -dtdV*(fA(i+1,nv) - fA(i,nv));

      #endif

      // I2b. Charge source terms

      #ifndef HIGH_ORDER // HO Source term computed later in Sources() 
        double qdt;

        #if DIVE_CONTROL == DIV_CLEANING
        qdt             = dt*vc(i,CRG); 
        rhs(i,PHI_GLM) += qdt;
        #elif DIVE_CONTROL == CONSTRAINED_TRANSPORT
        #if GEOMETRY == CARTESIAN
        qdt = dtdx1*(En(i+1) - En(i));
        #else
        qdt = dtdV*(Ax1(i+1,j,k)*En(i+1) - Ax1(i,j,k)*En(i));
        #endif
        #else
        #if GEOMETRY == CARTESIAN
        qdt = dtdx1*(vL(i+1,EXn) - vR(i,EXn));
        #else
        qdt = dtdV*(Ax1(i+1,j,k)*vL(i+1,EXn) - Ax1(i,j,k)*vR(i,EXn));
        #endif
        // qdt = dtdx1*(vc(i+1,EXn) - vc(i-1,EXn))*0.5;
        #endif // DIVE_CONTROL == DIV_CLEANING

        rhs(i,EX1) -= qdt*vc(i,VX1);
        rhs(i,EX2) -= qdt*vc(i,VX2);
        rhs(i,EX3) -= qdt*vc(i,VX3);

      #endif

      // I2c. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      src[MX1]           += dt*v[RHO]*g[IDIR];
      IF_ENERGY(src[ENG] += dt*0.5*(flux(i+1,RHO) + flux(i,RHO))*g[IDIR];)
      #endif

      // I2d. Add geometrical source terms

      #if GEOMETRY == CYLINDRICAL  // ** CYLINDRICAL GEOMETRY **

        double r = x1[i];
        TensorSource (state, run_config, T);
        src[MX1] += dt*T[JDIR][JDIR]/r;

        #if DIVB_CONTROL == DIV_CLEANING
        src[BX1] += dt*v[PSI_GLM]/r;
        #endif
        #if DIVE_CONTROL == DIV_CLEANING
        src[EX1] += dt*v[PHI_GLM]/r;
        #endif

        #if !COVARIANT_BX2
        src[BX2] -= dt*v[EX3]/r;
        #endif

        #if !COVARIANT_EX2
        src[EX2] += dt*v[BX3]/r;
        #endif

      #elif GEOMETRY == SPHERICAL

        double r = grid->rt[i]; // x1[i];
        TensorSource (state, run_config, T);
        src[MX1] += dt*(T[JDIR][JDIR] + T[KDIR][KDIR])/r;
        #if !INCLUDE_JDIR
        src[MX1] -= dt*T[JDIR][IDIR]*grid->cot[j]/r;
        #endif
        #if COVARIANT_MX2 == YES
        src[MX2] += dt*T[KDIR][KDIR]*grid->cot[j];
        #if !INCLUDE_JDIR
        src[MX2] -= dt*T[JDIR][JDIR]*grid->cot[j];
        #endif
        #else
        src[MX2] += dt*(T[KDIR][KDIR]*grid->cot[j] - T[JDIR][IDIR])/x1[i];
        #if !INCLUDE_JDIR
        src[MX2] -= dt*T[JDIR][JDIR]*grid->cot[j]/x1[i];
        #endif
        #endif

        #if DIVB_CONTROL == DIV_CLEANING
        src[BX1] += dt*2.0*v[PSI_GLM]/r;
        #if INCLUDE_JDIR
        #if COVARIANT_BX2 == YES
        src[BX2] += dt*v[PSI_GLM]*grid->cot[j]/(r*x1[i]);
        #else
        src[BX2] += dt*v[PSI_GLM]*grid->cot[j]/r;
        #endif
        #endif
        #endif
        #if DIVE_CONTROL == DIV_CLEANING
        src[EX1] += dt*2.0*v[PHI_GLM]/r;
        #if INCLUDE_JDIR
        #if COVARIANT_EX2 == YES
        src[EX2] += dt*v[PHI_GLM]*grid->cot[j]/(r*x1[i]);
        #else
        src[EX2] += dt*v[PHI_GLM]*grid->cot[j]/r;
        #endif
        #endif
        #endif

        #if !COVARIANT_BX2
        src[BX2] -= dt*v[EX3]/r;
        #endif

        #if !COVARIANT_EX2
        src[EX2] += dt*v[BX3]/r;
        #endif

      #endif

      // I2f. Transfrom to conservative and add source terms

      #if GEOMETRY != CARTESIAN
      OrthoToCons(src, metric);

      #pragma acc loop seq 
      NFLX_LOOP(nv) rhs(i,nv) += src[nv];   
      #endif

      // I2g. Compute inverse time step 

      data->C_dt(i,j,k) = 0.5*(cmax(i) + cmax(i+1))/dx1[i];

      }  // End loop on i  
    }}  // End loop on k,j 

  }else if (dir == JDIR){

    // ----------------------------------------------------
    // J. Compute rhs in the x2-direction
    // ----------------------------------------------------

    #if GEOMETRY != CARTESIAN

    // ----------------------------------------------------
    // J1. Modify flux to account for geometry or gravity
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) map (to:data[:1], run_config[:1], grid[:1])
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend+1; j++){
      long int offset = NVAR*ni*(i + nj*k);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      // J1a. Multiply flux x area
      //      Avoid singularities in fluxes of controvariant quantities 
      //      (like Bth, Bphi) by simplifying metric and area factors.

      #pragma acc loop seq
      NVAR_LOOP(nv) fA(j,nv) = flux(j,nv)*Ax2(i,j,k);

      #if GEOMETRY == CYLINDRICAL
      fA(j,MX2) *= x1[i];
      IF (COVARIANT_BX2, fA(j,BX2) /= x1[i];) 
      IF (COVARIANT_EX2, fA(j,EX2) /= x1[i];)
      #endif
      #if GEOMETRY == SPHERICAL
      IF (COVARIANT_MX2, fA(j,MX2) *= x1[i]);
      fA(j,MX3) *= x1[i]*sL[j];
      double A = DIM_EXPAND(dx1[i],  , *dx3[k]);
      IF (COVARIANT_BX2, fA(j,BX2) = flux(j,BX2)*A*sL[j]);
      IF (COVARIANT_EX2, fA(j,EX2) = flux(j,EX2)*A*sL[j]);
      fA(j,BX3) = flux(j,BX3)*A; 
      fA(j,EX3) = flux(j,EX3)*A; 
      #endif

      }
    }}

    #endif

    // ----------------------------------------------------
    // J2. Build rhs in the x2-dir
    // ----------------------------------------------------
    #pragma omp target teams distribute parallel for collapse(3) map (to:data[:1], grid[:1], Dts[:1])
    #pragma acc parallel loop collapse(3) present(data, grid, Dts)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){
      long int offset1 = ni*(i + nj*k);
      long int offset  = NVAR*offset1;

      Ary1D Bn   (&data->sweep->Bn[offset1],ni);
      Ary1D cmax(&Dts->cmax[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);

      Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset],ni,NVAR);
      Ary1D En(&data->sweep->En[offset1],ni);

      double dtdx2 = dt/dx2[j];
      double v[NVAR];
      double src[NVAR] = {0.0};  // Aggregate initialization

      State state;
      STATE_SET (state, v, NULL, NULL, x1[i], x2[j], x3[k]);

      Metric metric;
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(j,nv);

      // J2a. Conservative flux difference

      #if GEOMETRY == CARTESIAN

      #pragma acc loop seq
      NVAR_LOOP(nv) rhs(j,nv) = -dtdx2*(flux(j+1,nv) - flux(j,nv));

      #else

      double dtdV = dt/grid->dV(i,j,k);
      #pragma acc loop seq
      NVAR_LOOP(nv) rhs(j,nv) = -dtdV*(fA(j+1,nv) - fA(j,nv));

      #endif

      // J2b. Charge source terms

      #ifndef HIGH_ORDER  // HO Source term computed later in Sources() 
        double qdt;

        #if DIVE_CONTROL == DIV_CLEANING
        qdt             = dt*vc(j,CRG);   
//          qdt = dtdx1*(flux(j+1,PHI_GLM) - flux(j,PHI_GLM));  
        rhs(j,PHI_GLM) += qdt;
        #elif DIVE_CONTROL == CONSTRAINED_TRANSPORT
        #if GEOMETRY == CARTESIAN
        qdt = dtdx2*(En(j+1) - En(j));
        #else
        qdt = dtdV*(Ax2(i,j+1,k)*En(j+1) - Ax2(i,j,k)*En(j));
        #endif
        #else
        #if GEOMETRY == CARTESIAN
        qdt = dtdx2*(vL(j+1,EXn) - vR(j,EXn));
        #else
        qdt = dtdV*(Ax2(i,j+1,k)*vL(j+1,EXn) - Ax2(i,j,k)*vR(j,EXn));
        #endif
        #endif // DIVE_CONTROL == DIV_CLEANING

        rhs(j,EX1) -= qdt*vc(j,VX1);
        rhs(j,EX2) -= qdt*vc(j,VX2);
        rhs(j,EX3) -= qdt*vc(j,VX3);

      #endif

      // J2c. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      src[MX2] += dt*v[RHO]*g[JDIR];
      IF_ENERGY(src[ENG] += dt*0.5*(flux(j+1,RHO) + flux(j,RHO))*g[JDIR];)
      #endif 

      // J2d. Add geometrical source terms

      #if GEOMETRY == SPHERICAL
      #endif

      // J2f. Transfrom to conservative and add source terms

      #if GEOMETRY != CARTESIAN
      OrthoToCons(src, metric);

      #pragma acc loop seq 
      NFLX_LOOP(nv) rhs(j,nv) += src[nv];   
      #endif

      // J2g. Compute inverse time step 

      double dl = dx2[j]*metric.h2;
      data->C_dt(i,j,k) += 0.5*(cmax(j) + cmax(j+1))/dl;

      }  // End loop on j  
    }}

  }else if (dir == KDIR){

    // ----------------------------------------------------
    // K. Compute rhs in the x3-direction
    // ----------------------------------------------------

    #if GEOMETRY != CARTESIAN

    // ----------------------------------------------------
    // K1. Modify flux to account for geometry
    // ----------------------------------------------------
    
    #pragma omp target teams distribute parallel for collapse(3) map (to:data[:1], run_config[:1], grid[:1])
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend+1; k++)
    {
      long int offset = NVAR*ni*(i + nj*j);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      
      // K1a. Multiply flux x area.
      //      Avoid singularities in fluxes of controvariant quantities 
      //      (like Bth, Bphi) by simplifying metric and area factors.

      #pragma acc loop seq
      NVAR_LOOP(nv) fA(k,nv) = flux(k,nv)*Ax3(i,j,k);

      #if GEOMETRY == CYLINDRICAL
      fA(k,MX2) *= x1[i];
      double A  = DIM_EXPAND(dx1[i], *dx2[j], *1.0);
      IF (COVARIANT_BX2, fA(k,BX2) = A*flux(k,BX2));  // r*dr*dphi*(flux/r) = dr*dphi*flux
      IF (COVARIANT_EX2, fA(k,EX2) = A*flux(k,EX2));
      #endif
      #if GEOMETRY == SPHERICAL
      IF (COVARIANT_MX2, fA(k,MX2) *= x1[i]);
      fA(k,MX3) *= x1[i]*sC[j];
      double A  = DIM_EXPAND(dx1[i], *dx2[j], *1.0);
      IF (COVARIANT_BX2, fA(k,BX2) = flux(k,BX2)*A;)
      IF (COVARIANT_EX2, fA(k,EX2) = flux(k,EX2)*A;)
      fA(k,BX3) = flux(k,BX3)*A/sC[j]; 
      fA(k,EX3) = flux(k,EX3)*A/sC[j];
      #endif
      }
    }}

    #endif

    // ----------------------------------------------------
    //  K2. Build rhs in the x3-dir
    // ----------------------------------------------------
    
    #pragma omp target teams distribute parallel for collapse(3) map (to:data[:1], grid[:1], Dts[:1])
    #pragma acc parallel loop collapse(3) present(data, grid, Dts)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend; k++){
      long int offset = NVAR*ni*(i + nj*j);
      long int offset1 = ni*(i + nj*j);

      Ary1D Bn   (&data->sweep->Bn[offset1],ni);
      Ary1D cmax(&Dts->cmax[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);

      //Ary2D v(&data->sweep->vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset],ni,NVAR);
      Ary1D En(&data->sweep->En[offset1],ni);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);

      double dtdx3 = dt/dx3[k];
      double v[NVAR];
      double src[NVAR] = {0.0};  // Aggregate initialization

      State state;
      STATE_SET (state, v, NULL, NULL, x1[i], x2[j], x3[k]);

      Metric metric;
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif
      
      //#pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(k,nv);

      // K2a. Conservative flux difference
      
      #if GEOMETRY == CARTESIAN

      NVAR_LOOP(nv) rhs(k,nv) = -dtdx3*(flux(k+1,nv) - flux(k,nv));

      #else

      double dtdV = dt/grid->dV(i,j,k);
      //#pragma acc loop seq
      NVAR_LOOP(nv) rhs(k,nv) = -dtdV*(fA(k+1,nv) - fA(k,nv));

      #if DIVB_CONTROL == DIV_CLEANING
      #if GEOMETRY == SPHERICAL
      double inv_rgrad = dx2[j]/(grid->rt[i]*grid->dmu[j]);
      rhs(k,BXn) = -dtdx3*inv_rgrad*(flux(k+1,BXn) - flux(k,BXn));
      #else
      rhs(k,BXn) = -dtdx3*(flux(k+1,BXn) - flux(k,BXn));
      #endif
      #endif

      #if DIVE_CONTROL == DIV_CLEANING
      #if GEOMETRY == SPHERICAL
      double inv_rgrad_e = dx2[j]/(grid->rt[i]*grid->dmu[j]);
      rhs(k,EXn) = -dtdx3*inv_rgrad_e*(flux(k+1,EXn) - flux(k,EXn));
      #else
      rhs(k,EXn) = -dtdx3*(flux(k+1,EXn) - flux(k,EXn));
      #endif
      #endif
      #endif

      // K2b. Charge source terms

      #ifndef HIGH_ORDER  // HO Source term computed later in Sources() 
        double qdt;

        #if DIVE_CONTROL == DIV_CLEANING
        qdt             = dt*vc(k,CRG);
        // qdt = dtdx3*(flux(k+1,PHI_GLM) - flux(k,PHI_GLM)); 
        rhs(k,PHI_GLM) += qdt;
        #elif DIVE_CONTROL == CONSTRAINED_TRANSPORT
        #if GEOMETRY == CARTESIAN
        qdt = dtdx3*(En(k+1) - En(k));
        #else
        qdt = dtdV*(Ax3(i,j,k+1)*En(k+1) - Ax3(i,j,k)*En(k));
        #endif
        #else
        #if GEOMETRY == CARTESIAN
        qdt = dtdx3*(vL(k+1,EXn) - vR(k,EXn));
        #else
        qdt = dtdV*(Ax3(i,j,k+1)*vL(k+1,EXn) - Ax3(i,j,k)*vR(k,EXn));
        #endif
        #endif // DIVE_CONTROL == DIV_CLEANING

        rhs(k,EX1) -= qdt*vc(k,VX1);
        rhs(k,EX2) -= qdt*vc(k,VX2);
        rhs(k,EX3) -= qdt*vc(k,VX3);

      #endif

      // K2c. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      src[MX3] += dt*v[RHO]*g[KDIR];
      IF_ENERGY(src[ENG] += dt*0.5*(flux(k+1,RHO) + flux(k,RHO))*g[KDIR];)
      #endif  // (BODY_FORCE & VECTOR) 

      // K2f. Transfrom to conservative and add source terms

      #if GEOMETRY != CARTESIAN
      OrthoToCons(src, metric);

      //#pragma acc loop seq 
      NFLX_LOOP(nv) rhs(k,nv) += src[nv];   
      #endif

      // K2g. Compute inverse time step 

      double dl = dx3[k]*metric.h3;
      data->C_dt(i,j,k) += 0.5*(cmax(k) + cmax(k+1))/dl;

      }  // End loop on k  
    }}
  }
  mynvtxstop_();
}
