////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Compute the right hand side of the conservative 
         RHD/RMHD equations.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Dec 18, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// -- Set macro for include/exclude div(B) source terms -- 

#if (PHYSICS == RMHD) && (DIVB_CONTROL == EIGHT_WAVES)
  #define INCLUDE_DIVB_SOURCE YES 
#elif (PHYSICS == RMHD) && (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
  #define INCLUDE_DIVB_SOURCE YES 
#else
  #define INCLUDE_DIVB_SOURCE NO
#endif

#if (PHYSICS == RMHD) && (DIVB_CONTROL == EIGHT_WAVES)
////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void PowellSourceTerm(const double divB, const double* v, double* src) 
{
#if EIGHT_WAVES_REDUCED == NO

  double g_2 = 1.0 - (v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]);
  double vB  = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];

  src[MX1] = divB * (v[BX1]*g_2 + vB*v[VX1]);
  src[MX2] = divB * (v[BX2]*g_2 + vB*v[VX2]);
  src[MX3] = divB * (v[BX3]*g_2 + vB*v[VX3]);
  src[BX1] = divB * v[VX1];
  src[BX2] = divB * v[VX2];
  src[BX3] = divB * v[VX3];

  #if HAVE_ENERGY
  src[ENG] = divB * vB;
  #endif

#elif EIGHT_WAVES_REDUCED == YES

  src[BX1] = divB * v[VX1];
  src[BX2] = divB * v[VX2];
  src[BX3] = divB * v[VX3];

#endif
}
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Compute components of tensor sources
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void TensorSource(State &state, RunConfig& run_config, double T[][3])
{
  double *v   = state.v;
  double wg, wt, vmag2, inv_lor2, inv_lor;

  wg       = v[RHO]*Enthalpy(v,run_config);
  vmag2    = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  inv_lor2 = 1.0 - vmag2;
  inv_lor  = sqrt(inv_lor2);

  #if PHYSICS == RMHD
  double vB, Bmag2, ptot;
  vB    = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];

  wt   = wg + Bmag2*inv_lor2 + vB*vB;
  ptot = v[PRS] + 0.5*(Bmag2*inv_lor2 + vB*vB);
  #endif

  //  Compute T22 component in cylindrical coordinates 

  #if GEOMETRY == CYLINDRICAL

  #if PHYSICS == RMHD
  double bx2  = v[BX2]*inv_lor + vB*v[VX2]/inv_lor;
  T[JDIR][JDIR] = ptot   + wt/inv_lor2*v[VX2]*v[VX2] - bx2*bx2; 
  #elif PHYSICS == RHD
  T[JDIR][JDIR] = v[PRS] + wg/inv_lor2*v[VX2]*v[VX2];
  #endif       

  // Compute T21, T22 and T33 components in spherical coordinates 

  #elif GEOMETRY == SPHERICAL

  #if PHYSICS == RMHD
  double bx1    = v[BX1]*inv_lor + vB*v[VX1]/inv_lor;
  double bx2    = v[BX2]*inv_lor + vB*v[VX2]/inv_lor;
  double bx3    = v[BX3]*inv_lor + vB*v[VX3]/inv_lor;
  T[JDIR][IDIR] =        wt/inv_lor2*v[VX2]*v[VX1] - bx2*bx1; 
  T[JDIR][JDIR] = ptot + wt/inv_lor2*v[VX2]*v[VX2] - bx2*bx2; 
  T[KDIR][KDIR] = ptot + wt/inv_lor2*v[VX3]*v[VX3] - bx3*bx3; 
  #elif PHYSICS == RHD
  T[JDIR][IDIR] =          wg/inv_lor2*v[VX2]*v[VX1];
  T[JDIR][JDIR] = v[PRS] + wg/inv_lor2*v[VX2]*v[VX2];
  T[KDIR][KDIR] = v[PRS] + wg/inv_lor2*v[VX3]*v[VX3];
  #endif       

  #else
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Compute right hand side for RHD / RMHD equations.
  
  \param [in,out]  data  Pointer to Data structure 
  \param [in]      Dts   Pointer to timeStep structure 
  \param [in]      grid  Pointer to Grid structure
  \param [in]      dt    time increment
  \param [in]      box   RBox on which rhs has to be computed
*/
////////////////////////////////////////////////////////////////////////

#define RHS_VALUE(index,n) (fused ? local_rhs[n] : rhs(index,n))
template<bool fused>
static void RightHandSideImpl(Data *data, timeStep *Dts, const Grid *const grid,
                             double dt, RBox *box, int dir, const Ary4D &target)
 {
  int i, j, k, nv;
  int ni, nj;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  double *x1  = grid->xC[IDIR], *x2  = grid->xC[JDIR], *x3  = grid->xC[KDIR];
  double *dx1 = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];
  #if GEOMETRY != CARTESIAN
  double *x1m = grid->xL[IDIR];
  double A;
  const Ary3D &Ax1 = grid->A[IDIR];
  const Ary3D &Ax2 = grid->A[JDIR];
  const Ary3D &Ax3 = grid->A[KDIR];
  #endif

  #if GEOMETRY == SPHERICAL
  double *sC = grid->sC, *sL  = grid->sL;
  double *dmu = grid->dmu;
  #endif

  #if PHYSICS == RMHD
  double divB;
  double gradPSI;
  int    BXn = BX1 + dir;
  #endif

  RunConfig &run_config = *(data->run_config);

  // 0. Allocate memory

  static double *FA;

  #if GEOMETRY != CARTESIAN
  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];
  if (FA == NULL) {
    FA = ARRAY_1D(size_tot*NVAR, double);
  }
  #endif

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("RHS",GREEN);

  // ------------------------------------------------------
  // I. Compute rhs in the x1-direction
  // ------------------------------------------------------
     
  if (dir == IDIR){

    #if GEOMETRY != CARTESIAN

    // ----------------------------------------------------
    // I1. Modify flux to account for geometry or gravity
    // ----------------------------------------------------
    
    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){

      long int offset = NVAR*ni*(j + nj*k);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      // I1a. Multiply flux x area. 
      //      Avoid singularities in fluxes of controvariant quantities
      //      (like Bth, Bphi) by simplifying metric and area factors.

      #if GEOMETRY != CARTESIAN

        #pragma acc loop seq
        NVAR_LOOP(nv) fA(i,nv) = flux(i,nv)*Ax1(i,j,k);

        #if GEOMETRY == CYLINDRICAL
        fA(i,MX2) *= x1m[i];
        #if RADIATION
        fA(i,FR2) *= x1m[i];
        #endif
        #if PHYSICS == RMHD
        // Simplify area as follows:
        // A*f(BX2) = A*[flux/r] = r*dphi*dz*[flux/r] = dphi*dz*flux

        IF (COVARIANT_BX2, A = DIM_EXPAND(1.0, *dx2[j], *dx3[k]); 
                           fA(i,BX2) = flux(i,BX2)*A;)  
        #endif
        #endif  // GEOMETRY == CYLINDRICAL

        #if GEOMETRY == SPHERICAL
        IF (COVARIANT_MX2, fA(i,MX2) *= x1m[i]);
        fA(i,MX3) *= x1m[i]*sC[j];
        #if RADIATION
        IF (COVARIANT_MX2, fA(i,FR2) *= x1m[i]);
        fA(i,FR3) *= x1m[i]*sC[j];
        #endif
        #if PHYSICS == RMHD
        // Simplify area as follows:  
        // A*F(Bth)  = A*[flux/r] = r^2*dmu*dphi*[flux/r] = r*dmu*dphi*flux        
        // A*F(Bphi) = A*[flux/(r*s)] = r^2*dmu*dphi*[flux/(r*s)] 
        //                            = r*dmu*dphi*flux/s

        A   = DIM_EXPAND(x1m[i], *dmu[j], *dx3[k]);
        IF (COVARIANT_BX2, fA(i,BX2) = flux(i,BX2)*A); 
        fA(i,BX3)  = flux(i,BX3)*A/sC[j];  
        #endif
        #endif // GEOMETRY == SPHERICAL

      #endif // GEOMETRY != CARTESIAN 
    }}}
    #endif

    // ----------------------------------------------------
    // I2. Build rhs in the x1-dir
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, run_config, grid, target)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      long int offset1 = ni*(j + nj*k);
      long int offset  = NVAR*offset1;
      
      Ary1D Bn   (&data->sweep->Bn[offset1],   ni);
      Ary1D PSIn(&data->sweep->PSIn[offset1],ni);
      Ary1D cmax (&Dts->cmax[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);
 
      double v[NVAR];
      double local_rhs[NVAR];
      #if GEOMETRY != CARTESIAN
      double T[3][3];
      Metric metric;
      #endif
      State state;

      STATE_SET (state, v, NULL, NULL, x1[i], x2[j], x3[k]);
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(i,nv);

      // I2a. Conservative flux difference

      #if GEOMETRY == CARTESIAN         
        double dtdx1 = dt/dx1[i];
        #pragma acc loop seq
        NVAR_LOOP(nv)  RHS_VALUE(i,nv) = -dtdx1*(flux(i+1,nv) - flux(i,nv));

      #else

        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NVAR_LOOP(nv) RHS_VALUE(i,nv) = -dtdV*(fA(i+1,nv) - fA(i,nv));

      #endif

      // I2b. Add gravity

      #if (BODY_FORCE & VECTOR)
        double g[3];
        BodyForceVector(state, run_config, g);

        RHS_VALUE(i,MX1)   += dt*v[RHO]*g[IDIR];
        IF_ENERGY (RHS_VALUE(i,ENG) += dt*0.5*(flux(i+1,RHO) + flux(i,RHO))*g[IDIR];)
      #endif

      // I2c. Add geometrical source terms

      #if GEOMETRY == CYLINDRICAL  // ** CYLINDRICAL GEOMETRY ** 
      {
        double r     = x1[i];
        TensorSource (state, run_config, T);
        RHS_VALUE(i,MX1) += dt*T[JDIR][JDIR]/r;
        #if RADIATION
        {
          double Sr = EddTensor(v, FR2, FR2)*v[ENR];
          RHS_VALUE(i,FR1) += dt*Sr/r;
        }
        #endif
        #ifdef GLM_MHD
        RHS_VALUE(i,BX1) += dt*v[PSI_GLM]/r;
        #endif
        #if PHYSICS == RMHD && (!COVARIANT_BX2)
        {
          double Ez = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]); 
          RHS_VALUE(i,BX2) -= dt*Ez/r;
        }
        #endif
      }
      #elif GEOMETRY == SPHERICAL  // ** SPHERICAL GEOMETRY **
      {
        double r     = grid->rt[i]; // x1[i];
        TensorSource (state, run_config, T);
        RHS_VALUE(i,MX1) += dt*(T[JDIR][JDIR] + T[KDIR][KDIR])/r;
        #if !INCLUDE_JDIR
        RHS_VALUE(i,MX1) -= dt*T[JDIR][IDIR]*grid->cot[j]/r;
        #endif

        #if COVARIANT_MX2 == YES
        RHS_VALUE(i,MX2) += dt*T[KDIR][KDIR]*grid->cot[j];
        #if !INCLUDE_JDIR
        RHS_VALUE(i,MX2) -= dt*T[JDIR][JDIR]*grid->cot[j];
        #endif
        #else
        RHS_VALUE(i,MX2) += dt*(T[KDIR][KDIR]*grid->cot[j] - T[JDIR][IDIR])/x1[i];
        #if !INCLUDE_JDIR
        RHS_VALUE(i,MX2) -= dt*T[JDIR][JDIR]*grid->cot[j]/x1[i];
        #endif
        #endif

        #if RADIATION
        {
          double Pthth = EddTensor(v, iFRTH,  iFRTH)  * v[ENR];
          double Pphph = EddTensor(v, iFRPHI, iFRPHI) * v[ENR];
          RHS_VALUE(i,iFRR) += dt*(Pthth + Pphph)/r;
        }
        #endif

        #ifdef GLM_MHD
          RHS_VALUE(i,BX1) += dt*2.0*v[PSI_GLM]/r;
          #if INCLUDE_JDIR
          #if COVARIANT_BX2 
          RHS_VALUE(i,BX2) += dt*v[PSI_GLM]*grid->cot[j]/(r*x1[i]);
          #else
          RHS_VALUE(i,BX2) += dt*v[PSI_GLM]*grid->cot[j]/r;
          #endif
          #endif
        #endif

        #if PHYSICS == RMHD && (!COVARIANT_BX2)
        {
          double Ez = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]); 
          RHS_VALUE(i,BX2) -= dt*Ez/r;
        }
        #endif

      }
      #endif // GEOMETRY == SPHERICAL 

      // I2e. Add div(B) source terms

      #if INCLUDE_DIVB_SOURCE
        double src[NVAR] = {0.0};  // Aggregate initialization
        
        #if GEOMETRY == CARTESIAN
        divB = (Bn(i+1) - Bn(i))/dx1[i];
        #else
        divB = (Ax1(i+1,j,k)*Bn(i+1) - Ax1(i,j,k)*Bn(i))/grid->dV(i,j,k);
        #endif
        #if DIVB_CONTROL == EIGHT_WAVES 
        PowellSourceTerm(dt*divB, v, src);
        #elif (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
        gradPSI = (PSIn(i+1) - PSIn(i))/dx1[i];
        GLM_ExtendedSource (dt*divB, dt*gradPSI, v, src, dir);
        #endif
        #if GEOMETRY != CARTESIAN
        OrthoToCons(src, metric);
        #endif
        
        #pragma acc loop seq 
        NFLX_LOOP(nv) RHS_VALUE(i,nv) -= src[nv];
      #endif // INCLUDE_DIVB_SOURCE

      // I2h. Compute inverse time step 

      data->C_dt(i,j,k) = 0.5*(cmax(i) + cmax(i+1))/dx1[i];
      if constexpr (fused) {
        #pragma acc loop seq
        NFLUID_LOOP(nv) target(i,j,k,nv) = local_rhs[nv];
      }
    }}}

  }else if (dir == JDIR){

    // ----------------------------------------------------
    // J. Compute rhs in the x2-direction
    // ----------------------------------------------------

    #if GEOMETRY != CARTESIAN

    // ----------------------------------------------------
    // J1. Modify flux to account for geometry or gravity
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, grid)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend+1; j++){

      long int offset = NVAR*ni*(i + nj*k);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      // J1a. Multiply flux x area
      //      Avoid singularities in fluxes of controvariant quantities 
      //      (like Bth, Bphi) by simplifying metric and area factors.
      
      #if GEOMETRY != CARTESIAN

        #pragma acc loop seq
        NVAR_LOOP(nv) fA(j,nv) = flux(j,nv)*Ax2(i,j,k);

        #if GEOMETRY == CYLINDRICAL
        fA(j,MX2) *= x1[i];
        #if RADIATION
        fA(j,FR2) *= x1[i];
        #endif
        #if PHYSICS == RMHD
        IF (COVARIANT_BX2, fA(j,BX2) /= x1[i];) // A*F(BX2) = A*[flux/r]  
                                                //          = dr*dz*[flux/r] 
        #endif
        #endif  // GEOMETRY == CYLINDRICAL

        #if GEOMETRY == SPHERICAL
        IF (COVARIANT_MX2, fA(j,MX2) *= x1[i]);
        fA(j,MX3) *= x1[i]*sL[j];
        #if RADIATION
        IF (COVARIANT_MX2, fA(j,FR2) *= x1[i]);
        fA(j,FR3) *= x1[i]*sL[j];
        #endif
        #if PHYSICS == RMHD
        // Simplify area as follows:  
        // A*F(Bth)  = A*[flux/r]     = r*dr*s*dphi*flux/r     = s*dr*dphi*flux   
        // A*F(Bphi) = A*[flux/(r*s)] = r*dr*s*dphi*flux/(r*s) =   dr*dphi*flux          

        A = DIM_EXPAND(dx1[i],  , *dx3[k]);
        IF (COVARIANT_BX2, fA(j,BX2) = flux(j,BX2)*A*sL[j]);
        fA(j,BX3) = flux(j,BX3)*A; 
        #endif
        #endif  // GEOMETRY == SPHERICAL

      #endif 

    }}}
    #endif

    // ----------------------------------------------------
    // J2. Build rhs in the x2-dir
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid, target)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){      

      long int offset1 = ni*(i + nj*k);
      long int offset  = NVAR*offset1;
      
      Ary1D Bn   (&data->sweep->Bn[offset1],ni);
      Ary1D PSIn (&data->sweep->PSIn[offset1],ni);
      Ary1D cmax (&Dts->cmax[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);

      double v[NVAR];
      double local_rhs[NVAR];
      State state;
      Metric metric;
      
      STATE_SET (state, v, NULL, NULL, x1[i], x2[j], x3[k]);
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      double dl   = dx2[j]*metric.h2;

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(j,nv);

      // J2a. Conservative flux difference

      #if GEOMETRY == CARTESIAN
        double dtdx2 = dt/dx2[j];
        #pragma acc loop seq
        NVAR_LOOP(nv) RHS_VALUE(j,nv) = -dtdx2*(flux(j+1,nv) - flux(j,nv));

      #else

        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NVAR_LOOP(nv) RHS_VALUE(j,nv) = -dtdV*(fA(j+1,nv) - fA(j,nv));

      #endif

      // J2b. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      RHS_VALUE(j,MX2) += dt*v[RHO]*g[JDIR];
      IF_ENERGY(RHS_VALUE(j,ENG) += dt*0.5*(flux(j+1,RHO) + flux(j,RHO))*g[JDIR];)
      #endif 

      // J2c. Add geometrical source terms

      #if RADIATION && (GEOMETRY == SPHERICAL)
      {
        double ct    = grid->cot[j];
        double Pphph = EddTensor(v, iFRPHI, iFRPHI) * v[ENR];

        #if COVARIANT_MX2 == YES
        RHS_VALUE(j,iFRTH) += dt*Pphph*ct;
        #else
        double Prth = EddTensor(v, iFRR, iFRTH) * v[ENR];
        RHS_VALUE(j,iFRTH) += dt*(Pphph*ct - Prth)/x1[i];
        #endif
      }
      #endif

      // J2e. Add div(B) source terms

      #if INCLUDE_DIVB_SOURCE
        double src[NVAR] = {0.0};  // Aggregate initialization
        
        #if GEOMETRY == CARTESIAN
        divB = (Bn(j+1) - Bn(j))/dx2[j];
        #else
        divB = (Ax2(i,j+1,k)*Bn(j+1) - Ax2(i,j,k)*Bn(j))/grid->dV(i,j,k);
        #endif
        #if DIVB_CONTROL == EIGHT_WAVES 
        PowellSourceTerm(dt*divB, v, src);
        #elif (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
        gradPSI = (PSIn(j+1) - PSIn(j))/dl;
        GLM_ExtendedSource (dt*divB, dt*gradPSI, v, src);
        #endif
        #if GEOMETRY != CARTESIAN
        OrthoToCons(src, metric);
        #endif
        #pragma acc loop seq 
        NFLX_LOOP(nv) RHS_VALUE(j,nv) -= src[nv];
      #endif // INCLUDE_DIVB_SOURCE

      // J2f. Compute inverse time step 

      double q = 1.0;
      #if (RING_AVERAGE > 1) && (GEOMETRY == CYLINDRICAL)
      q = 1.0/grid->ring_av_csize[i];
      #endif
      data->C_dt(i,j,k) += q*0.5*(cmax(j) + cmax(j+1))/dl;
      if constexpr (fused) {
        #pragma acc loop seq
        NFLUID_LOOP(nv) target(i,j,k,nv) += local_rhs[nv];
      }
    }}}

  }else if (dir == KDIR){

    // ----------------------------------------------------
    // K. Compute rhs in the x3-direction
    // ----------------------------------------------------

    #if GEOMETRY != CARTESIAN

    // ----------------------------------------------------
    // K1. Modify flux to account for geometry
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, grid)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend+1; k++){

      long int offset = NVAR*ni*(i + nj*j);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      // K1a. Multiply flux x area.  
      //      Avoid singularities in fluxes of controvariant quantities 
      //      (like Bth, Bphi) by simplifying metric and area factors.

      #if GEOMETRY != CARTESIAN
        #pragma acc loop seq
        NVAR_LOOP(nv) fA(k,nv) = flux(k,nv)*Ax3(i,j,k);

        #if GEOMETRY == CYLINDRICAL
        fA(k,MX2) *= x1[i];
        #if RADIATION
        fA(k,FR2) *= x1[i];
        #endif
        #if PHYSICS == RMHD
        // Simplify area as follows:
        // r*dr*dphi*(flux/r) = dr*dphi*flux

        IF (COVARIANT_BX2, A = DIM_EXPAND(dx1[i], *dx2[j], *1.0);
                           fA(k,BX2) = A*flux(k,BX2);)
        #endif
        #endif  // GEOMETRY == CYLINDRICAL

        #if GEOMETRY == SPHERICAL 
        IF (COVARIANT_MX2, fA(k,MX2) *= x1[i]);
        fA(k,MX3) *= x1[i]*sC[j];
        #if RADIATION
        IF (COVARIANT_MX2, fA(k,FR2) *= x1[i]);
        fA(k,FR3) *= x1[i]*sC[j];
        #endif
        #if PHYSICS == RMHD
        // Simplify area as follows:  
        // A*F(Bth)  = A*[flux/r]     = r*dr*dth*[flux/r]     = dr*dth*flux
        // A*F(Bphi) = A*[flux/(r*s)] = r*dr*dth*[flux/(r*s)] = dr*dth*flux/s        

        A  = DIM_EXPAND(dx1[i], *dx2[j], *1.0);
        IF (COVARIANT_BX2, fA(k,BX2) = flux(k,BX2)*A;) 
        fA(k,BX3) = flux(k,BX3)*A/sC[j]; 
        #endif
        #endif  // GEOMETRY == SPHERICAL

      #endif
    }}}
    #endif

    // ----------------------------------------------------
    // K2. Build rhs in the x3-dir
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(data, run_config, grid, target)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend; k++){

      long int offset = NVAR*ni*(i + nj*j);
      long int offset1 = ni*(i + nj*j);

      Ary1D Bn   (&data->sweep->Bn[offset1],ni);
      Ary1D PSIn (&data->sweep->PSIn[offset1],ni);
      Ary1D cmax (&Dts->cmax[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);

      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);

      double dtdx3 = dt/dx3[k];
      double v[NVAR];
      double local_rhs[NVAR];

      State state;
      Metric metric;

      STATE_SET (state, v, NULL, NULL, x1[i], x2[j], x3[k]);
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      double dl   = dx3[k]*metric.h3;

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(k,nv);

      // K2a. Conservative flux difference
      
      #if GEOMETRY == CARTESIAN

        #pragma acc loop seq
        NVAR_LOOP(nv) RHS_VALUE(k,nv) = -dtdx3*(flux(k+1,nv) - flux(k,nv));

      #else
      {
        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NVAR_LOOP(nv) RHS_VALUE(k,nv) = -dtdV*(fA(k+1,nv) - fA(k,nv));

        #ifdef GLM_MHD
        #if GEOMETRY == SPHERICAL
        double inv_rgrad = dx2[j]/(grid->rt[i]*grid->dmu[j]);
        RHS_VALUE(k,BXn) = -dtdx3*inv_rgrad*(flux(k+1,BXn) - flux(k,BXn));
        #else
        RHS_VALUE(k,BXn) = -dtdx3*(flux(k+1,BXn) - flux(k,BXn));
        #endif

        #endif
      }
      #endif

      // K2b. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      RHS_VALUE(k,MX3) += dt*v[RHO]*g[KDIR];
      IF_ENERGY(RHS_VALUE(k,ENG) += dt*0.5*(flux(k+1,RHO) + flux(k,RHO))*g[KDIR];)
      #endif  // (BODY_FORCE & VECTOR) 

      // K2e. Add div(B) source terms

      #if INCLUDE_DIVB_SOURCE
        double src[NVAR] = {0.0};  // Aggregate initialization
      
        #if GEOMETRY == CARTESIAN
        divB = (Bn(k+1) - Bn(k))/dx3[k];
        #else
        divB = (  Ax3(i,j,k+1)*Bn(k+1) 
                - Ax3(i,j,k)*Bn(k))/grid->dV(i,j,k);
        #endif
        #if DIVB_CONTROL == EIGHT_WAVES 
        PowellSourceTerm(dt*divB, v, src);
        #elif (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
        gradPSI = (PSIn(k+1) - PSIn(k))/dl;
        GLM_ExtendedSource (dt*divB, dt*gradPSI, v, src);
        #endif
        #if GEOMETRY != CARTESIAN
        OrthoToCons(src, metric);
        #endif

        #pragma acc loop seq 
        NFLX_LOOP(nv) RHS_VALUE(k,nv) -= src[nv];
      #endif // INCLUDE_DIVB_SOURCE

      // K2f. Compute inverse time step

      double q = 1.0;
      #if (RING_AVERAGE > 1) && (GEOMETRY == SPHERICAL)
      q = 1.0/grid->ring_av_csize[j];
      #endif
      data->C_dt(i,j,k) += q*0.5*(cmax(k) + cmax(k+1))/dl;
      if constexpr (fused) {
        #pragma acc loop seq
        NFLUID_LOOP(nv) target(i,j,k,nv) += local_rhs[nv];
      }
    }}}

  }

mynvtxstop_();
}
#undef RHS_VALUE

// Preserve implicit CPU RHSTarget callers as well as the legacy sweep output.
void RightHandSide(Data *data, timeStep *Dts, const Grid *const grid,
                   double dt, RBox *box, const int dir)
{
#if PCP_FOFC == YES
  if (fofc_cache::rhs_target) {
    RightHandSideImpl<true>(data, Dts, grid, dt, box, dir, *fofc_cache::rhs_target);
    return;
  }
#endif
  RightHandSideImpl<false>(data, Dts, grid, dt, box, dir, data->Uc);
}
#if PCP_FOFC == YES
void RightHandSideFused(Data *data, timeStep *Dts, const Grid *const grid,
                       double dt, RBox *box, int dir, const Ary4D &target)
{
  RightHandSideImpl<true>(data, Dts, grid, dt, box, dir, target);
}
#endif
