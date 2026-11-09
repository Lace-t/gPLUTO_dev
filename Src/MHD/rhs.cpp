////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Compute the right hand side of the conservative 
         HD/MHD equations.

  This function constructs the one-dimensional right hand sides of 
  the conservative MHD or HD equations in the direction given by dir 
  in different geometries.
  The right hand side in the \c d direction is computed as a two-point
  flux difference term plus a source term:  
  \f[
      R_{\bf i}^{(d)} = 
      -\frac{\Delta t}{\Delta V_{\bf i}}\
      \Big[  (AF)_{{\bf i} + \frac{1}{2}\mathbf{\hat e}_d} 
            -(AF)_{{\bf i} - \frac{1}{2}\mathbf{\hat e}_d} \Big]
         + \Delta t S_{\bf i}
  \f] 
   where \f${\bf i} = (i,j,k)\f$ while   
   - \f$ A        \f$: interface areas
   - \f$ {V} \f$: cell volume
   - \f$ F   \f$: interface fluxes
   - \f$ \Delta t \f$: time step
   - \f$ S \f$: source term including geometrical terms and 
                       body forces.

  See also the \ref RHS_page.
 
  \b References
     - "PLUTO: A Numerical Code for Computational Astrophysics."
       Mignone et al, ApJS (2007) 170, 228
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical Fluid
        Dynamics" 
       Mignone et al, ApJS (2012) 198, 7M
   
  \author A. Mignone  (andrea.mignone@unito.it)\n
          M. Rosssazza (marco.rossazza@unito.it)
  \date   Dec 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// -- Set macro for include/exclude div(B) source terms -- 

#if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
  #define INCLUDE_DIVB_SOURCE YES 
#elif (PHYSICS == MHD) && (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
  #define INCLUDE_DIVB_SOURCE YES 
#else
  #define INCLUDE_DIVB_SOURCE NO
#endif

// -- Compute the dyad BB -- 

#if PHYSICS == MHD
  #if BACKGROUND_FIELD == YES
    #define TOT_BB(v, b0, a, b) (v[a]*(v[b] + b0[b-BX1]) + b0[a-BX1]*v[b])
  #else
    #define TOT_BB(v, b0, a, b) (v[a]*v[b])
  #endif
#else
  #define TOT_BB(v, b0, a, b)  0.0
#endif

#if GEOMETRY == CYLINDRICAL
  #define iMPHI  MX2
  #define CYL_RAD(x1, s) fabs(x1)
#elif GEOMETRY == SPHERICAL 
  #define iMPHI  MX3
  #define CYL_RAD(x1, s) fabs((x1)*(s))
#endif

#ifndef iMPHI
 #define iMPHI MX2  // -- for Cartesian coordinates -- 
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Compute total (thermal + magnetic) pressure
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline double TotalPressure(const State &state)
{
  double ptot;
  double *v   = state.v;

  #if EOS == IDEAL 
  ptot = v[PRS];
  #elif EOS == ISOTHERMAL
  ptot = state.a2*v[RHO];
  #endif

  #if PHYSICS == MHD
  ptot += 0.5*DOT_PRODUCT(v + BX1, v + BX1);
  #if BACKGROUND_FIELD == YES
  double *Bgf = state.Bgf;
  ptot += DOT_PRODUCT(Bgf, v + BX1);
  #endif
  #endif
  return ptot;
}


#if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
////////////////////////////////////////////////////////////////////////
/*!
  Add Powell source terms to src[]
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void PowellSourceTerm(State &state, double* src) 
{
  double divB = state.divB;
  double *v   = state.v;
  #if BACKGROUND_FIELD == YES
  double *B0 = state.Bgf;
  double Bt[3];
  Bt[IDIR] = v[BX1] + B0[IDIR];
  Bt[JDIR] = v[BX2] + B0[JDIR];
  Bt[KDIR] = v[BX3] + B0[KDIR];
  #else
  double *Bt = v + BX1;
  #endif

  src[MX1] -= divB * Bt[IDIR];
  src[MX2] -= divB * Bt[JDIR];
  src[MX3] -= divB * Bt[KDIR];
  src[BX1] -= divB * v[VX1];
  src[BX2] -= divB * v[VX2];
  src[BX3] -= divB * v[VX3];

#if HAVE_ENERGY
  src[ENG] -= divB * DOT_PRODUCT(v + VX1, Bt);
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
  #if GEOMETRY != CARTESIAN
  const double *v   = state.v;
  #endif
  #if BACKGROUND_FIELD == YES
  const double *Bgf = state.Bgf;
  #endif
  #ifdef iVPHI
  double vphi = v[iVPHI] + state.w;
  #endif

  #if EOS == ISOTHERMAL  
  SoundSpeed2 (state, run_config);
  #else
  #endif

  #if GEOMETRY != CARTESIAN
  double ptot = TotalPressure(state);
  #endif

  #if GEOMETRY == CYLINDRICAL
  T[JDIR][JDIR] = ptot + v[RHO]*vphi*vphi - TOT_BB(v, Bgf, BX2, BX2);
  #elif GEOMETRY == SPHERICAL
  T[JDIR][JDIR] = ptot + v[RHO]*v[VX2]*v[VX2] - TOT_BB(v, Bgf, BX2, BX2);
  T[JDIR][IDIR] =        v[RHO]*v[VX2]*v[VX1] - TOT_BB(v, Bgf, BX2, BX1);
  T[KDIR][KDIR] = ptot + v[RHO]*vphi*vphi     - TOT_BB(v, Bgf, BX3, BX3);
  #elif GEOMETRY == CARTESIAN
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Compute right hand side for HD / MHD equations.
  
  \param[in,out] data  Pointer to Data structure 
  \param[in]     Dts   Pointer to timeStep structure 
  \param[in]     grid  Pointer to Grid structure
  \param[in]     dt    time increment
  \param[in]     box   RBox on which rhs has to be computed
  \param[in]     dir   sweeping direction
*/
////////////////////////////////////////////////////////////////////////
void RightHandSide (Data *data, timeStep *Dts, const Grid * const grid, 
                    double dt, RBox *box, const int dir)
{
  int i, j, k, nv;
  int ni, nj;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  double *x1  = grid->xC[IDIR], *x2  = grid->xC[JDIR], *x3  = grid->xC[KDIR];
  double *dx1 = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];
  double *x1m = grid->xL[IDIR], *x2m = grid->xL[JDIR], *x3m = grid->xL[KDIR];
  #if GEOMETRY == SPHERICAL
  double *sC = grid->sC, *sL  = grid->sL;
  #if PHYSICS == MHD
  double *dmu = grid->dmu;
  #endif
  #endif
  #if INCLUDE_DIVB_SOURCE
  double divB;
  #endif
  #if GEOMETRY != CARTESIAN
  const Ary3D &Ax1 = grid->A[IDIR];
  const Ary3D &Ax2 = grid->A[JDIR];
  const Ary3D &Ax3 = grid->A[KDIR];
  #endif
  RunConfig &run_config = *(data->run_config);
  #if (defined FARGO) || (ROTATING_FRAME == YES)
  double OmegaZ = run_config.OmegaZ;
  #endif
  #if (BODY_FORCE == POTENTIAL) & (HAVE_ENERGY)
  double phi_c;
  #endif
  static double *FA;
  #ifdef FARGO
  double **wA = FARGO_Velocity();
  #endif

  // 0. Allocate memory

  #if GEOMETRY != CARTESIAN
  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];
  if (FA == NULL) {
    FA   = ARRAY_1D(size_tot*NVAR, double);
  }
  #endif

  ARRAY_STRIDES (grid, dir, ni, nj);

  mynvtxstart_("RHS",GREEN);

  // ------------------------------------------------------
  // I. Compute rhs in the x1-direction
  // ------------------------------------------------------

  if (dir == IDIR){

    #if    (BODY_FORCE & POTENTIAL) || (GEOMETRY != CARTESIAN) \
        || (ROTATING_FRAME == YES)  || (defined FARGO)

    // ----------------------------------------------------
    // I1. Modify flux to account for geometry or gravity
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) //map(present: d[:1])
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){
    
      long int offset = NVAR*ni*(j + nj*k);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      State state;
      STATE_SET (state, NULL, NULL, NULL, x1m[i], x2[j], x3[k]);

      // I1a. Modify flux for ROTATING_FRAME or FARGO:
      //
      //      F(ENG)  += rho*vr*w^2/2 + w*rho*vphi*vr = F(RHO)*w^2/2 + w*F(MPHI)
      //      F(MPHI) += w*F(RHO) 

      #if (defined FARGO) || (ROTATING_FRAME == YES)
      double wm = 0.0;
      IF_FARGO(wm += 0.5*(FARGO_SPEED(wA, i-1,j,k) + FARGO_SPEED(wA, i,j,k)));
      IF (ROTATING_FRAME, double R = CYL_RAD(grid->rL[i], grid->sC[j]);
                          wm += OmegaZ*R;)
      IF_ENERGY  (flux(i,ENG) += wm*(0.5*wm*flux(i,RHO) + flux(i,iMPHI));)
      flux(i,iMPHI) += wm*flux(i,RHO);
      #endif

      // I1c. Add grav. potential to energy flux

      #if BODY_FORCE & POTENTIAL
      #pragma omp declare target
      data->phi_m(i,j,k) = BodyForcePotential(state, run_config);
      #pragma omp end declare target
      IF_ENERGY(flux(i,ENG) += flux(i,RHO)*data->phi_m(i,j,k);)
      #endif

      //  I1d. Multiply flux x area. 
      //       Avoid singularities in fluxes of controvariant quantities 
      //       (like Bth, Bphi) by simplifying metric and area factors.

      #if GEOMETRY != CARTESIAN

        #pragma acc loop seq
        NVAR_LOOP(nv) fA(i,nv) = flux(i,nv)*Ax1(i,j,k);

        #if GEOMETRY == CYLINDRICAL
        fA(i,MX2) *= fabs(x1m[i]);
        #if PHYSICS == MHD
        // Simplify area as follows:
        // A*f(BX2) = A*[flux/r] = r*dphi*dz*[flux/r] = dphi*dz*flux  
        double A;
        IF (COVARIANT_BX2, A = DIM_EXPAND(1.0, *dx2[j], *dx3[k]); 
                           fA(i,BX2) = flux(i,BX2)*A;)  
        #endif
        #endif

        #if GEOMETRY == SPHERICAL
        IF (COVARIANT_MX2, fA(i,MX2) *= fabs(x1m[i]));
        fA(i,MX3) *= fabs(x1m[i]*sC[j]);
        #if PHYSICS == MHD
        // Simplify area as follows:
        // A*F(Bth)  = A*[flux/r] = r^2*dmu*dphi*[flux/r] = r*dmu*dphi*flux        
        // A*F(Bphi) = A*[flux/(r*s)] = r^2*dmu*dphi*[flux/(r*s)] 
        //                            = r*dmu*dphi*flux/s      

        double A = DIM_EXPAND(x1m[i], *dmu[j], *dx3[k]);
        IF (COVARIANT_BX2, fA(i,BX2) = flux(i,BX2)*fabs(A)); 
        fA(i,BX3)  = flux(i,BX3)*fabs(A/sC[j]);  
        #endif
        #endif

      #endif // GEOMETRY != CARTESIAN 
    }}}
    #endif

    // ----------------------------------------------------
    // I2. Build rhs in the x1-dir
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){   
    
      long int offset1 = ni*(j + nj*k);
      long int offset  = NVAR*offset1;

      Ary1D Bn(&data->sweep->Bn[offset1],ni);
      Ary1D PSIn(&data->sweep->PSIn[offset1],ni);
      Ary1D cmax(&Dts->cmax[offset1],ni);
      #if RESISTIVITY
      long int offset2 = dir*offset1;
      Ary1D Dmax(&Dts->Dmax[offset2],ni);
      #endif
      Ary1D press(&data->sweep->press[offset1],ni);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary2D fA(&FA[offset],ni,NVAR);
      #endif
      #if SELF_GRAV
      Ary1D g_c(&data->sweep->g_c[offset1],ni);
      #endif

      double dtdx1 = dt/dx1[i];
      double v[NVAR];
      double Bgf[3]; // For BackgroundField, if needed
      double src[NVAR] = {0.0};  // Aggregate initialization
      #if GEOMETRY != CARTESIAN
      double T[3][3];
      Metric metric;
      #endif
      State state;

      STATE_SET (state, v, NULL, Bgf, x1[i], x2[j], x3[k]);
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(i,nv);

      #if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
      BackgroundField (state, run_config); 
      #endif

      // I2a. Conservative flux difference

      #if GEOMETRY == CARTESIAN

        #pragma acc loop seq
        NFLUID_LOOP(nv)  rhs(i,nv) = -dtdx1*(flux(i+1,nv) - flux(i,nv));
        
        #if (defined FARGO)
        double w = wA[k][i];
        rhs(i,MX2) -= w*rhs(i,RHO);
        IF_ENERGY(rhs(i,ENG) -= w*(rhs(i,MX2) + 0.5*w*rhs(i,RHO));)
        #endif

      #else

        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NFLUID_LOOP(nv) rhs(i,nv) = -dtdV*(fA(i+1,nv) - fA(i,nv));

      #endif

      // I2b. Modify rhs to enforce conservation

      #if GEOMETRY != CARTESIAN        
        #if (defined FARGO) || (ROTATING_FRAME == YES)
        double R = CYL_RAD(grid->rC[i], grid->sC[j]);
        double w = 0.0;
        IF_FARGO (w += FARGO_SPEED(wA,i,j,k));
        IF (ROTATING_FRAME, w += R*OmegaZ);

        rhs(i, iMPHI) -= w*rhs(i,RHO)*R;
        IF_ENERGY(rhs(i,ENG) -= 0.5*w*w*rhs(i,RHO) + (w/R)*rhs(i,iMPHI);)

        state.w = w;
        #endif // (defined FARGO) || (ROTATING_FRAME == YES) 
      #endif

      // I2c. Add source term: gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      src[MX1]            += dt*v[RHO]*g[IDIR];
      IF_ENERGY (src[ENG] += dt*0.5*(flux(i+1,RHO) + flux(i,RHO))*g[IDIR];)
      #endif

      #if (BODY_FORCE & POTENTIAL)
      src[MX1]           -= dtdx1*v[RHO]*(data->phi_m(i+1,j,k) - data->phi_m(i,j,k));
      IF_ENERGY(phi_c     = BodyForcePotential(state, run_config);
                src[ENG] -= phi_c*rhs(i,RHO); )
      #endif

      #if (SELF_GRAV)
      rhs(i,MX1)   += dt*v[RHO]*g_c(i);
      IF_ENERGY (rhs(i,ENG) += dt*0.5*(flux(i+1,RHO) + flux(i,RHO))*g_c(i);)
      #endif

      // I2d. Add geometrical source terms

      #if GEOMETRY == CYLINDRICAL    // ** CYLINDRICAL GEOMETRY ** 
        double r     = x1[i];
        TensorSource (state, run_config, T);

        src[MX1] += dt*T[JDIR][JDIR]/r;
        #ifdef GLM_MHD
        src[BX1] += dt*v[PSI_GLM]/r;
        #endif

        #if PHYSICS == MHD && (!COVARIANT_BX2)
        {
          double Ez = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]); 
          rhs(i,BX2) -= dt*Ez/r;
        }
        #endif

      #elif GEOMETRY == SPHERICAL    // ** SPHERICAL GEOMETRY **

        double r     = grid->rt[i]; // x1[i];
        TensorSource (state, run_config, T);
        src[MX1] += dt*(T[JDIR][JDIR] + T[KDIR][KDIR])/r;
        #if !INCLUDE_JDIR
        src[MX1] -= dt*T[JDIR][IDIR]*grid->cot[j]/r;
        #endif

        #if COVARIANT_MX2
        src[MX2] += dt*T[KDIR][KDIR]*grid->cot[j]/r;
        #if !INCLUDE_JDIR
        src[MX2] -= dt*T[JDIR][JDIR]*grid->cot[j]/r;
        #endif
        #else
        src[MX2] += dt*(T[KDIR][KDIR]*grid->cot[j] - T[JDIR][IDIR])/r;
        #if !INCLUDE_JDIR
        src[MX2] -= dt*T[JDIR][JDIR]*grid->cot[j]/r;
        #endif
        #endif

        #ifdef GLM_MHD
        src[BX1] += dt*2.0*v[PSI_GLM]/r;
        #if INCLUDE_JDIR
        src[BX2] += dt*v[PSI_GLM]*grid->cot[j]/r;
        #endif
        #endif

        #if PHYSICS == MHD && (!COVARIANT_BX2)
        {
          double Ez = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]); 
          src[BX2] -= dt*Ez/r;
        }
        #endif

      #endif // GEOMETRY == SPHERICAL

      // I2e. Add div(B) source terms

      #if INCLUDE_DIVB_SOURCE
        #if GEOMETRY == CARTESIAN
        state.divB = divB = (Bn(i+1) - Bn(i))/dx1[i];
        #else
        state.divB = divB = (Ax1(i+1,j,k)*Bn(i+1) - Ax1(i,j,k)*Bn(i))/grid->dV(i,j,k);
        #endif
        #if DIVB_CONTROL == EIGHT_WAVES 
        state.divB *= dt;
        PowellSourceTerm(state, src);

        #elif (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
        double gradPSI = (PSIn(i+1) - PSIn(i))/dx1[i];
        GLM_ExtendedSource (dt*divB, dt*gradPSI, v, src, dir);
        #endif
      #endif // INCLUDE_DIVB_SOURCE

      //  I2f. Transfrom to conservative and add source terms

      #if GEOMETRY != CARTESIAN
      OrthoToCons(src, metric);
      #endif

      #pragma acc loop seq 
      NFLX_LOOP(nv) rhs(i,nv) += src[nv];   

      // I2g. Compute inverse time step (we assume metric.h1 = 1)
      
      data->C_dt(i,j,k) = 0.5*(cmax(i) + cmax(i+1))/dx1[i];
      #if RESISTIVITY
        data->C_dtpar[IDIR](i,j,k) = 0.5*(Dmax(i) + Dmax(i+1))/(dx1[i]*dx1[i]);
      #endif
    }}}

  }else if (dir == JDIR){

    // ----------------------------------------------------
    // J. Compute rhs in the x2-direction
    // ----------------------------------------------------

    #if (BODY_FORCE & POTENTIAL) || (GEOMETRY != CARTESIAN)

    // ----------------------------------------------------
    // J1. Modify flux to account for geometry or gravity
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend+1; j++){

      long int offset = NVAR*ni*(i + nj*k);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      State state;

      STATE_SET (state, NULL, NULL, NULL, x1[i], x2m[j], x3[k]);

      // J1a. Modify flux for ROTATING_FRAME or FARGO:
      //
      //      F(ENG)  += rho*vth*w^2/2 + w*rho*vphi*vth = F(RHO)*w^2/2  + w*F(MPHI)
      //      F(MPHI) += w*F(RHO) 

      #if GEOMETRY == SPHERICAL
      #if (defined FARGO) || (ROTATING_FRAME == YES)
      double wm = 0.0;
      IF_FARGO(            wm += 0.5*(wA[j-1][i] + wA[j][i]);)
      IF (ROTATING_FRAME,  double R  = grid->rC[i]*grid->sL[j];
                           wm += OmegaZ*R;)
      IF_ENERGY  (flux(j,ENG) += wm*(0.5*wm*flux(j,RHO) + flux(j,iMPHI));)
      flux(j,iMPHI) += wm*flux(j,RHO);
      #endif
      #endif

      // J1c. Add grav. potential to energy flux

      #if BODY_FORCE & POTENTIAL
      data->phi_m(i,j,k) = BodyForcePotential(state, run_config);
      IF_ENERGY(flux(j,ENG) += flux(j,RHO)*data->phi_m(i,j,k);)
      #endif
    
      // J1d. Multiply flux x area
      //      Avoid singularities in fluxes of controvariant quantities 
      //      (like Bth, Bphi) by simplifying metric and area factors.
      
      #if GEOMETRY != CARTESIAN

        #pragma acc loop seq
        NVAR_LOOP(nv) fA(j,nv) = flux(j,nv)*Ax2(i,j,k);

        #if GEOMETRY == CYLINDRICAL
        fA(j,MX2) *= fabs(x1[i]);
        #if PHYSICS == MHD
        IF (COVARIANT_BX2, fA(j,BX2) /= fabs(x1[i])); // A*F(BX2) = A*[flux/r] 
                                                      //          = dr*dz*[flux/r] 
        #endif
        #endif

        #if GEOMETRY == SPHERICAL
        IF (COVARIANT_MX2, fA(j,MX2) *= fabs(x1[i]));
        fA(j,MX3) *= fabs(x1[i]*sL[j]);
        #if PHYSICS == MHD
        // A*F(Bth) = A*[flux/r] = r*dr*s*dphi*flux/r 
                                  // =   dr*s*dphi*flux   
            // A*F(Bphi) = A*[flux/(r*s)] = r*dr*s*dphi*flux/(r*s) 
                                      // =   dr*dphi*flux           
        double A = DIM_EXPAND(dx1[i],  , *dx3[k]);
        IF (COVARIANT_BX2, fA(j,BX2) = flux(j,BX2)*fabs(A*sL[j]));
        fA(j,BX3) = flux(j,BX3)*fabs(A); 
        #endif
        #endif

      #endif 
    }}}
    #endif

    // ----------------------------------------------------
    // J2. Build rhs in the x2-dir
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){  

      long int offset1 = ni*(i + nj*k);
      long int offset  = NVAR*offset1;

      Ary1D Bn(&data->sweep->Bn[offset1],ni);
      Ary1D PSIn(&data->sweep->PSIn[offset1],ni);
      Ary1D cmax(&Dts->cmax[offset1],ni);
      #if RESISTIVITY
      long int offset2 = dir*offset1;
      Ary1D Dmax(&Dts->Dmax[offset2],ni);
      #endif
      Ary1D press(&data->sweep->press[offset1],ni);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary2D fA(&FA[offset],ni,NVAR);
      #endif
      #if SELF_GRAV
      Ary1D g_c(&data->sweep->g_c[offset1],ni);
      #endif

      double v[NVAR];
      double Bgf[3]; // For BackgroundField, if needed
      double src[NVAR] = {0.0};  // Aggregate initialization
      State state;
      Metric metric;

      STATE_SET (state, v, NULL, Bgf, x1[i], x2[j], x3[k]);
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      double dl   = dx2[j]*metric.h2;

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(j,nv);

      #if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
      BackgroundField (state, run_config); 
      #endif

      // J2a. Conservative flux difference

      #if GEOMETRY == CARTESIAN
        double dtdx2 = dt/dx2[j];
        
        #pragma acc loop seq
        NFLUID_LOOP(nv) rhs(j,nv) = -dtdx2*(flux(j+1,nv) - flux(j,nv));

      #else 

        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NFLUID_LOOP(nv) rhs(j,nv) = -dtdV*(fA(j+1,nv) - fA(j,nv));

      #endif  // GEOMETRY == CARTESIAN

      // J2b. Modify rhs to enforce conservation

      #if GEOMETRY == SPHERICAL
        #if (defined FARGO) || (ROTATING_FRAME == YES)
        double R = CYL_RAD(grid->rC[i], grid->sC[j]);
        double w = 0.0;
        IF_FARGO (w += FARGO_SPEED(wA,i,j,k));
        IF (ROTATING_FRAME, w += R*OmegaZ);

        rhs(j,MX3) -= w*rhs(j,RHO)*R;
        IF_ENERGY(rhs(j,ENG)  -= 0.5*w*w*rhs(j,RHO) + (w/R)*rhs(j,iMPHI);)

        state.w = w;
        #endif
      #endif

      // J2c. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      src[MX2] += dt*v[RHO]*g[JDIR];
      IF_ENERGY(src[ENG] += dt*0.5*(flux(j+1,RHO) + flux(j,RHO))*g[JDIR];)
      #endif 

      #if (BODY_FORCE & POTENTIAL)
      double dtdl = dt/dl;
      src[MX2]           -= dtdl*v[RHO]*(data->phi_m(i,j+1,k) - data->phi_m(i,j,k));
      IF_ENERGY(phi_c     = BodyForcePotential(state, run_config);
                src[ENG] -= phi_c*rhs(j,RHO);)
      #endif

      #if (SELF_GRAV)
      rhs(j,MX2)   += dt*v[RHO]*g_c(j);
      IF_ENERGY (rhs(j,ENG) += dt*0.5*(flux(j+1,RHO) + flux(j,RHO))*g_c(j);)
      #endif

      // J2d. Add geometrical source terms

      #if GEOMETRY == SPHERICAL
      #endif

      // J2e. Add div(B) source terms

      #if INCLUDE_DIVB_SOURCE
        #if GEOMETRY == CARTESIAN
        divB = state.divB = (Bn(j+1) - Bn(j))/dx2[j];
        #else
        divB = state.divB = (Ax2(i,j+1,k)*Bn(j+1) - Ax2(i,j,k)*Bn(j))/grid->dV(i,j,k);
        #endif
        #if DIVB_CONTROL == EIGHT_WAVES 
        state.divB *= dt;
        PowellSourceTerm(state, src);
        #elif (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
        double gradPSI = (PSIn(j+1) - PSIn(j))/dl;
        GLM_ExtendedSource (dt*divB, dt*gradPSI, v, src,dir);
        #endif
      #endif // INCLUDE_DIVB_SOURCE

      // J2f. Transfrom to conservative and add source terms

      #if GEOMETRY != CARTESIAN
      OrthoToCons(src, metric);
      #endif

      #pragma acc loop seq 
      NFLX_LOOP(nv) rhs(j,nv) += src[nv];   

      // J2g. Compute inverse time step 

      double q = 1.0;
      #if (RING_AVERAGE > 1) && (GEOMETRY == CYLINDRICAL)
      q = 1.0/grid->ring_av_csize[i];
      #endif
      data->C_dt(i,j,k) += q*0.5*(cmax(j) + cmax(j+1))/dl;
      #if RESISTIVITY
        data->C_dtpar[JDIR](i,j,k) = q*0.5*(Dmax(j) + Dmax(j+1))/(dl*dl);
      #endif
    }}}
    
  }else if (dir == KDIR){

    // ----------------------------------------------------
    // K. Compute rhs in the x3-direction
    // ----------------------------------------------------

    #if (BODY_FORCE & POTENTIAL) || (GEOMETRY != CARTESIAN) || (defined FARGO)

    // ----------------------------------------------------
    // K1. Modify flux to account for geometry, gravity 
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend+1; k++){

      long int offset = NVAR*ni*(i + nj*j);
      Ary2D fA(&FA[offset],ni,NVAR);
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);

      State state;
      STATE_SET (state, NULL, NULL, NULL, x1[i], x2[j], x3m[k]);

      // K1a. Modify flux for ROTATING_FRAME or FARGO

      #if (GEOMETRY != SPHERICAL) && (defined FARGO)
      double wm = 0.5*(wA[k-1][i] + wA[k][i]);
      IF_ENERGY(flux(k,ENG) += wm*(0.5*wm*flux(k,RHO) + flux(k,iMPHI));)
      flux(k,iMPHI) += wm*flux(k,RHO);
      #endif

      // K1c. Add grav. potential to energy flux

      #if BODY_FORCE & POTENTIAL
      data->phi_m(i,j,k) = BodyForcePotential(state, run_config);
      IF_ENERGY(flux(k,ENG) += flux(k,RHO)*data->phi_m(i,j,k);)
      #endif
    
      // K1d. Multiply flux x area 
      //      Avoid singularities in fluxes of controvariant quantities 
      //      (like Bth, Bphi) by simplifying metric and area factors.

      #if GEOMETRY != CARTESIAN

        #pragma acc loop seq
        NVAR_LOOP(nv) fA(k,nv) = flux(k,nv)*Ax3(i,j,k);

        #if GEOMETRY == CYLINDRICAL
        fA(k,MX2) *= fabs(x1[i]);
        #if PHYSICS == MHD
        // Simplify area as follows:
        // r*dr*dphi*(flux/r) = dr*dphi*flux
        double A;
        IF (COVARIANT_BX2, A = DIM_EXPAND(dx1[i], *dx2[j], *1.0);
                           fA(k,BX2) = A*flux(k,BX2);)
        #endif
        #endif

        #if GEOMETRY == SPHERICAL 
        IF (COVARIANT_MX2, fA(k,MX2) *= fabs(x1[i]));
        fA(k,MX3) *= fabs(x1[i]*sC[j]);
        #if PHYSICS == MHD
        // A*F(Bth)  = A*[flux/r] = r*dr*dth*[flux/r]
                                  // =   dr*dth*flux
            // A*F(Bphi) = A*[flux/(r*s)] = r*dr*dth*[flux/(r*s)]
                                      // =   dr*dth*flux/s         
        double A  = DIM_EXPAND(dx1[i], *dx2[j], *1.0);
        IF (COVARIANT_BX2, fA(k,BX2) = flux(k,BX2)*A;) 
        fA(k,BX3) = flux(k,BX3)*fabs(A/sC[j]); 
        #endif
        #endif

      #endif
    }}}
    #endif

    // ----------------------------------------------------
    // K2. Build rhs in the x3-dir
    // ----------------------------------------------------

    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, run_config, grid)
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
    for (k = kbeg; k <= kend; k++){

      long int offset1 = ni*(i + nj*j);
      long int offset  = NVAR*offset1;

      Ary1D cmax(&Dts->cmax[offset1],ni);
      #if RESISTIVITY
      long int offset2 = dir*offset1;
      Ary1D Dmax(&Dts->Dmax[offset2],ni);
      #endif
      Ary1D Bn(&data->sweep->Bn[offset1],ni);
      Ary1D PSIn(&data->sweep->PSIn[offset1],ni);
      Ary1D press(&data->sweep->press[offset1],ni);
      #if GEOMETRY != CARTESIAN
      Ary2D fA(&FA[offset],ni,NVAR);
      #endif
      Ary2D flux(&data->sweep->flux[offset],ni,NVAR);
      Ary2D rhs(&data->sweep->rhs[offset],ni,NVAR);
      Ary2D vc(&data->sweep->vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep->vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep->vR[offset],ni,NVAR);
      #if SELF_GRAV
      Ary1D g_c(&data->sweep->g_c[offset1],ni);
      #endif
      
      double v[NVAR];
      double Bgf[3]; // For BackgroundField, if needed
      double src[NVAR] = {0.0};  // Aggregate initialization

      State state;
      Metric metric;

      STATE_SET (state, v, NULL, Bgf, x1[i], x2[j], x3[k]);
      #if GEOMETRY != CARTESIAN
      MetricSet (metric, grid->rC[i], grid->sC[j]);
      #endif

      double dl   = dx3[k]*metric.h3;

      #pragma acc loop seq
      NVAR_LOOP(nv) v[nv] = vc(k,nv);

      #if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
      BackgroundField (state, run_config); 
      #endif

      // K2a. Conservative flux difference
      
      #if GEOMETRY == CARTESIAN
        double dtdx3 = dt/dx3[k];

        #pragma acc loop seq
        NFLUID_LOOP(nv) rhs(k,nv) = -dtdx3*(flux(k+1,nv) - flux(k,nv));

        #if (GEOMETRY == CARTESIAN) && (defined FARGO)
        double w = wA[k][i];
        rhs(k,MX2) -= w*rhs(k,RHO);
        IF_ENERGY(rhs(k,ENG) -= w*(rhs(k,MX2) + 0.5*w*rhs(k,RHO));)
        #endif

      #else

        double dtdV = dt/grid->dV(i,j,k);
        #pragma acc loop seq
        NFLUID_LOOP(nv) rhs(k,nv) = -dtdV*(fA(k+1,nv) - fA(k,nv));

      #endif

      // K2b. modify rhs to enforce conservation (FARGO only)
      //      (solid body rotations are not included the
      //      velocity depends on the cylindrical radius only)

      #if (GEOMETRY == CYLINDRICAL) && (defined FARGO)
      double R = CYL_RAD(grid->rC[i], grid->sC[j]);
      double w = wA[k][i];
      rhs(k,MX2) -= w*rhs(k,RHO)*R;
      IF_ENERGY(rhs(k,ENG) -= w*(rhs(k,MX2)/R + 0.5*w*rhs(k,RHO));)
      #endif

      // K2c. Add gravity

      #if (BODY_FORCE & VECTOR)
      double g[3];
      BodyForceVector(state, run_config, g);

      src[MX3] += dt*v[RHO]*g[KDIR];
      IF_ENERGY(src[ENG] += dt*0.5*(flux(k+1,RHO) + flux(k,RHO))*g[KDIR];)
      #endif  // (BODY_FORCE & VECTOR) 

      #if (BODY_FORCE & POTENTIAL)
      double dtdl = dt/dl;
      src[MX3]           -= dtdl*v[RHO]*(data->phi_m(i,j,k+1) - data->phi_m(i,j,k));
      IF_ENERGY(phi_c     = BodyForcePotential(state, run_config);
                src[ENG] -= phi_c*rhs(k,RHO);)
      #endif      

      #if (SELF_GRAV)
      rhs(k,MX3)   += dt*v[RHO]*g_c(k);
      IF_ENERGY (rhs(k,ENG) += dt*0.5*(flux(k+1,RHO) + flux(k,RHO))*g_c(k);)
      #endif

      // K2e. Add div(B) source terms

      #if INCLUDE_DIVB_SOURCE
        #if GEOMETRY == CARTESIAN
        divB = state.divB = (Bn(k+1) - Bn(k))/dx3[k];
        #else
        divB = state.divB = (Ax3(i,j,k+1)*Bn(k+1) - Ax3(i,j,k)*Bn(k))/grid->dV(i,j,k);
        #endif
        #if DIVB_CONTROL == EIGHT_WAVES 
        state.divB *= dt;
        PowellSourceTerm(state, src);
        #elif (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
        double gradPSI = (PSIn(k+1) - PSIn(k))/dl;
        GLM_ExtendedSource (dt*divB, dt*gradPSI, v, src, dir);
        #endif
      #endif // INCLUDE_DIVB_SOURCE

      // K2f. Transfrom to conservative and add  source terms

      #if GEOMETRY != CARTESIAN
      OrthoToCons(src, metric);
      #endif

      #pragma acc loop seq 
      NFLX_LOOP(nv) rhs(k,nv) += src[nv];   

      // K2g. Compute inverse time step 

      double q = 1.0;
      #if (RING_AVERAGE > 1) && (GEOMETRY == SPHERICAL)
      q = 1.0/grid->ring_av_csize[j];
      #endif
      data->C_dt(i,j,k) += q*0.5*(cmax(k) + cmax(k+1))/dl;
      #if RESISTIVITY
        data->C_dtpar[KDIR](i,j,k) = q*0.5*(Dmax(k) + Dmax(k+1))/(dl*dl);
      #endif
    }}}
  }
  
mynvtxstop_();
}
