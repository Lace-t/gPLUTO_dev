////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Set labels, indexes and prototypes for the MHD module.

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the MHD module.

  \author M. Rossazza (marco.rossazza@unito.it) 
          A. Mignone  (andrea.mignone@unito.it)
  \date   Oct 23, 2025
*/
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/* Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.              */
////////////////////////////////////////////////////////////////////////

#define RHO  0
#define MX1  1
#define MX2  2
#define MX3  3
#define ENG  4

#define VX1  MX1
#define VX2  MX2
#define VX3  MX3
#define PRS  ENG

#define NFLX (4 + HAVE_ENERGY)

////////////////////////////////////////////////////////////////////////
/* Now define more convenient and user-friendly 
     pointer labels for geometry setting                              */
////////////////////////////////////////////////////////////////////////

#if GEOMETRY == CYLINDRICAL 

 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3

 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3

#endif

#if GEOMETRY == SPHERICAL 

 #define iVR    VX1
 #define iVTH   VX2
 #define iVPHI  VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3

#endif

#ifndef COVARIANT_MX2 
 #define COVARIANT_MX2 NO
#endif

// -----------------------------------------------------
// Function prototype
// -----------------------------------------------------

#pragma omp declare target
#pragma acc routine seq
uint32_t RHD_EnergySolve  (double *, double *, RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
uint32_t RHD_EntropySolve (double *, double *, RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
uint32_t RHD_PressureFix  (double *, double *, RunConfig &);
#pragma omp end declare target

void LF_Solver_T       (Data *, timeStep *, Grid *, RBox *, int);
void HLL_Solver_T      (Data *, timeStep *, Grid *, RBox *, int);
void HLLC_Solver_T     (Data *, timeStep *, Grid *, RBox *, int);
void TwoShock_Solver_T (Data *, timeStep *, Grid *, RBox *, int);

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Include additional header files                                   */
////////////////////////////////////////////////////////////////////////

#include "eos.hpp"
#include "vel_limiter.hpp"

////////////////////////////////////////////////////////////////////////
/* Inline functions                                                   */
////////////////////////////////////////////////////////////////////////

// Function prototype
template<int dir> 
inline void Flux (State &, double *, double &);

// Explicit template instantiations
template void Flux<IDIR>(State &, double *, double &);
template void Flux<JDIR>(State &, double *, double &);
template void Flux<KDIR>(State &, double *, double &);

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void Flux (State &state, double *flux, double &p)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  constexpr int VXn = VX1 + dir;
  constexpr int MXn = MX1 + dir;
  const double *v = state.v;
  const double *u = state.u;

  flux[RHO] = u[RHO]*v[VXn];
  flux[MX1] = u[MX1]*v[VXn];
  flux[MX2] = u[MX2]*v[VXn];
  flux[MX3] = u[MX3]*v[VXn];
  flux[ENG] = u[MXn];

  flux[MXn] += v[PRS];
  
  p = v[PRS];
}

// Function prototype
template<int dir> 
inline void MaxSignalSpeed (State &, double &, double &);

// Explicit template instantiations
template void MaxSignalSpeed<IDIR>(State &, double &, double &);
template void MaxSignalSpeed<JDIR>(State &, double &, double &);
template void MaxSignalSpeed<KDIR>(State &, double &, double &);

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void MaxSignalSpeed (State &state, double &cmin, double &cmax)
/*!
  Compute the maximum and minimum characteristic velocities for the
  RMHD equations, cmin= v - cf, cmax = v + cf
  
*/
////////////////////////////////////////////////////////////////////////
{
  double vx, vt2, vel2;
  double sroot;
  const double *v  = state.v;
  const double cs2 = state.a2;

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  vx   = v[VXn];
  vt2  = v[VXt]*v[VXt] + v[VXb]*v[VXb];
  vel2 = vx*vx + vt2;

  sroot = sqrt(cs2*(1.0 - vx*vx - vt2*cs2)*(1.0 - vel2));

  cmax = (vx*(1.0 - cs2) + sroot)/(1.0 - vel2*cs2);
  cmin = (vx*(1.0 - cs2) - sroot)/(1.0 - vel2*cs2);
}

// Function prototype 
template<int dir> 
inline void HLL_Speed (State &, State &, double &, double &);

// Explicit template instantiations
template void HLL_Speed<IDIR> (State &, State &, double &, double &);
template void HLL_Speed<JDIR> (State &, State &, double &, double &);
template void HLL_Speed<KDIR> (State &, State &, double &, double &);

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void HLL_Speed (State &stateL, State &stateR, double &SL, double &SR)
/*!
  Compute leftmost (SL) and rightmost (SR) speeds for the Riemann fan.
  
*/
////////////////////////////////////////////////////////////////////////
{
  double sl_min, sl_max;
  double sr_min, sr_max;

  // ----------------------------------------------
  // 1. Compute speeds using DAVIS Estimate
  // ----------------------------------------------

  MaxSignalSpeed<dir> (stateL, sl_min, sl_max);
  MaxSignalSpeed<dir> (stateR, sr_min, sr_max);

  SL = MIN(sl_min, sr_min);
  SR = MAX(sl_max, sr_max);
}

// Function prototype
template<int dir> 
inline void PrimEigenvectors (State &, RunConfig &, 
                              double LL[][NFLX], double RR[][NFLX], 
                              double *);
// Explicit template instantiations
template void PrimEigenvectors<IDIR> (State &, RunConfig &, 
                                      double LL[][NFLX], double RR[][NFLX],
                                      double *);
template void PrimEigenvectors<JDIR> (State &, RunConfig &, 
                                      double LL[][NFLX], double RR[][NFLX],
                                      double *);
template void PrimEigenvectors<KDIR> (State &, RunConfig &, 
                                      double LL[][NFLX], double RR[][NFLX],
                                      double *);
////////////////////////////////////////////////////////////////////////
/*!
  Provide eigenvectors eigenvalues of the relativistic
  equations in primitive form.
  
  \param [in]   state       a State structure
  \param [in]   run_config  reference to a RunConfig structure
  \param [out]  LL          left eigenvectors
  \param [out]  RR          right eigenvectors
  \param [out]  lambda      eigenvalues
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir> 
inline void  PrimEigenvectors (State &state, RunConfig &run_config,
                               double LL[][NFLX], double RR[][NFLX],
                               double *lambda)
{
  int    k,l,nv;
  double cs, g, g2, eta;
  double D, g2_1d;
  double a,b;
  double u, v, w, rp,rm, v2tan;
  double cs2, h;
  double *vp = state.v;

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  // 1. Reset eigenvectors

  for (nv = 0; nv < NFLX; nv++){
  for (l = 0; l < NFLX; l++){
    LL[nv][l] = 0.0;
    RR[nv][l] = 0.0;
  }}

  // 2a. Compute thermodynamics and other recurring
  //    quantities

  SoundSpeed2 (state,run_config);
  cs2 = state.a2;
  h   = state.h;

  u = vp[VXn]; // Three-velocities
  v = vp[VXt];
  w = vp[VXb];

  cs    = sqrt(cs2);
  v2tan = v*v + w*w;
  eta   = sqrt(1.0 - u*u - cs2*v2tan);
  g2_1d = 1.0/(1.0 - u*u);
  g2    = u*u + v2tan;
  a     = 1.0 - g2*cs2;
  g2    = 1.0/(1.0 - g2);
  g     = sqrt(g2);
  D     = vp[RHO]*g;                  // Lab - Density 
  rm    = u*(1.0 - cs2) - cs*eta/g;
  rp    = u*(1.0 - cs2) + cs*eta/g;

  // 2b. Define eigenvalues 

  lambda[0] = rm/a;
  lambda[1] = rp/a;

  for (k = 2; k < NFLX; k++) lambda[k] = u;

  // 3. Compute Right eigenvectors. 
  
                                        // Equivalent to
  rp = cs*rp/(vp[RHO]*(u*cs - g*eta));  // rp = -cs*(g*eta*u + cs)/(g*D*(1.0 - u*u));
  rm = cs*rm/(vp[RHO]*(u*cs + g*eta));  // rm =  cs*(g*eta*u - cs)/(g*D*(1.0 - u*u));

  a = cs*eta/D;
  b = h*cs2;

  // 3a. Right eigenvector, lambda = u - cs 

  RR[RHO][0]  = 1.0;
  RR[VXn][0] = -a;
  RR[VXt][0] = v*rm;
  RR[VXb][0] = w*rm;
  RR[PRS][0] = b;

  // 3b. Right eigenvector, lambda = u + cs 

  RR[RHO][1] = 1.0;
  RR[VXn][1] = a;
  RR[VXt][1] = v*rp;
  RR[VXb][1] = w*rp;
  RR[PRS][1] = b;

  // 3c. Right eigenvector, lambda = u

  RR[RHO][2] = 1.0;
  RR[VXt][3] = 1.0;
  RR[VXb][4] = 1.0;

  // 4. Left eigenvectors

  a = 0.5/a;
  b = 0.5/b;

  // 4a. Left eigenvector, lambda = u - cs 

  LL[0][VXn] = -a;
  LL[0][PRS] =  b;

  // 4b. Left eigenvector, lambda = u + cs

  LL[1][VXn] = a;
  LL[1][PRS] = b;

  // 4c. Left eigenvector, lambda = u

  LL[2][RHO] =  1.0;
  LL[2][PRS] = -2.0*b;

  // 4d. Left eigenvector, lambda = u (transverse vel)

  LL[3][VXn] = u*v*g2_1d;
  LL[3][PRS] = v*g2_1d/(g2*vp[RHO]*h);

  LL[4][VXn] = u*w*g2_1d;
  LL[4][VXb] = 1.0;
  LL[4][PRS] = w*g2_1d/(g2*vp[RHO]*h);

}

// Function prototype
template<int dir> 
inline void PrimToChar (double Lp[][NFLX], double *, double *);

// Explicit template instantiations
template void PrimToChar<IDIR> (double Lp[][NFLX], double *, double *);
template void PrimToChar<JDIR> (double Lp[][NFLX], double *, double *);
template void PrimToChar<KDIR> (double Lp[][NFLX], double *, double *);

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void PrimToChar (double L[][NFLX], double *vp, double *w)
/*!
  Compute the matrix-vector product between the L matrix
  (containing primitive left eigenvectors) and the vector v.
  For efficiency purpose, multiplication is done explicitly, so that
  only nonzero entries of the left primitive eigenvectors are considered.
  
  \param [in]  L     Left eigenvectors
  \param [in]  v     primitive variables
  \param [out] w     characteristic variables
*/
////////////////////////////////////////////////////////////////////////
{
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  w[0] = L[0][VXn]*vp[VXn] + L[0][PRS]*vp[PRS];
  w[1] = L[1][VXn]*vp[VXn] + L[1][PRS]*vp[PRS];
  w[2] = vp[RHO] + L[2][PRS]*vp[PRS];
  w[3] = L[3][VXn]*vp[VXn] + vp[VXt] + L[3][PRS]*vp[PRS];
  w[4] = L[4][VXn]*vp[VXn] + vp[VXb] + L[4][PRS]*vp[PRS];

  // 2. For passive scalars, the characteristic variable is 
  //    equal to the primitive one, since  l = r = (0,..., 1 , 0 ,....)

  #if (NSCL > 0) || RADIATION
  int nv;
  #endif

  #if NSCL > 0
  NSCL_LOOP(nv) w[nv] = vp[nv];
  #endif

  #if RADIATION
  NRAD_LOOP(nv) w[nv] = vp[nv];
  #endif
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void PrimToCons (double *v, double *u, RunConfig &run_config)
/*!
  \param [in]  v           array of primitive variables 
  \param [out] u           array of conservative variables 
  \param [in]  run_config  reference to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
{
  int nv;
  double  rhoh_g2, scrh, g;
  double  beta_fix = 0.9999;
  double  h = Enthalpy (v, run_config);

  g = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];

  if (g >= 1.0){
    g = beta_fix/sqrt(g);
    v[VX1] *= g;
    v[VX2] *= g;
    v[VX3] *= g;

    g = beta_fix*beta_fix;
  }

  g    = 1.0/(1.0 - g);
  scrh = rhoh_g2 = v[RHO]*h*g;
  g    = sqrt(g);

  u[RHO] = v[RHO]*g;
  u[MX1] = scrh*v[VX1];
  u[MX2] = scrh*v[VX2];
  u[MX3] = scrh*v[VX3];
  u[ENG] = rhoh_g2 - v[PRS];

  #if RADIATION
  NRAD_LOOP(nv) u[nv] = v[nv];
  #endif

  #if NSCL > 0 
  NSCL_LOOP(nv) u[nv] = v[nv]*u[RHO];
  #endif
}

#if GEOMETRY != CARTESIAN
////////////////////////////////////////////////////////////////////////
/*!
  Convert covariant conservative to orthonormal
*/
////////////////////////////////////////////////////////////////////////
#if GEOMETRY != CARTESIAN
#pragma acc routine seq
inline void ConsToOrtho (double *u, Metric &metric)
{
  #if GEOMETRY == CYLINDRICAL
  u[MX2] *= metric.h2_inv;
  #if RADIATION
  u[FR2] *= metric.h2_inv;
  #endif
  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2_inv);
  u[MX3] *= metric.h3_inv;
  #if RADIATION
  IF (COVARIANT_MX2, u[FR2] *= metric.h2_inv);
  u[FR3] *= metric.h3_inv;
  #endif
  #endif
}
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Convert orthonormal conservative to covariant (for cons. variables)
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void OrthoToCons (double *u, Metric &metric)
{
  #if GEOMETRY == CYLINDRICAL
  u[MX2] *= metric.h2;
  #if RADIATION
  u[FR2] *= metric.h2;
  #endif
  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2);
  u[MX3] *= metric.h3;
  #if RADIATION
  IF (COVARIANT_MX2, u[FR2] *= metric.h2);
  u[FR3] *= metric.h3;
  #endif
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Convert orthonormal conservative to covariant (for fluxes)
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void OrthoToCons (Ary2D &flux, const int i, Metric &metric)
{
  #if GEOMETRY == CYLINDRICAL
  flux(i,MX2) *= metric.h2;
  #if RADIATION
  flux(i,FR2) *= metric.h2;
  #endif
  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, flux(i,MX2) *= metric.h2);
  flux(i,MX3) *= metric.h3;
  #if RADIATION
  IF (COVARIANT_MX2, flux(i,FR2) *= metric.h2);
  flux(i,FR3) *= metric.h3;
  #endif
  #endif
}
#endif  // GEOMETRY != CARTESIAN 
