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
#if HAVE_ENERGY
  #define ENG  4
  #define PRS  ENG
#endif

#define VX1  MX1
#define VX2  MX2
#define VX3  MX3

#define NFLX (4 + HAVE_ENERGY)

////////////////////////////////////////////////////////////////////////
/* Now define more convenient and user-friendly 
     pointer labels for geometry setting                              */
////////////////////////////////////////////////////////////////////////

#if GEOMETRY == CYLINDRICAL 
  #define iVR    VX1
  #define iMR    MX1

  #define iVPHI  VX2
//  #define iMPHI  MX2

  #define iVZ    VX3
  #define iMZ    MX3
#endif

#if GEOMETRY == SPHERICAL 
  #define iVR    VX1
  #define iMR    MX1

  #define iVTH   VX2
  #define iMTH   MX2

  #define iVPHI  VX3
//  #define iMPHI  MX3
#endif


#ifndef COVARIANT_MX2 
 #define COVARIANT_MX2  NO
#endif

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping                                              */
////////////////////////////////////////////////////////////////////////

void LF_Solver_T   (Data *, timeStep *, Grid *, RBox *, int);
void Roe_Solver_T  (Data *, timeStep *, Grid *, RBox *, int);
void HLL_Solver_T  (Data *, timeStep *, Grid *, RBox *, int);
void HLLC_Solver_T (Data *, timeStep *, Grid *, RBox *, int);
void TwoShock_Solver_T (Data *, timeStep *, Grid *, RBox *, int);

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Include additional header files                                   */
////////////////////////////////////////////////////////////////////////

#include "eos.hpp"

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
/*!
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void Flux (State &state, double *flux, double &p)
{
  constexpr int MXn = MX1 + dir;
  constexpr int VXn = VX1 + dir;
  const double *u    = state.u;
  const double *v    = state.v;

  flux[RHO] = u[MXn];
  flux[MX1] = u[MX1]*v[VXn];
  flux[MX2] = u[MX2]*v[VXn];
  flux[MX3] = u[MX3]*v[VXn];

  #if HAVE_ENERGY
  flux[ENG] = (u[ENG] + v[PRS])*v[VXn];
  p = v[PRS];
  #elif EOS == ISOTHERMAL
  p = state.a2*v[RHO];
  #endif

  flux[MXn] += p;
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
  Compute the maximum signal speed cmax in a given direction.
*/
////////////////////////////////////////////////////////////////////////
{
  constexpr int VXn = VX1 + dir;
  const double a = sqrt(state.a2);

  cmax = state.v[VXn] + a;
  cmin = state.v[VXn] - a;
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
#pragma omp declare target
#pragma acc routine seq
template<int dir> 
inline void  PrimEigenvectors (State &state, RunConfig &run_config,
                               double LL[][NFLX], double RR[][NFLX],
                               double *lambda)
/*!
  Provide left and right eigenvectors and corresponding
  eigenvalues for the primitive form of the HD equations
  (adiabatic, pvte & isothermal cases).
  
  \param [in]   state       a State structure
  \param [in]   run_config  reference to a RunConfig structure
  \param [out]  LL          left eigenvectors
  \param [out]  RR          right eigenvectors
  \param [out]  lambda      eigenvalues
  
  Wave names and their order are defined as enumeration constants in
  mod_defs.h.
  
  Advection modes associated with passive scalars are simple cases
  for which lambda = u (entropy mode) and l = r = (0, ... , 1, 0, ...).
  For this reason they are NOT defined here and must be treated
  seperately.
  
  \b References:
  - "Riemann Solvers and Numerical Methods for Fluid Dynamics", \n
  Toro, 1997 Springer-Verlag, Eq. [3.18]
*/
////////////////////////////////////////////////////////////////////////
{
  int     k,l,nv;
  double  rhocs, rho_cs, cs, cs2;
  double *v = state.v;

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  // ----------------------------------------------
  // 1. Reset eigenvectors
  // ----------------------------------------------

  for (nv = 0; nv < NFLX; nv++){
  for (l = 0; l < NFLX; l++){
    LL[nv][l] = 0.0;
    RR[nv][l] = 0.0;
  }}

  SoundSpeed2 (state, run_config);
  cs2 = state.a2;
  cs = sqrt(cs2);

  rhocs  = v[RHO]*cs;
  rho_cs = v[RHO]/cs;

  for (nv = 0; nv < NFLX; nv++){
  for (l = 0; l < NFLX; l++){
    LL[nv][l] = 0.0;
    RR[nv][l] = 0.0;
  }}

  // ----------------------------------------------
  // 2. Compute RIGHT eigenvectors
  // ----------------------------------------------

  lambda[0]  = v[VXn] - cs;  //  lambda = u - c   
  RR[RHO][0] =  0.5*rho_cs;
  RR[VXn][0] = -0.5;
  #if HAVE_ENERGY
  RR[PRS][0] =  0.5*rhocs;
  #endif

  lambda[1]  = v[VXn] + cs;  //  lambda = u + c   
  RR[RHO][1] = 0.5*rho_cs;
  RR[VXn][1] = 0.5;
  #if HAVE_ENERGY
  RR[PRS][1] = 0.5*rhocs;
  #endif

  for (k = 2; k < NFLX; k++) lambda[k] = v[VXn];  //  lambda = u   

  #if HAVE_ENERGY
  RR[RHO][2] = 1.0;
  RR[VXt][3] = 1.0;
  RR[VXb][4] = 1.0;
  #elif EOS == ISOTHERMAL
  RR[VXt][2] = 1.0;
  RR[VXb][3] = 1.0;
  #endif

  // ----------------------------------------------
  // 3. Compute LEFT eigenvectors
  // ----------------------------------------------

  LL[0][VXn] = -1.0;
  #if HAVE_ENERGY
  LL[0][PRS] =  1.0/rhocs;
  #elif EOS == ISOTHERMAL
  LL[0][RHO] =  1.0/rho_cs;
  #endif

  LL[1][VXn] = 1.0;
  #if HAVE_ENERGY
  LL[1][PRS] = 1.0/rhocs;
  #elif EOS == ISOTHERMAL
  LL[1][RHO] =  1.0/rho_cs;
  #endif

  #if HAVE_ENERGY
  LL[2][RHO] = 1.0;
  LL[3][VXt] = 1.0;
  LL[4][VXb] = 1.0;
  LL[2][PRS] = -1.0/cs2;
  #elif EOS == ISOTHERMAL
  LL[2][VXt] = 1.0;
  LL[3][VXb] = 1.0;
  #endif

}
#pragma omp end declare target


// Function prototype
template<int dir> 
inline void PrimToChar (double Lp[][NFLX], double *, double *);

// Explicit template instantiations
template void PrimToChar<IDIR> (double Lp[][NFLX], double *, double *);
template void PrimToChar<JDIR> (double Lp[][NFLX], double *, double *);
template void PrimToChar<KDIR> (double Lp[][NFLX], double *, double *);

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
template<int dir>
inline void PrimToChar (double L[NFLX][NFLX], double *v, double *w)
/*!
  Compute the matrix-vector multiplcation between the
  the L matrix (containing primitive left eigenvectors)
  and the vector v. The result containing the characteristic
  variables is stored in the vector w = L.v
  
  For efficiency purpose, multiplication is done
  explicitly, so that only nonzero entries
  of the left primitive eigenvectors are considered.
  
  \param [in]  L     Left eigenvectors
  \param [in]  v     primitive variables
  \param [out] w     characteristic variables
*/
////////////////////////////////////////////////////////////////////////
{
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  #if HAVE_ENERGY
  w[0] = L[0][VXn]*v[VXn] + L[0][PRS]*v[PRS];
  w[1] = L[1][VXn]*v[VXn] + L[1][PRS]*v[PRS];
  w[2] = v[RHO] + L[2][PRS]*v[PRS];
  w[3] = v[VXt];
  w[4] = v[VXb];
  #elif EOS == ISOTHERMAL
  w[0] = L[0][RHO]*v[RHO] + L[0][VXn]*v[VXn];
  w[1] = L[1][RHO]*v[RHO] + L[1][VXn]*v[VXn];
  w[2] = v[VXt];
  w[3] = v[VXb];
  #endif

  // ----------------------------------------------
  // For passive scalars, the characteristic
  // variable is equal to the primitive one,
  // since  l = r = (0,..., 1 , 0 ,....)
  // ----------------------------------------------

  #if (NSCL > 0) || RADIATION
  int nv;
  #endif

  #if NSCL > 0
  NSCL_LOOP(nv) w[nv] = v[nv];
  #endif

  #if RADIATION
  NRAD_LOOP(nv) w[nv] = v[nv];
  #endif
  
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void PrimToCons (double *v, double *u, RunConfig &run_config)
/*!
  Convert an array of primitive quantities to conservative
  
  \param [in]  v           array of primitive variables 
  \param [out] u           array of conservative variables 
  \param [in]  run_config  reference to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
{

  u[RHO] = v[RHO];

  u[MX1] = v[RHO]*v[VX1];
  u[MX2] = v[RHO]*v[VX2];
  u[MX3] = v[RHO]*v[VX3];

  #if HAVE_ENERGY
  double kin;
  kin    = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  kin    = v[RHO]*kin*0.5;
  u[ENG] = kin + v[PRS]/(run_config.gamma - 1.0);
  #else
  #endif

  #if (NSCL > 0) || RADIATION
  int nv;
  #endif

  #if NSCL > 0
  NSCL_LOOP(nv) u[nv] = v[RHO]*v[nv];
  #endif

  #if RADIATION
  NRAD_LOOP(nv) u[nv] = v[nv];
  #if IRRADIATION
  u[FIR] = v[FIR];
  #endif
  #endif
}
#pragma omp end declare target

#if GEOMETRY != CARTESIAN
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void ConsToOrtho (double *u, Metric &metric)
/*!
  Convert covariant conservative to orthonormal (v2)
*/
////////////////////////////////////////////////////////////////////////
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
#pragma omp end declare target
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void OrthoToCons (double *u, Metric &metric)
/*!
  Convert orthonormal conservative to covariant (for cons. variables)
*/
////////////////////////////////////////////////////////////////////////
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
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void OrthoToCons (Ary2D &flux, const int i, Metric &metric)
/*!
  Convert orthonormal conservative to covariant (for fluxes)
*/
////////////////////////////////////////////////////////////////////////
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
#pragma omp end declare target
#endif  // GEOMETRY != CARTESIAN 

//#define CENTROID(r, dr)  ( (r) + (dr)*(dr)/(12.0*(r)) )
#define CENTROID(r, dr)  (r)

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Additional.hppeader files                                         */
////////////////////////////////////////////////////////////////////////
