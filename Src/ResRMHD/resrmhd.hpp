////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Set labels, indexes and prototypes for the MHD module.

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the MHD module.

  \author M. Rossazza (marco.rossazza@unito.it) 
          A. Mignone  (andrea.mignone@unito.it)
  \date    Oct 23, 2025
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
#define BX1  4
#define BX2  5
#define BX3  6

#define VX1  MX1
#define VX2  MX2
#define VX3  MX3

#if HAVE_ENERGY
  #define ENG  7
  #define PRS  ENG
#endif

#define  EX1   8
#define  EX2   9
#define  EX3  10

#if DIVB_CONTROL == DIV_CLEANING
  #define PSI_GLM  11
#endif

#if DIVE_CONTROL == DIV_CLEANING
  #define  PHI_GLM  12
  #define  CRG      13
#endif

// Define NFluxes and NVAR

#define NFLX (10 +   (DIVB_CONTROL == DIV_CLEANING)  \
                 + 2*(DIVE_CONTROL == DIV_CLEANING) + HAVE_ENERGY)

#ifndef NVAR
  #define NSCL (NTRACER + 0)
  #define NVAR (NFLX + NSCL)
#endif

// Check compatibility 

#if DIVE_CONTROL == DIV_CLEANING && DIVB_CONTROL != DIV_CLEANING
  #error Must set also DIVB_CONTROL == DIV_CLEANING
#endif

#if DIVE_CONTROL == CONSTRAINED_TRANSPORT && DIVB_CONTROL != CONSTRAINED_TRANSPORT
  #error Must set also DIVB_CONTROL == CONSTRAINED_TRANSPORT
#endif

////////////////////////////////////////////////////////////////////////
/* Label the different waves in increasing order
    following the number of vector components.

    IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are
               present only in the MHD-GLM formulation.
               We keep them at the END of the enumeration
               so we can skip them in unnecessary loops.
               Please do NOT change them !                            */
////////////////////////////////////////////////////////////////////////

enum KWAVES {
 KFASTM, KFASTP
 #if HAVE_ENERGY
  , KENTRP
 #endif

 #if DIVB_CONTROL != DIV_CLEANING
  , KDIVB
 #endif

 , KSLOWM, KSLOWP, KALFVM, KALFVP

 #if DIVB_CONTROL == DIV_CLEANING
  , KPSI_GLMM, KPSI_GLMP
 #endif
};

/*! \name Vector Potential Labels
    These may only be used in the STARTUP / INIT  functions.
    They're convenient in obtaining a discretization that preserve
    the divergence-free condition (for staggered field) or if you simply
    wish to initialize the magnetic field from the vector potential.     */
/**@{ */
#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)

#define AX  AX1
#define AY  AX2
#define AZ  AX3
/**@} */

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

 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3

 #define iER    EX1
 #define iEPHI  EX2
 #define iEZ    EX3

#endif

#if GEOMETRY == SPHERICAL

 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3

 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3

 #define iER    EX1
 #define iETH   EX2
 #define iEPHI  EX3

#endif

#ifndef COVARIANT_MX2 
 #define COVARIANT_MX2 NO
#endif

#ifndef COVARIANT_BX2 
 #define COVARIANT_BX2 YES
#endif

#ifndef COVARIANT_EX2 
 #define COVARIANT_EX2 YES
#endif

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping                                              */
////////////////////////////////////////////////////////////////////////

void   ComputeCurrent(double *, double *);
void   ComputeElectricField(double *, double *);

void   IMEXImplicitUpdate (Data *, double, Ary4D &, Ary4D &, Grid *, RBox *);
void   IMEXExplicitSource (Data *, Ary4D & , Grid *);

#pragma acc routine vector
void   Powell (double **, double **, double **, Grid *, double **, int, int);

#pragma omp declare target
#pragma acc routine seq
double Resistive_eta(const State &, const RunConfig &);
#pragma omp end declare target

void MHLL_Solver_T   (Data *, timeStep *, Grid *, RBox *, int);
void MHLLC_Solver_T  (Data *, timeStep *, Grid *, RBox *, int);
void LF_Solver_T     (Data *, timeStep *, Grid *, RBox *, int);

#pragma omp declare target
#pragma acc routine seq
uint32_t  RHD_EntropySolve (double *, double *, RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
uint32_t  RHD_EnergySolve  (double *, double *, RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
uint32_t  RMHD_EnergySolve  (double *, double *, RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
uint32_t  RHD_PressureFix  (double *, double *, RunConfig &);
#pragma omp end declare target

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
inline void Flux (State &, double *, double &, RunConfig &);

// Explicit template instantiations
template void Flux<IDIR>(State &, double *, double &, RunConfig &);
template void Flux<JDIR>(State &, double *, double &, RunConfig &);
template void Flux<KDIR>(State &, double *, double &, RunConfig &);

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void Flux (State &state, double *fx, double &pt, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double w, wg2, Bmag2, Emag2;
  #ifdef CRG
  double eta, eta_1, Ev, vxB;
  #endif
  double g, g2;
  #ifdef GLM_MHD
  double glm_ch = run_config.glm_ch;
  #else
  #endif

  constexpr int MXn = MX1 + dir;

  const double *v = state.v;
  const double *u = state.u;
  const double  h = state.h;

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, EX1, EXn, EXt, EXb);

  g     = u[RHO]/v[RHO];
  g2    = g*g;

  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  Emag2 = v[EX1]*v[EX1] + v[EX2]*v[EX2] + v[EX3]*v[EX3];

  pt  = v[PRS] + 0.5*Bmag2 + 0.5*Emag2;
  w   = v[RHO]*h;
  wg2 = w*g2;

  fx[RHO]  = u[RHO]*v[VXn];
  fx[MX1] = wg2*v[VX1]*v[VXn] - v[BX1]*v[BXn] - v[EX1]*v[EXn];
  fx[MX2] = wg2*v[VX2]*v[VXn] - v[BX2]*v[BXn] - v[EX2]*v[EXn];
  fx[MX3] = wg2*v[VX3]*v[VXn] - v[BX3]*v[BXn] - v[EX3]*v[EXn];
  fx[BXn] =  0.0;
  fx[BXt] = -v[EXb];
  fx[BXb] =  v[EXt];
  fx[EXn] =  0.0;
  fx[EXt] =  v[BXb];
  fx[EXb] = -v[BXt];

  fx[ENG] = u[MXn];
  #if RMHD_REDUCED_ENERGY
  fx[ENG] -= fx[RHO];
  #endif

  fx[MXn] += pt;

  #if DIVB_CONTROL == DIV_CLEANING
  fx[BXn] = v[PSI_GLM];
  fx[PSI_GLM] = glm_ch*glm_ch*v[BXn];
  #endif

  #if DIVE_CONTROL == DIV_CLEANING
  fx[EXn] = v[PHI_GLM];
  fx[PHI_GLM] = glm_ch*glm_ch*v[EXn];

  eta     = state.eta;
  eta_1   = 1.0/eta;
  Ev      = v[EX1]*v[VX1] + v[EX2]*v[VX2] + v[EX3]*v[VX3];
  vxB     = v[VXt]*v[BXb] - v[VXb]*v[BXt];
  fx[CRG] = g*eta_1*(v[EXn] + vxB - Ev*v[VXn]) + u[CRG]*v[VXn];
  #endif
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void  PrimEigenvectors (double *v, double LL[][NFLX], double RR[][NFLX],
                               double *lambda, Grid *grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  printf ("! PrimToChar(): not available\n");

}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void PrimToChar (double Lp[][NFLX], double *v, double *w)
/*! 
*/
////////////////////////////////////////////////////////////////////////
{
  printf ("! PrimToChar(): not available\n");

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
  double  rhoh;
  double  vel2, Bmag2, Emag2;
  double  g, g2;
  double  h = Enthalpy(v, run_config);

  vel2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  g2   = 1.0/(1.0 - vel2);
  g    = sqrt(g2);

  rhoh  = v[RHO]*h;
  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  Emag2 = v[EX1]*v[EX1] + v[EX2]*v[EX2] + v[EX3]*v[EX3];

  u[RHO] = v[RHO]*g;

  u[BX1] = v[BX1];
  u[BX2] = v[BX2];
  u[BX3] = v[BX3];

  u[EX1] = v[EX1];
  u[EX2] = v[EX2];
  u[EX3] = v[EX3];

  #ifdef CRG
  u[CRG] = v[CRG];
  #endif

  u[MX1] = rhoh*g2*v[VX1] + v[EX2]*v[BX3] - v[EX3]*v[BX2];
  u[MX2] = rhoh*g2*v[VX2] + v[EX3]*v[BX1] - v[EX1]*v[BX3];
  u[MX3] = rhoh*g2*v[VX3] + v[EX1]*v[BX2] - v[EX2]*v[BX1];

  #if RMHD_REDUCED_ENERGY == YES
  u[ENG] = rhoh*g2 - v[PRS] - u[RHO] + 0.5*(Emag2 + Bmag2);
  #else
  u[ENG] = rhoh*g2 - v[PRS] + 0.5*(Emag2 + Bmag2);
  #endif

  #ifdef PSI_GLM
  u[PSI_GLM] = v[PSI_GLM];
  #endif
  #ifdef PHI_GLM
  u[PHI_GLM] = v[PHI_GLM];
  #endif

  #if NSCL > 0
  int nv;
  NSCL_LOOP(nv) u[nv] = u[RHO]*v[nv];
  #endif
}

////////////////////////////////////////////////////////////////////////
template<int dir> void MaxSignalSpeed_RHD (double *, double, double *, double *);


#pragma acc routine seq
template<int dir>
inline void MaxSignalSpeed_RHD (double *v, double a2, double *cmin, double *cmax)
/*!
  Compute the maximum and minimum characteristic velocities for the
  MHD equations, cmin= v - cf, cmax = v + cf
  
  \param [in]  state  pointer to a State structure
  \param [out] cmin  1-D array containing the leftmost characteristic
  \param [out] cmin  1-D array containing the rightmost characteristic
  \param [in]  beg   starting index of computation
  \param [in]  end   final   index of computation
*/
////////////////////////////////////////////////////////////////////////
{
  double vx, vt2, vel2, arg, sroot;

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);

  vx   = v[VXn];
  vt2  = v[VXt]*v[VXt] + v[VXb]*v[VXb];
  vel2 = vx*vx + vt2;

  arg  = a2*(1.0 - vx*vx - vt2*a2)*(1.0 - vel2);
  sroot = sqrt(arg);

  if (arg < 0.0) printf ("! MaxSignalSpeed_RHD(): cannot take square root\n");

  *cmax = (vx*(1.0 - a2) + sroot)/(1.0 - vel2*a2);
  *cmin = (vx*(1.0 - a2) - sroot)/(1.0 - vel2*a2);
}


#define ROE_ESTIMATE    NO
#define DAVIS_ESTIMATE  YES
#ifndef CHECK_WAVES
 #define CHECK_WAVES   NO
#endif

////////////////////////////////////////////////////////////////////////
template<int dir> 
void HLL_Speed_RHD (double *, double, double *, double,
                    double*, double*);

#pragma acc routine seq
template<int dir>
inline void HLL_Speed_RHD (double *vL, double a2L, double *vR, double a2R,
                           double* SL, double* SR)
/*!
  Compute leftmost (SL) and rightmost (SR) speeds for the Riemann fan.
  
  \param [in]  stateL   pointer to a state structure for the left state
  \param [in]  stateR   pointer to a state structure for the right state
  \param [out] SL       the (estimated) leftmost speed of the Riemann fan
  \param [out] SR       the (estimated) rightmost speed of the Riemann fan
  \param [in]  beg      starting index of computation
  \param [in]  end      final index of computation
  
  Switches:
  
  ROE_ESTIMATE (YES/NO), DAVIS_ESTIMATE (YES/NO)
  
  These switches set how the wave speed estimates are
  going to be computed. Only one can be set to 'YES', and
  the rest of them must be set to 'NO'
  
  ROE_ESTIMATE:    b_m = \min(0, \min(u_R - c_R, u_L - c_L, u_{roe} - c_{roe}))
  b_m = \min(0, \min(u_R + c_R, u_L + c_L, u_{roe} + c_{roe}))
  
  where u_{roe} and c_{roe} are computed using Roe averages.
  
  DAVIS_ESTIMATE:  b_m = \min(0, \min(u_R - c_R, u_L - c_L))
  b_m = \min(0, \min(u_R + c_R, u_L + c_L))
*/
////////////////////////////////////////////////////////////////////////
{
  double sl_min, sl_max;
  double sr_min, sr_max;

  #if DAVIS_ESTIMATE == YES
  MaxSignalSpeed_RHD<dir> (vL, a2L, &sl_min, &sl_max);
  MaxSignalSpeed_RHD<dir> (vR, a2R, &sr_min, &sr_max);

  *SL = MIN(sl_min, sr_min);
  *SR = MAX(sl_max, sr_max);
  #endif  // DAVIS_ESTIMATE == YES 
}
#undef DAVIS_ESTIMATE
#undef ROE_ESTIMATE

#if GEOMETRY != CARTESIAN
////////////////////////////////////////////////////////////////////////
/*!
  Convert covariant conservative to orthonormal
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void ConsToOrtho (double *u, Metric &metric)
{
  #if GEOMETRY == CYLINDRICAL
  u[MX2] *= metric.h2_inv;

  u[BX2] *= metric.h2;
  u[EX2] *= metric.h2;

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2_inv;)
  u[MX3] *= metric.h3_inv;

  IF (COVARIANT_BX2, u[BX2] *= metric.h2;)
  u[BX3] *= metric.h3;

  IF (COVARIANT_EX2, u[EX2] *= metric.h2;)
  u[EX3] *= metric.h3;
  #endif
}

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

  u[BX2] *= metric.h2_inv;
  u[EX2] *= metric.h2_inv;

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2;)
  u[MX3] *= metric.h3;

  IF (COVARIANT_BX2, u[BX2] *= metric.h2_inv;)
  u[BX3] *= metric.h3_inv;
  IF (COVARIANT_EX2, u[EX2] *= metric.h2_inv;)
  u[EX3] *= metric.h3_inv;
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
  flux(i,BX2) *= metric.h2_inv;
  flux(i,EX2) *= metric.h2_inv;

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, flux(i,MX2) *= metric.h2;)
  flux(i,MX3) *= metric.h3;

  IF (COVARIANT_BX2, flux(i,BX2) *= metric.h2_inv;)
  flux(i,BX3) *= metric.h3_inv;

  IF (COVARIANT_EX2, flux(i,EX2) *= metric.h2_inv;)
  flux(i,EX3) *= metric.h3_inv;
  #endif
}
#endif  // GEOMETRY != CARTESIAN 

#pragma acc routine seq
int QuarticSolve (double, double, double, double, double *);

void Source (Data *, Grid *, double);

#ifndef RMHD_REDUCED_ENERGY
  #define RMHD_REDUCED_ENERGY    NO  /**< By turning RMHD_REDUCED_ENERGY to YES,
                                          we let PLUTO evolve the total energy
                                          minus the mass density contribution. */
#endif

#ifndef USE_NEWTON_BROYDEN   // YES /NO
 #define USE_NEWTON_BROYDEN YES
#endif

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Additional.hppeader files                                         */
////////////////////////////////////////////////////////////////////////

#if DIVB_CONTROL == EIGHT_WAVES
 void POWELL_DIVB_SOURCE(const Sweep *, int, int, Grid *);
 void HLL_DIVB_SOURCE (const Sweep *, double **, int, int, Grid *);
#elif DIVB_CONTROL == DIV_CLEANING
  #include "MHD/GLM/glm.hpp"
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
 #include "MHD/CT/ct.hpp"
 #include "MHD/CT/ct_store.hpp"
#endif
