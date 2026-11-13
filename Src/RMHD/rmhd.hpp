////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Set labels, indexes and prototypes for the RMHD module.

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
#define ENG  7

#define VX1  MX1
#define VX2  MX2
#define VX3  MX3
#define PRS  ENG

#if DIVB_CONTROL == DIV_CLEANING
  #define PSI_GLM  8
  #define DIV_COMP 1
#else
  #define DIV_COMP 0
#endif

#define NFLX (8 + DIV_COMP)

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

#endif


#ifndef COVARIANT_MX2 
 #define COVARIANT_MX2 NO
#endif
#ifndef COVARIANT_BX2 
 #define COVARIANT_BX2 YES
#endif

////////////////////////////////////////////////////////////////////////
/* Module-specific symbolic constants (switches)                      */
////////////////////////////////////////////////////////////////////////
   
#ifndef RMHD_FAST_EIGENVALUES
  #define RMHD_FAST_EIGENVALUES  NO /**< If set to YES, use approximate (and
                                        //  faster) expressions when computing the
                                        //  fast magnetosonic speed, see Sect. 3.3
                                        //  of Del Zanna,  A&A (2006), 473.
                                        //  Solutions of quartic equation is avoided
                                        //  and replace with upper bounds provided by
                                        //  quadratic equation. */
#endif

#ifndef RMHD_REDUCED_ENERGY
  #define RMHD_REDUCED_ENERGY    YES  /**< By turning RMHD_REDUCED_ENERGY to YES,
                                      //    we let PLUTO evolve the total energy
                                      //    minus the mass density contribution. */
#endif

// NR enables primitive recovery and the interface limiter; FOFC is optional.
// NR and FOFC both require stage admission; criterion has no separate switch.
#ifndef PCP_NR
#define PCP_NR NO
#endif

#ifndef PCP_FOFC
#define PCP_FOFC NO
#endif

#ifndef PCP_PRESSURE_FIX
#define PCP_PRESSURE_FIX YES
#endif

// Energy admission distance: factor * DBL_EPSILON * max(1, abs(E)).
// Kept separate from the checker cancellation band and density floor.
#ifndef PCP_ADMISSIBILITY_EPS_FACTOR
#define PCP_ADMISSIBILITY_EPS_FACTOR 256.0
#endif

// Relative retreat toward the fixed-B anchor after a limiter repair.
#ifndef PCP_LIMIT_BACKOFF
#define PCP_LIMIT_BACKOFF 1.e-6
#endif

#ifndef PCP_STAGE_LOCAL_PASSES
#define PCP_STAGE_LOCAL_PASSES 4
#endif

#ifndef CT_RMHD_CONSISTENT
#define CT_RMHD_CONSISTENT NO
#endif

#define FLAG_PCP_FOFC (1u << 29)
#define FLAG_PCP_CT_FAIL (1u << 30)
#if PCP_STAGE_LOCAL_PASSES < 1
#error "PCP retry limits must be nonnegative (local passes positive)"
#endif
#if (PCP_NR == YES) && ((EOS != IDEAL) || ENTROPY_SWITCH || RADIATION)
#error "PCP_NR requires ideal EOS without entropy switch or radiation"
#endif
#if (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)

#if (TIME_STEPPING != RK2 && TIME_STEPPING != RK3) || defined(HIGH_ORDER) || RADIATION || defined(FARGO) || RING_AVERAGE > 1
#error "PCP RunStage supports FV RK2/RK3 without radiation/FARGO/ring average"
#endif

// Stage-local retries leave particle/module scheduling in rk_step. Whole-step
// retries require additional state beyond U/V/B for these coupled modules.
#if FALLBACK_MODE == FALLBACK_RECOMPUTE && ((PARTICLES != NO) || (COOLING != NO))
#error "Enhanced RMHD FALLBACK_RECOMPUTE needs particle/cooling rollback"
#endif

#if FALLBACK_MODE == FALLBACK_BLENDING
#error "PCP stage recovery requires FALLBACK_MODE=NO or FALLBACK_RECOMPUTE"
#endif

#endif
#if CT_RMHD_CONSISTENT && (EOS != IDEAL || ENTROPY_SWITCH || RADIATION)
#error "CT_RMHD_CONSISTENT requires ideal EOS without entropy/radiation"
#endif

#pragma omp declare target
#pragma acc routine seq
int RMHD_PCP_Check(const double *);
#pragma acc routine seq
int RMHD_PCP_Limit(const double *, double *);
#pragma acc routine seq
int RMHD_PCP_NR(const double *, double *, RunConfig &);
#pragma acc routine seq
void RMHD_PCP_SetFallback(const double *, double *, RunConfig &);
#pragma acc routine seq
int RMHD_PCP_Recover(double *, double *, RunConfig &, int);
#pragma acc routine seq
int CT_RMHD_CorrectCell(double *, double *, const double *, RunConfig &);
#pragma omp end declare target
void RMHD_PCP_LimitFaces(Data *, Grid *, RBox *, int);
#if PCP_FOFC == YES
// FOFC specialization of the shared LF kernel in tvdlf.cpp.
void PCP_LF_Solver(Data *, timeStep *, Grid *, RBox *, int);
#endif
void RMHD_PCP_SetStageTime(Data *, timeStep *, int);
uint32_t RMHD_PCP_FinishStage(Data *, Grid *, RBox *);
uint32_t RMHD_PCP_RunStage(Data *, timeStep *, double, const Ary4D &,
                           Ary3D (&)[6], Ary4D &, Ary3D (&)[6], double, double,
                           RBox *, Grid *);

#pragma omp declare target
#pragma acc routine seq
uint32_t  RMHD_EntropySolve (double *, double *, RunConfig &);
#pragma omp end declare target
#pragma omp declare target
#pragma acc routine seq
uint32_t  RMHD_EnergySolve  (double *, double *, RunConfig &);
#pragma omp end declare target
#pragma omp declare target
#pragma acc routine seq
uint32_t  RMHD_PressureFix  (double *, double *, RunConfig &);
#pragma omp end declare target

void QuarticPrintCoeffs(double b, double c, double data, double e);

#pragma omp declare target
#pragma acc routine vector
void Powell (double **, double **, double **, Grid *, double **, int, int);
#pragma omp end declare target

void HLL_Solver_T    (Data *, timeStep *, Grid *, RBox *, int);
void HLLC_Solver_T   (Data *, timeStep *, Grid *, RBox *, int);
void HLLD_Solver_T   (Data *, timeStep *, Grid *, RBox *, int);
void LF_Solver_T     (Data *, timeStep *, Grid *, RBox *, int);
void GFORCE_Solver_T (Data *, timeStep *, Grid *, RBox *, int);

void Source (Data *, Grid *, double);

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Include additional header files                                    */
////////////////////////////////////////////////////////////////////////

#if DIVB_CONTROL == DIV_CLEANING
  #include "MHD/GLM/glm.hpp"
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
 #include "MHD/CT/ct.hpp"
 #include "MHD/CT/ct_store.hpp"
#endif

#include "eos.hpp"
#include "vel_limiter.hpp"
#if PCP_FOFC == YES
#if DIVB_CONTROL == CONSTRAINED_TRANSPORT && CT_EMF_AVERAGE != ARITHMETIC && CT_EMF_AVERAGE != CT_CONTACT && CT_EMF_AVERAGE != CT_FLUX && CT_EMF_AVERAGE != UCT_HLL
#error "PCP LLF correction supports arithmetic/contact/flux/UCT_HLL CT EMFs"
#endif
#endif

////////////////////////////////////////////////////////////////////////
/* Inline functions                                                    */
////////////////////////////////////////////////////////////////////////

// Function prototype
template<int dir> 
inline void Flux (State &, double *, double &, RunConfig &);

// Explicit template instantiations
template void Flux<IDIR>(State &, double *, double &, RunConfig &);
template void Flux<JDIR>(State &, double *, double &, RunConfig &);
template void Flux<KDIR>(State &, double *, double &, RunConfig &);

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
template<int dir>
inline void Flux (State &state, double *flux, double &pt, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  constexpr int MXn = MX1 + dir;
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  const double *u = state.u;
  const double *v = state.v;
  const double h  = state.h;

  double vB, b2, wt, wtg2, Bmag2;
  double b[4], g, g2, g_2;
  
  g     = u[RHO]/v[RHO];
  g2    = g*g;
  g_2   = 1.0/g2;

  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  vB = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  b[IDIR] = g*(v[BX1]*g_2 + vB*v[VX1]);
  b[JDIR] = g*(v[BX2]*g_2 + vB*v[VX2]);
  b[KDIR] = g*(v[BX3]*g_2 + vB*v[VX3]);
  b2 = Bmag2*g_2 + vB*vB;

  pt   = v[PRS] + 0.5*b2;
  wt   = v[RHO]*h + b2;
  wtg2 = wt*g2;

  flux[RHO] = u[RHO]*v[VXn];
  flux[MX1] = wtg2*v[VX1]*v[VXn] - b[IDIR]*b[dir];
  flux[MX2] = wtg2*v[VX2]*v[VXn] - b[JDIR]*b[dir];
  flux[MX3] = wtg2*v[VX3]*v[VXn] - b[KDIR]*b[dir];
  flux[BXn] = 0.0;
  flux[BXt] = v[VXn]*v[BXt] - v[BXn]*v[VXt];
  flux[BXb] = v[VXn]*v[BXb] - v[BXn]*v[VXb];

  flux[ENG] = u[MXn];
  #if RMHD_REDUCED_ENERGY
  flux[ENG] -= flux[RHO];
  #endif

  flux[MXn] += pt;

  #ifdef GLM_MHD
  double glm_ch = run_config.glm_ch;
  flux[BXn]     = v[PSI_GLM];
  flux[PSI_GLM] = glm_ch*glm_ch*v[BXn];
  #else
  #endif
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
inline void SortArray(double *z, int n)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int    i, j;
  double f;

  for (j = 1; j < n; j++){
    f = z[j];
    i = j - 1;
    while(i >= 0 && z[i] > f) {
       z[i+1] = z[i];
       i--;
    }
    z[i+1] = f;
  }
}

////////////////////////////////////////////////////////////////////////
inline int CubicSolve (double b, double c, double d, double z[])
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double b2;
  double Q, R;
  double sQ, arg, theta, cs, sn, p;
  const double one_3  = 1.0/3.0;
  const double one_27 = 1.0/27.0;
  const double sqrt3  = sqrt(3.0);

  b2 = b*b;

  // ----------------------------------------------
  // the expression for f should be
  // Q = c - b*b/3.0; however, to avoid negative
  // round-off making h > 0.0 or g^2/4 - h < 0.0
  // we let c --> c(1- 1.1e-16)
  // ----------------------------------------------

  Q  = b2*one_3 - c*(1.0 - 1.e-16);    // = 3*Q, with Q given by Eq. [5.6.10] 
  R  = b*(2.0*b2 - 9.0*c)*one_27 + d;  // = 2*R, with R given by Eq. [5.6.10] 

Q = MAX(Q, 0.0);
/*
if (fabs(Q) < 1.e-18){
  printLog ("CubicSolve() very small Q = %12.6e\n",Q);
  QuitPLUTO(" ", 1);
}
if (Q < 0.0){
  printLog ("! CubicSolve(): Q = %8.3 < 0 \n",Q);
  QuitPLUTO(" ", 1);
}
*/

  // -------------------------------------------------------
  // We assume the cubic *always* has 3 real root for
  // which R^2 > Q^3.
  // It follows that Q is always > 0
  // -------------------------------------------------------

  sQ  = sqrt(Q)/sqrt3;
  arg = -1.5*R/(Q*sQ);

//  this is to prevent unpleseant situation
    // where both g and i are close to zero       

  arg = MAX(-1.0, arg);
  arg = MIN( 1.0, arg);

  theta = acos(arg)*one_3;     // Eq. [5.6.11], note that  pi/3 < theta < 0 

  cs = cos(theta);              //   > 0   
  sn = sqrt3*sin(theta);    //   > 0   
  p  = -b*one_3;

  z[0] = -sQ*(cs + sn) + p;
  z[1] = -sQ*(cs - sn) + p;
  z[2] =  2.0*sQ*cs    + p;

  return(0);
}

////////////////////////////////////////////////////////////////////////
inline int QuarticSolve (double b, double c, double d, double e, double *z)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int    ifail;
  double b2, sq;
  double a2, a1, a0, u[4];
  double p, q, r;
  // const  double three_256 = 3.0/256.0;
  // const  double one_64 = 1.0/64.0;
  double sz1, sz2, sz3;

  b2 = b*b;

  // --------------------------------------------------------------
  // 1) Compute cubic coefficients using the method outlined in
  // http://eqworld.ipmnet.ru/En/solutions/ae/ae0108.pdf
  // --------------------------------------------------------------

  p = c - b2*0.375;
  q = d + b2*b*0.125 - b*c*0.5;
  r = e - 3.0*b2*b2/256.0 + b2*c/16.0 - 0.25*b*d;

  a2 = 2.0*p;
  a1 = p*p - 4.0*r;
  a0 = -q*q;

  ifail = CubicSolve(a2, a1, a0, u);
  if (ifail != 0) return 1;

  u[0] = MAX(u[0],0.0);
  u[1] = MAX(u[1],0.0);
  u[2] = MAX(u[2],0.0);

  if (u[0] != u[0] || u[1] != u[1] || u[2] != u[2]) return 1;

  sq  = -0.5*DSIGN(q);
  sz1 = sq*sqrt(u[0]);
  sz2 = sq*sqrt(u[1]);
  sz3 = sq*sqrt(u[2]);

  z[0] = -0.25*b + sz1 + sz2 + sz3;
  z[1] = -0.25*b + sz1 - sz2 - sz3;
  z[2] = -0.25*b - sz1 + sz2 - sz3;
  z[3] = -0.25*b - sz1 - sz2 + sz3;
  SortArray(z,4);
/*
  if (debug_print){
    printLog ("Quartic roots = %f  %f  %f  %f; q = %8.3e\n",z[0],z[1],z[2],z[3],q);
    CheckSolutions(b,c,d,e,z);
  }
*/
  return 0;
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
inline void PrimToChar (double Lp[][NFLX], double *v, double *w, Grid *grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  printf ("! PrimToChar(): not available\n");
}

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void PrimToCons (double *v, double *u, RunConfig &run_config)
/*!
  \param [in]  v           array of primitive variables 
  \param [out] u           array of conservative variables 
  \param [in]  run_config  refernece to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
{
  double  vel2, Bmag2;
  double ex, ey, ez, e2;
  double  vx1, vx2, vx3;
  double  g, g2, wt;
  #if EOS == IDEAL
  double gmmr = run_config.gamma/(run_config.gamma - 1.0);
  #endif
  double h = Enthalpy(v, run_config);

  vx1  = v[VX1];
  vx2  = v[VX2];
  vx3  = v[VX3];
  vel2 = vx1*vx1 + vx2*vx2 + vx3*vx3;
  g2 = 1.0/(1.0 - vel2);
  g  = sqrt(g2);

  // Use B x (v x B) for magnetic momentum.  Expanding this as
  // B^2 v - (v.B) B cancels large terms for flow parallel to a strong B.
  ex = vx2*v[BX3] - vx3*v[BX2];
  ey = vx3*v[BX1] - vx1*v[BX3];
  ez = vx1*v[BX2] - vx2*v[BX1];
  e2 = ex*ex + ey*ey + ez*ez;
  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  wt    = v[RHO]*h*g2;

  // -------------------------------------------------------
  // Convert from primitive (v) to conservative (u)
  // -------------------------------------------------------

  u[RHO] = g*v[RHO];
  u[MX1] = wt*vx1 + v[BX2]*ez - v[BX3]*ey;
  u[MX2] = wt*vx2 + v[BX3]*ex - v[BX1]*ez;
  u[MX3] = wt*vx3 + v[BX1]*ey - v[BX2]*ex;

  u[BX1] = v[BX1];
  u[BX2] = v[BX2];
  u[BX3] = v[BX3];

  #if RMHD_REDUCED_ENERGY == YES
  #if EOS == IDEAL
  u[ENG]  =   v[PRS]*(g2*gmmr - 1.0) + u[RHO]*g2*vel2/(g + 1.0)
            + 0.5*(Bmag2 + e2);
  #elif EOS == TAUB
  wt     =   v[PRS]/v[RHO];
  u[ENG] =   v[PRS]*(g2*2.5 - 1.0)
           + u[RHO]*g2*(2.25*wt*wt + vel2)/(g*sqrt(1.0 + 2.25*wt*wt) + 1.0)
           + 0.5*(Bmag2 + e2);
  #endif
  #else
  u[ENG]  = v[RHO]*h*g2 - v[PRS] + 0.5*(Bmag2 + e2);
  #endif  // RMHD_REDUCED_ENERGY 

  #if RADIATION
  u[ENR] = v[ENR];
  u[FR1] = v[FR1];
  u[FR2] = v[FR2];
  u[FR3] = v[FR3];
  #endif

  #if NSCL > 0
  int nv;
  NSCL_LOOP(nv) u[nv] = u[RHO]*v[nv];
  #endif

  #ifdef GLM_MHD
   u[PSI_GLM] = v[PSI_GLM];
  #endif
}
#pragma omp end declare target

// Function prototype
template<int dir> 
inline void MaxSignalSpeed (State &, double &, double &);

// Explicit template instantiations
template void MaxSignalSpeed<IDIR>(State &, double &, double &);
template void MaxSignalSpeed<JDIR>(State &, double &, double &);
template void MaxSignalSpeed<KDIR>(State &, double &, double &);

#if RMHD_FAST_EIGENVALUES == NO  
////////////////////////////////////////////////////////////////////////
/*!
  Compute the maximum and minimum characteristic velocities for the
  RMHD equations, cmin= v - cf, cmax = v + cf
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void MaxSignalSpeed (State &state, double &cmin, double &cmax)
{
  int    err = 0;

  constexpr int VXn = VX1 + dir;
  constexpr int BXn = BX1 + dir;

  double lambda[NFLX];
  double  vB, vB2, Bm2, w, w_1;
  double  eps2, one_m_eps2, a2_w;
  double  vx, vx2, u02;
  double  a4, a3, a2, a1, a0;
  double  scrh1, scrh2, scrh;
  double  b2, del, z[4];
  double *v  = state.v;
  double cs2 = state.a2;
  double h   = state.h;

  scrh   = fabs(v[VXn])/sqrt(cs2);

  vB  = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  u02 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  Bm2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];

  if (u02 < 1.e-14) {

    // 1. If total velocity is = 0, the eigenvalue equation
    //    reduces to a biquadratic:
    //
    //      x^4 + a1*x^2 + a0 = 0
    //
    //    with 
    //
    //      a0 = cs^2*bx*bx/W > 0,
    //      a1 = - a0 - eps^2 < 0.

    w_1  = 1.0/(v[RHO]*h + Bm2);
    eps2 = cs2 + Bm2*w_1*(1.0 - cs2);
    a0   = cs2*v[BXn]*v[BXn]*w_1;
    a1   = - a0 - eps2;
    del  = a1*a1 - 4.0*a0;
    del  = MAX(del, 0.0);
    del  = sqrt(del);

    lambda[KFASTP] =  sqrt(0.5*(-a1 + del));
    lambda[KFASTM] = -lambda[KFASTP];
    lambda[KSLOWP] =  sqrt(0.5*(-a1 - del));
    lambda[KSLOWM] = -lambda[KSLOWP];

  } else { 

    // 2. Total velocity is != 0

    vB2 = vB*vB;
    u02 = 1.0/(1.0 - u02);
    b2  = Bm2/u02 + vB2;
    w_1 = 1.0/(v[RHO]*h + b2);
    vx  = v[VXn];
    vx2 = vx*vx;

    if (fabs(v[BXn]) < 1.e-14){
  
      // 2a. Zero normal component

      w     = v[RHO]*h;
      scrh1 = b2 + cs2*(w - vB2); // wc_s^2 + Q 
      scrh2 = w*(1.0 - cs2)*u02;
      del   = 4.0*scrh1*(scrh2*(1.0 - vx2) + scrh1);

      a2  = scrh1 + scrh2;
      a1  = -2.0*vx*scrh2;

      lambda[KFASTP] = 0.5*(-a1 + sqrt(del))/a2;
      lambda[KFASTM] = 0.5*(-a1 - sqrt(del))/a2;
      lambda[KSLOWP] = lambda[KSLOWM] = v[VXn];

    }else{

      // 2b. General case

      scrh1 = v[BXn]/u02 + vB*vx;  // -- this is bx/u0 -- 
      scrh2 = scrh1*scrh1;

      a2_w       = cs2*w_1;
      eps2       = (cs2*v[RHO]*h + b2)*w_1;
      one_m_eps2 = u02*v[RHO]*h*(1.0 - cs2)*w_1;

      //  Define coefficients for the quartic

      scrh = 2.0*(a2_w*vB*scrh1 - eps2*vx);
      a4 = one_m_eps2 - a2_w*vB2 + eps2;
      a3 = - 4.0*vx*one_m_eps2 + scrh;
      a2 =   6.0*vx2*one_m_eps2 + a2_w*(vB2 - scrh2) + eps2*(vx2 - 1.0);
      a1 = - 4.0*vx*vx2*one_m_eps2 - scrh;
      a0 = vx2*vx2*one_m_eps2 + a2_w*scrh2 - eps2*vx2;

      if (a4 < 1.e-12){
        err = 1;
      }

      scrh = 1.0/a4;

      a3 *= scrh;
      a2 *= scrh;
      a1 *= scrh;
      a0 *= scrh;

      int iflag = QuarticSolve(a3, a2, a1, a0, z);

      if (iflag){
        err = 1;
      }

      lambda[KFASTM] = z[0];
      lambda[KFASTP] = z[3];
      lambda[KSLOWP] = z[2];
      lambda[KSLOWM] = z[1];
    }
  } // end if (u02 > 1.e-14)
  cmax = lambda[KFASTP];
  cmin = lambda[KFASTM];
}

#elif RMHD_FAST_EIGENVALUES == YES  

////////////////////////////////////////////////////////////////////////
/*!
  Compute an approximate expression for the fast magnetosonic speed
  using the upper-bound estimated outlined by
  Leismann et al. (A&A 2005, 436, 503), Eq. [27].
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void MaxSignalSpeed (State &state, double &cmin, double &cmax)
{
  int    err = 0;

  constexpr int VXn = VX1 + dir;

  double  vel2, vB, Bmag2;
  double  vx, vx2, vB2, b2, ca2;
  double  om2, vl, scrh;
  double  *v  = state.v;
  double  cs2 = state.a2;
  double  h   = state.h;

  vel2  = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  vB    = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];  

  vx  = v[VXn];
  vx2 = vx*vx;
  vB2 = vB*vB;
  b2  = Bmag2*(1.0 - vel2) + vB2;
  ca2 = b2/(v[RHO]*h + b2);

  om2  = cs2 + ca2 - cs2*ca2;
  vl   = vx*(1.0 - om2)/(1.0 - vel2*om2);
  scrh = om2*(1.0 - vel2)*((1.0 - vel2*om2) - vx2*(1.0 - om2));
  scrh = sqrt(scrh)/(1.0 - vel2*om2);

  cmin = vl - scrh;
  cmax = vl + scrh;

  if (fabs(cmin) > 1.0 || fabs(cmax) > 1.0){
    err = 1;
  }
  if (isnan(cmin) || isnan(cmax)){
    err = 1;
  }

}
#endif // RMHD_FAST_EIGENVALUES == YES

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

  // Compute speeds using DAVIS Estimate

  MaxSignalSpeed<dir> (stateL, sl_min, sl_max);
  MaxSignalSpeed<dir> (stateR, sr_min, sr_max);

  SL = MIN(sl_min, sr_min);
  SR = MAX(sl_max, sr_max);
}

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
  u[MX2] *= metric.h2_inv;                    // mphi = U/R
  IF (COVARIANT_BX2, u[BX2] *= metric.h2;)    // Bx2  = U*R
  #if RADIATION
  u[FR2] *= metric.h2_inv;
  #endif

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2_inv;)
  u[MX3] *= metric.h3_inv;

  IF (COVARIANT_BX2, u[BX2] *= metric.h2;)
  u[BX3] *= metric.h3;
  #if RADIATION
  IF (COVARIANT_MX2, u[FR2] *= metric.h2_inv;)
  u[FR3] *= metric.h3_inv;
  #endif
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
  u[MX2] *= metric.h2;                           // U = mphi*R
  IF (COVARIANT_BX2, u[BX2] *= metric.h2_inv;)   // U = Bphi/R
  #if RADIATION
  u[FR2] *= metric.h2;
  #endif

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2;)
  u[MX3] *= metric.h3;

  IF (COVARIANT_BX2, u[BX2] *= metric.h2_inv;)
  u[BX3] *= metric.h3_inv;
  #if RADIATION
  IF (COVARIANT_MX2, u[FR2] *= metric.h2;)
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
  IF (COVARIANT_BX2, flux(i,BX2) *= metric.h2_inv;)
  #if RADIATION
  flux(i,FR2) *= metric.h2;
  #endif

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, flux(i,MX2) *= metric.h2;)
  flux(i,MX3) *= metric.h3;

  IF (COVARIANT_BX2, flux(i,BX2) *= metric.h2_inv;)
  flux(i,BX3) *= metric.h3_inv;
  #if RADIATION
  IF (COVARIANT_MX2, flux(i,FR2) *= metric.h2;)
  flux(i,FR3) *= metric.h3;
  #endif
  #endif
}
#endif // GEOMETRY != CARTESIAN

////////////////////////////////////////////////////////////////////////
// FOFC stage cache. PCP_FOFC owns the interface; the backend owns storage.
////////////////////////////////////////////////////////////////////////
#if PCP_FOFC == YES
namespace fofc_cache {
extern bool active;
extern int pass;
extern Ary4D *rhs_target;

// Unsupported configurations keep the original reconstruction/recovery path.
// These are implementation limits, not additional algorithm switches.
inline constexpr bool reusable() {
#if defined(_OPENACC)
#if GEOMETRY == CARTESIAN && DIVB_CONTROL == DIV_CLEANING && \
    !defined(HIGH_ORDER) && !RADIATION && RESISTIVITY == NO
  return true;
#else
  return false;
#endif
#elif defined(_OPENMP) || defined(HIGH_ORDER)
  return false;
#else
#if !(INTERNAL_BOUNDARY == YES && INTERNAL_BOUNDARY_REFLECT == YES) && \
    !RADIATION && !SELF_GRAV
  return true;
#else
  return false;
#endif
#endif
}

struct Stage {
  Stage();
  ~Stage();
  Stage(const Stage &) = delete;
  Stage &operator=(const Stage &) = delete;
};
struct RHSTarget {
  explicit RHSTarget(Ary4D &u);
  ~RHSTarget();
  RHSTarget(const RHSTarget &) = delete;
  RHSTarget &operator=(const RHSTarget &) = delete;
private:
  Ary4D *previous;
};
bool stream_enabled();
bool rk_enabled();
inline bool restore() { return active && pass > 0; }
inline bool save() { return active && pass == 0; }
void transfer(Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir, bool save);
}

void RightHandSideFused(Data *data, timeStep *Dts, const Grid *const grid,
                       double dt, RBox *box, int dir, const Ary4D &target);
#endif
