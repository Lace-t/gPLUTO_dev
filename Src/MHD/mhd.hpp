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
#define BX1  4
#define BX2  5
#define BX3  6
#if HAVE_ENERGY
  #define ENG  7
  #define PRS  ENG
#endif

#define VX1  MX1
#define VX2  MX2
#define VX3  MX3

#if DIVB_CONTROL == DIV_CLEANING
  #define PSI_GLM  (7 + HAVE_ENERGY)
#endif

#define NFLX (7 + HAVE_ENERGY + (DIVB_CONTROL == DIV_CLEANING))

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
  #define iMR    MX1
  #define iBR    BX1

  #define iVPHI  VX2
  #define iBPHI  BX2

  #define iVZ    VX3
  #define iMZ    MX3
  #define iBZ    BX3
#endif

#if GEOMETRY == SPHERICAL
  #define iVR    VX1
  #define iMR    MX1
  #define iBR    BX1
  
  #define iVTH   VX2
  #define iMTH   MX2
  #define iBTH   BX2
 
  #define iVPHI  VX3
  #define iBPHI  BX3
#endif


#ifndef COVARIANT_MX2 
 #define COVARIANT_MX2 NO
#endif
#ifndef COVARIANT_BX2 
 #define COVARIANT_BX2 YES
#endif

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping                                              */
////////////////////////////////////////////////////////////////////////

void ConvertTo4vel (Ary2D &, int);
void GFORCE_Solver_T (Data *, timeStep *, Grid *, RBox *, int);
void HLLC_Solver_T   (Data *, timeStep *, Grid *, RBox *, int);
void HLLD_Solver_T   (Data *, timeStep *, Grid *, RBox *, int);
void HLL_Solver_T    (Data *, timeStep *, Grid *, RBox *, int);
void LF_Solver_T     (Data *, timeStep *, Grid *, RBox *, int);
void Roe_Solver_T    (Data *, timeStep *, Grid *, RBox *, int);

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Include additional header files                                   */
////////////////////////////////////////////////////////////////////////

#if DIVB_CONTROL == DIV_CLEANING
  #include "MHD/GLM/glm.hpp"
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
 #include "MHD/CT/ct.hpp"
 #include "MHD/CT/ct_store.hpp"
#endif
#if RESISTIVITY != NO
 #include "MHD/Resistivity/res.hpp"
#endif

#include "eos.hpp"

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
/*! Compute flux for the MHD equations.
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
template<int dir>
inline void Flux (State &state, double *flux, double &pt, RunConfig &run_config)
{
  constexpr int MXn = MX1 + dir;
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  const double *u = state.u;
  const double *v = state.v;
  #if DIVB_CONTROL == DIV_CLEANING
  double glm_ch = run_config.glm_ch;
  #else
  #endif

  #if HAVE_ENERGY
  double vB = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  pt  = v[PRS] + 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
  #elif EOS == ISOTHERMAL
  pt  = state.a2*v[RHO] + 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
  #else
  #error ! Flux(): not defined for this EoS
  #endif

  // 1a. Compute ideal MHD fluxes

  flux[RHO] = u[MXn];
  flux[MX1] = u[MX1]*v[VXn] - v[BX1]*v[BXn];
  flux[MX2] = u[MX2]*v[VXn] - v[BX2]*v[BXn];
  flux[MX3] = u[MX3]*v[VXn] - v[BX3]*v[BXn];
  flux[BX1] = v[BX1]*v[VXn] - v[BXn]*v[VX1];
  flux[BX2] = v[BX2]*v[VXn] - v[BXn]*v[VX2];
  flux[BX3] = v[BX3]*v[VXn] - v[BXn]*v[VX3];
  #if HAVE_ENERGY
  flux[ENG] = (u[ENG] + pt)*v[VXn] - v[BXn]*vB;
  #endif

  // 1b. Add background field
  //     Note, ptot = p + (B1)^2/2 + B1.B0

  #if BACKGROUND_FIELD == YES
  const double *Bbck = state.Bgf;
  double B0B1 = Bbck[IDIR]*v[BX1] + Bbck[JDIR]*v[BX2] + Bbck[KDIR]*v[BX3];
  pt += B0B1;

  flux[MX1] -= Bbck[BXn-BX1]*v[BX1] + v[BXn]*Bbck[BX1-BX1];
  flux[MX2] -= Bbck[BXn-BX1]*v[BX2] + v[BXn]*Bbck[BX2-BX1];
  flux[MX3] -= Bbck[BXn-BX1]*v[BX3] + v[BXn]*Bbck[BX3-BX1];

  flux[BXt] += v[VXn]*Bbck[BXt-BX1] - Bbck[BXn-BX1]*v[VXt];
  flux[BXb] += v[VXn]*Bbck[BXb-BX1] - Bbck[BXn-BX1]*v[VXb];
  #if HAVE_ENERGY
  flux[ENG] += B0B1*v[VXn] - Bbck[BXn-BX1]*vB;
  #endif
  #endif

  flux[MXn] += pt;

  #ifdef GLM_MHD
  flux[BXn]     = v[PSI_GLM];
  flux[PSI_GLM] = glm_ch*glm_ch*v[BXn];
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

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
template<int dir>
inline void MaxSignalSpeed (State &state, double &cmin, double &cmax)
{
  constexpr int VXn = VX1+dir;

  double gpr, cf, Bmag2, Btmag2;
  double b1, b2, b3;
  const double *v = state.v;
  const double a2 = state.a2;

  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);

  b1 = b2 = b3 = 0.0;

  gpr = a2*v[RHO];

  b1 = v[BXn];
  b2 = v[BXt];
  b3 = v[BXb];

  #if BACKGROUND_FIELD == YES
  const double *bgf = state.Bgf;
  b1 += bgf[BXn - BX1];
  b2 += bgf[BXt - BX1];
  b3 += bgf[BXb - BX1];
  #endif

  Btmag2 = b2*b2 + b3*b3;
  Bmag2  = b1*b1 + Btmag2;

  cf = gpr - Bmag2;
  cf = gpr + Bmag2 + sqrt(cf*cf + 4.0*gpr*Btmag2);
  cf = sqrt(0.5*cf/v[RHO]);

  cmax = v[VXn] + cf;
  cmin = v[VXn] - cf;
}
#pragma omp end declare target

// Function prototype 
template<int dir> 
inline void HLL_Speed (State &, State &, double &, double &);

// Explicit template instantiations
template void HLL_Speed<IDIR> (State &, State &, double &, double &);
template void HLL_Speed<JDIR> (State &, State &, double &, double &);
template void HLL_Speed<KDIR> (State &, State &, double &, double &);

////////////////////////////////////////////////////////////////////////
/*!
  Compute leftmost (SL) and rightmost (SR) speeds for the Riemann fan.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void HLL_Speed (State &stateL, State &stateR, double &SL, double &SR)
{
  double sl_min, sl_max;
  double sr_min, sr_max;

  // Compute speeds using DAVIS Estimate

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
/*
 * Provide left and right eigenvectors and corresponding 
 * eigenvalues for the PRIMITIVE form of the MHD equations 
 * (adiabatic & isothermal cases).
 *
  \param [in]   state       a State structure
  \param [in]   run_config  reference to a RunConfig structure
  \param [out]  LL          left eigenvectors
  \param [out]  RR          right eigenvectors
  \param [out]  lambda      eigenvalues
 *
 *  Wave names and their order are defined as enumeration constants in 
 *  mod_defs.h. 
 *  Notice that, the characteristic decomposition may differ 
 *  depending on the way div.B is treated.
 * 
 *  Advection modes associated with passive scalars are simple cases 
 *  for which lambda = u (entropy mode) and l = r = (0, ... , 1, 0, ...). 
 *  For this reason they are NOT defined here and must be treated 
 *  seperately.
 *
 * \b References:
 *    - "Notes on the Eigensystem of Magnetohydrodynamics", \n
 *      P.L. Roe, D.S. Balsara,
 *      SIAM Journal on Applied Mathematics, 56, 57 (1996)
 *
 *    - "A solution adaptive upwind scheme for ideal MHD", \n
 *      K. G. Powell, P. L. Roe, and T. J. Linde,
 *      Journal of Computational Physics, 154, 284-309, (1999).
 *
 *    - "A second-order unsplit Godunov scheme for cell-centered MHD:
 *       the CTU-GLM scheme"\n
 *      Mignone \& Tzeferacos, JCP (2010) 229, 2117
 *
 *    - "ATHENA: A new code for astrophysical MHD",
 *      J. Stone, T. Gardiner, 
 *      ApJS, 178, 137 (2008)
 *
 * 
 * The derivation of the isothermal eigenvectors follows the
 * consideration given in roe.c.
*/
////////////////////////////////////////////////////////////////////////
#define  SQRT_HALF  (0.70710678118654752440)
#pragma acc routine seq
template<int dir> 
inline void  PrimEigenvectors (State &state, RunConfig &run_config,
                               double LL[][NFLX], double RR[][NFLX],
                               double *lambda)
{
  int k,nv,l;
  double scrh0, scrh1, scrh2, scrh3;
  double u, a, ca2, cf2, cs2;
  double cs, ca, cf;
  double bx2, by2, bz2, bt2, b2;
  double alpha_f, alpha_s, beta_y, beta_z;
  double rho, sqrt_rho, tau, S;
  double a2;
  double *v = state.v;

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);

  // 1a. Compute relevant quantities

  SoundSpeed2 (state, run_config);
  a2 = state.a2;

  u   = v[VXn];
  rho = v[RHO];
  tau = 1.0/rho;
  sqrt_rho = sqrt(rho);

  bx2 = v[BXn]*v[BXn]*tau;        //  Bx^2/rho 
  by2 = v[BXt]*v[BXt]*tau;        //  By^2/rho 
  bz2 = v[BXb]*v[BXb]*tau;        //  Bz^2/rho 
  bt2 = by2 + bz2;                    //  Bt^2/rho 
  b2  = bx2 + bt2;                   

  scrh1 = a2 - b2;
  scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*bt2);

  // 1b. Obtain fast, slow, Alfven & sound speeds

  cf2 = 0.5*(a2 + b2 + scrh0);
  cs2 = a2*bx2/cf2;
  ca2 = bx2;

  cf = sqrt(cf2);
  cs = sqrt(cs2);
  ca = sqrt(ca2);
  a  = sqrt(a2);

  // 2a. Compute alpha_f, alpha_s. Handle degeneracies 

  if (cf == cs) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  }else{
    alpha_f = a2*(cf2-bx2)/cf2/scrh0; // a2-cs2= cf2-b2, see Roe, Balsara Eq. 4 
    alpha_f = MAX(0.0, alpha_f);
    alpha_s = MAX(0.0, 1.0 - alpha_f);

    alpha_f = sqrt(alpha_f);
    alpha_s = sqrt(alpha_s);
  }

  // 2b. Compute beta_y and beta_z, handle degen.

  scrh0 = sqrt(bt2*rho);

  if (scrh0 > 1.e-9) {
    beta_y = v[BXt] / scrh0;
    beta_z = v[BXb] / scrh0;
  } else {
    beta_z = beta_y = SQRT_HALF;
//    beta_z = beta_y = 0.0;
  }

  S = (v[BXn] >= 0.0 ? 1.0 : -1.0);

  // 3. Define primitive right and left eigenvectors RR and LL
  //    for all of the 8 waves.
  //    Left eigenvectors for fast & slow waves can be defined in 
  //    terms of right eigenvectors.

  for (nv = 0; nv < NFLX; nv++) {
  for (l = 0; l < NFLX; l++){
    LL[nv][l] = 0.0;
    RR[nv][l] = 0.0;
  }}

  // 3a-. Fast wave,  (u - c_f)

  k = KFASTM;
  lambda[k] = u - cf;
  scrh0 = alpha_s*cs*S;
  scrh1 = alpha_s*sqrt_rho*a;
  scrh2 = 0.5 / a2;
  scrh3 = scrh2*tau;

  RR[RHO][k] = rho*alpha_f;
  RR[VXn][k] = -cf*alpha_f;
  RR[VXt][k] = scrh0*beta_y;
  RR[VXb][k] = scrh0*beta_z;

  RR[BXt][k] = scrh1*beta_y;
  RR[BXb][k] = scrh1*beta_z;

  #if HAVE_ENERGY
  double scrh4 = alpha_f*a2*rho;
  RR[PRS][k] = scrh4;
  #endif

  LL[k][RHO] = 0.0;
  #if EOS == ISOTHERMAL
  LL[k][RHO] = 0.5*alpha_f/rho;
  #endif
  LL[k][VXn] = RR[VXn][k]*scrh2;
  LL[k][VXt] = RR[VXt][k]*scrh2;
  LL[k][VXb] = RR[VXb][k]*scrh2;
  LL[k][BXt] = RR[BXt][k]*scrh3;
  LL[k][BXb] = RR[BXb][k]*scrh3;
  #if HAVE_ENERGY
  LL[k][PRS] = alpha_f*scrh3;
  #endif

  // 3a+. Fast wave,  (u + c_f)

  k = KFASTP;
  lambda[k] = u + cf;
  RR[RHO][k] = RR[RHO][KFASTM];
  RR[VXn][k] = -RR[VXn][KFASTM];
  RR[VXt][k] = -RR[VXt][KFASTM];
  RR[VXb][k] = -RR[VXb][KFASTM];
  RR[BXt][k] = RR[BXt][KFASTM];
  RR[BXb][k] = RR[BXb][KFASTM];
  #if HAVE_ENERGY
  RR[PRS][k] = RR[PRS][KFASTM];
  #endif

  #if EOS == ISOTHERMAL
  LL[k][RHO] = LL[KFASTM][RHO];
  #endif
  LL[k][VXn] = -LL[KFASTM][VXn];
  LL[k][VXt] = -LL[KFASTM][VXt];
  LL[k][VXb] = -LL[KFASTM][VXb];
  LL[k][BXt] = LL[KFASTM][BXt];
  LL[k][BXb] = LL[KFASTM][BXb];
  #if HAVE_ENERGY
  LL[k][PRS] = LL[KFASTM][PRS];
  #endif

  // 3b. Entropy wave, (u)

  #if HAVE_ENERGY
  k = KENTRP;
  lambda[k] = u;
  RR[RHO][k] = 1.0;
  LL[k][RHO] = 1.0;
  LL[k][PRS] = - 1.0/a2;
  #endif

  // 3c. Magnetic flux, u

#if DIVB_CONTROL != DIV_CLEANING
  k = KDIVB;
  #if 1 // DIVB_CONTROL == EIGHT_WAVES
  lambda[k] = u;
  RR[BXn][k] = 1.0;
  LL[k][BXn] = 1.0;
  #else
  lambda[k] = 0.0;
  #endif
#endif  //  DIVB_CONTROL != DIV_CLEANING 

  // 3d-. Slow wave (u - c_s)

  k = KSLOWM;
  lambda[k] = u - cs;
  scrh0 = alpha_f*cf*S;
  scrh1 = alpha_f*sqrt_rho*a;

  RR[RHO][k] = v[RHO]*alpha_s;
  RR[VXn][k] = -cs*alpha_s;
  RR[VXt][k] = -scrh0*beta_y;
  RR[VXb][k] = -scrh0*beta_z;
  RR[BXt][k] = -scrh1*beta_y;
  RR[BXb][k] = -scrh1*beta_z;
  #if HAVE_ENERGY
  scrh4 = alpha_s*a2*rho;
  RR[PRS][k] = scrh4;
  #endif

  #if EOS == ISOTHERMAL
  LL[k][RHO] = 0.5*alpha_s/rho;
  #endif
  LL[k][VXn] = RR[VXn][k]*scrh2;
  LL[k][VXt] = RR[VXt][k]*scrh2;
  LL[k][VXb] = RR[VXb][k]*scrh2;
  LL[k][BXt] = RR[BXt][k]*scrh3;
  LL[k][BXb] = RR[BXb][k]*scrh3;
  #if HAVE_ENERGY
  LL[k][PRS] = alpha_s*scrh3;
  #endif

  // 3d+. Slow wave (u + c_s)

  k = KSLOWP;
  lambda[k] = u + cs;

  RR[RHO][k] =  RR[RHO][KSLOWM];
  RR[VXn][k] = -RR[VXn][KSLOWM];
  RR[VXt][k] = -RR[VXt][KSLOWM];
  RR[VXb][k] = -RR[VXb][KSLOWM];
  RR[BXt][k] = RR[BXt][KSLOWM];
  RR[BXb][k] = RR[BXb][KSLOWM];
  #if HAVE_ENERGY
  RR[PRS][k] = scrh4;
  #endif

  #if EOS == ISOTHERMAL
  LL[k][RHO] = LL[KSLOWM][RHO];
  #endif
  LL[k][VXn] = -LL[KSLOWM][VXn];
  LL[k][VXt] = -LL[KSLOWM][VXt];
  LL[k][VXb] = -LL[KSLOWM][VXb];
  LL[k][BXt] = LL[KSLOWM][BXt];
  LL[k][BXb] = LL[KSLOWM][BXb];
  #if HAVE_ENERGY
  LL[k][PRS] = LL[KSLOWM][PRS];
  #endif

  // 3e-. Alfven wave (u - c_a)

  k = KALFVM;
  lambda[k] = u - ca;
  scrh2 = beta_y*SQRT_HALF;
  scrh3 = beta_z*SQRT_HALF;

  RR[VXt][k] = -scrh3;
  RR[VXb][k] =  scrh2;
  RR[BXt][k] = -scrh3*sqrt_rho*S;
  RR[BXb][k] =  scrh2*sqrt_rho*S;
  LL[k][VXt] = RR[VXt][k];
  LL[k][VXb] = RR[VXb][k];
  LL[k][BXt] = RR[BXt][k]*tau;
  LL[k][BXb] = RR[BXb][k]*tau;

  // 3e+ Alfven wave (u + c_a)

  k = KALFVP;
  lambda[k] = u + ca;
  RR[VXt][k] =   RR[VXt][KALFVM];
  RR[VXb][k] =   RR[VXb][KALFVM];
  RR[BXt][k] = - RR[BXt][KALFVM];
  RR[BXb][k] = - RR[BXb][KALFVM];

  LL[k][VXt] =   LL[KALFVM][VXt];
  LL[k][VXb] =   LL[KALFVM][VXb];
  LL[k][BXt] = - LL[KALFVM][BXt];
  LL[k][BXb] = - LL[KALFVM][BXb];

  // -------------------------------------------------------------------
  // HotFix: when B=0 in a 2D plane and only Bz != 0, enforce zero
  // jumps in BXt. This is useful with CharTracing since 2 equal and
  // opposite jumps may accumulate due to round-off.
  // It also perfectly consistent, since no (Bx,By) can be generated.
  // -------------------------------------------------------------------

  #if DIMENSIONS <= 2
  if (v[BXn] == 0 && v[BXt] == 0.0) for (k = 0; k < NFLX; k++) RR[BXt][k] = 0.0;
  #endif

  #ifdef GLM_MHD
  double glm_ch = run_config.glm_ch;

  // 3g-. GLM wave,  -glm_ch

  k = KPSI_GLMM;
  lambda[k] = -glm_ch;
  RR[BXn][k]      =  1.0;
  RR[PSI_GLM][k] = -glm_ch;

  LL[k][BXn]      =  0.5;
  LL[k][PSI_GLM] = -0.5/glm_ch;

  // 3g+. GLM wave,  +glm_ch

  k = KPSI_GLMP;
  lambda[k] = glm_ch;
  RR[BXn][k]      =  1.0;
  RR[PSI_GLM][k] =  glm_ch;

  LL[k][BXn]     =  0.5;
  LL[k][PSI_GLM] =  0.5/glm_ch;

  #endif

  // 4. Verify eigenvectors consistency:
  //
  //    - verify orthonormality, L.R = R.L = I
  //     -2. check that A = L.Lambda.R, where A is
  //     the Jacobian dF/dU
  
#if CHECK_EIGENVECTORS == YES
  {
    int ip,jp,kp;
    double dA;
    static double **A, **ALR, **I;
  
    if (A == NULL){
      A   = ARRAY_2D(NFLX, NFLX, double);
      ALR = ARRAY_2D(NFLX, NFLX, double);
      I   = ARRAY_2D(NFLX, NFLX, double);
    }

    // 4a. Verify orthonormality
  
    for (ip = 0; ip < NFLX; ip++){
      for (jp = 0; jp < NFLX; jp++){
      I[ip][jp] = 0.0;
      for (kp = 0; kp < NFLX; kp++){
        I[ip][jp] += LL[ip][kp]*RR[kp][jp];
      }
    }}

    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
      #if (DIVB_CONTROL == NO) || (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
      if (ip == KDIVB || jp == KDIVB) continue;
      #endif
      if ( (ip == jp && fabs(I[ip][jp]-1.0) > 1.e-8) ||
           (ip != jp && fabs(I[ip][jp])>1.e-8) ) {
        printf ("! PrimEigenvectors: Eigenvectors not orthogonal (dir = %d)\n",
               dir);
        ShowMatrix (I, NFLX, 1.e-9);
        printf ("!   i,j = %d, %d  %12.6e \n",ip,jp,a);
        abort();
      }
    }}

    // 4b. Verify A = L.Lambda.R 

    for (ip = 0; ip < NFLX; ip++){
      for (jp = 0; jp < NFLX; jp++) A[ip][jp] = ALR[ip][jp] = 0.0;
    }

    // 4b1. Construct A (Jacobian of primitive form of the equations)

    #if HAVE_ENERGY
  
     A[RHO][RHO] = v[VXn]; A[RHO][VXn] = v[RHO];
     A[VXn][VXn] = v[VXn]; A[VXn][BXt] =  v[BXt]*tau; A[VXn][BXb] = v[BXb]*tau; A[VXn][PRS] = tau;
     A[VXt][VXt] = v[VXn]; A[VXt][BXt] = -v[BXn]*tau; 
     A[VXb][VXb] = v[VXn]; A[VXb][BXb] = -v[BXn]*tau; 
     A[BXt][VXn] = v[BXt]; A[BXt][VXt] = -v[BXn]; A[BXt][BXn] = -v[VXt]; A[BXt][BXt] = v[VXn];
     A[BXb][VXn] = v[BXb]; A[BXb][VXb] = -v[BXn]; A[BXb][BXn] = -v[VXb]; A[BXb][BXb] = v[VXn];
     A[PRS][VXn] = a2*v[RHO]; A[PRS][PRS] =  v[VXn];
  
    #elif EOS == ISOTHERMAL
  
     A[RHO][RHO] = v[VXn] ; A[RHO][VXn] = v[RHO];
     A[VXn][VXn] = v[VXn] ; A[VXn][BXt] =  v[BXt]*tau; A[VXn][BXb] = v[BXb]*tau; 
     A[VXn][RHO] = tau*a2;
     A[VXt][VXt] = v[VXn] ; A[VXt][BXt] = -v[BXn]*tau; 
     A[VXb][VXb] = v[VXn] ; A[VXb][BXb] = -v[BXn]*tau; 
     A[BXt][VXn] = v[BXt] ; A[BXt][VXt] = -v[BXn]; A[BXt][BXn] = -v[VXt]; A[BXt][BXt] = v[VXn];
     A[BXb][VXn] = v[BXb] ; A[BXb][VXb] = -v[BXn]; A[BXb][BXn] = -v[VXb]; A[BXb][BXb] = v[VXn];
  
    #endif
  
    #ifdef GLM_MHD
    A[BXn][PSI_GLM] = 1.0;
    A[PSI_GLM][BXn] = glm_ch*glm_ch;
    #endif
  
    // 4b2. ALR = L.Lambda.R

    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
      ALR[ip][jp] = 0.0;
      for (kp = 0; kp < NFLX; kp++){
        ALR[ip][jp] += RR[ip][kp]*lambda[kp]*LL[kp][jp];
      }
    }}
  
    // 4b3. Verify

    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
  
      // NOTE: if the standard 7-wave formulation is adopted, the column 
      //       and the row corresponding to B(normal) do not exist. 
      //       However, PLUTO uses a primitive matrix with 8  entries and 
      //       the B(normal) column contains two entries (-vy, -vz) which  
      //       cannot be recovered using a 7x7 system.
  
      if (jp == BXn) continue;
      dA = (ALR[ip][jp] - A[ip][jp]);
      if (fabs(dA) > 1.e-5){
        printf ("! PrimEigenvectors(): inconsistent eigenvectors (dir = %d)\n",
                    dir);
        printf ("! A[%d][%d] = %16.9e, R.Lambda.L[%d][%d] = %16.9e\n",
                  ip,jp, A[ip][jp], ip,jp,ALR[ip][jp]);
        printf ("\n\n A = \n");
        ShowMatrix(A, NFLX, 1.e-8);
        printf ("\n\n R.Lambda.L = \n");
        ShowMatrix(ALR, NFLX, 1.e-8);
        abort();
      }
    }}  
  
  } // End block CHECK_EIGENVECTORS
#endif  // CHECK_EIGENVECTORS == YES

}
#undef SQRT_HALF

// Function prototype
template<int dir> 
inline void PrimToChar (double Lp[][NFLX], double *, double *);

// Explicit template instantiations
template void PrimToChar<IDIR> (double Lp[][NFLX], double *, double *);
template void PrimToChar<JDIR> (double Lp[][NFLX], double *, double *);
template void PrimToChar<KDIR> (double Lp[][NFLX], double *, double *);

////////////////////////////////////////////////////////////////////////
/*!
  Compute the matrix-vector multiplcation between the
  the L matrix (containing primitive left eigenvectors)
  and the vector v. The result containing the characteristic
  variables is stored in the vector w = L.v
  
  For efficiency purpose, multiplication is done
  explicitly, so that only nonzero entries
  of the left primitive eigenvectors are considered.
  
  \param [in]  Lp    Left eigenvectors
  \param [in]  v     primitive variables
  \param [out] w     characteristic variables
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir> 
inline void PrimToChar (double Lp[NFLX][NFLX], double *v, double *w)
{
  int k;
  double wv, wB;
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);

  // 1a. Fast waves

  k  = KFASTM;
  wv = Lp[k][VXn]*v[VXn] + Lp[k][VXt]*v[VXt] + Lp[k][VXb]*v[VXb];
  #if HAVE_ENERGY
  wB = Lp[k][PRS]*v[PRS] + Lp[k][BXt]*v[BXt] + Lp[k][BXb]*v[BXb];
  #elif EOS == ISOTHERMAL
  wB = Lp[k][RHO]*v[RHO] + Lp[k][BXt]*v[BXt] + Lp[k][BXb]*v[BXb];
  #endif
  w[KFASTM] =  wv + wB;
  w[KFASTP] = -wv + wB;

  // 1b. Entropy mode

  #if HAVE_ENERGY
  k = KENTRP;
  w[KENTRP] = Lp[k][RHO]*v[RHO] + Lp[k][PRS]*v[PRS];
  #endif

  #ifndef GLM_MHD
  k = KDIVB;
  #if 1 // (DIVB_CONTROL == EIGHT_WAVES)
  w[KDIVB] = v[BXn];
  #else
  w[KDIVB] = 0.0;
  #endif
  #endif  // ndef GLM_MHD 

  // 1c. Slow waves

  k = KSLOWM;
  wv = Lp[k][VXn]*v[VXn] + Lp[k][VXt]*v[VXt] + Lp[k][VXb]*v[VXb];

  #if HAVE_ENERGY
  wB = Lp[k][PRS]*v[PRS] + Lp[k][BXt]*v[BXt] + Lp[k][BXb]*v[BXb];
  #elif EOS == ISOTHERMAL
  wB = Lp[k][RHO]*v[RHO] + Lp[k][BXt]*v[BXt] + Lp[k][BXb]*v[BXb];
  #endif

  w[KSLOWM] =  wv + wB;
  w[KSLOWP] = -wv + wB;

  // 1d. Alfven waves

  k = KALFVM;
  wv = Lp[k][VXt]*v[VXt] + Lp[k][VXb]*v[VXb];
  wB = Lp[k][BXt]*v[BXt] + Lp[k][BXb]*v[BXb];
  w[KALFVM] = wv + wB;
  w[KALFVP] = wv - wB;

  // 1e. GLM waves

  #ifdef GLM_MHD
  k = KPSI_GLMP;
  w[KPSI_GLMM] = Lp[k][BXn]*v[BXn] - Lp[k][PSI_GLM]*v[PSI_GLM];
  w[KPSI_GLMP] = Lp[k][BXn]*v[BXn] + Lp[k][PSI_GLM]*v[PSI_GLM];
  #endif

  // 2. For passive scalars, the characteristic variable is 
  //    equal to the primitive one, since  l = r = (0,..., 1 , 0 ,....)

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

////////////////////////////////////////////////////////////////////////
/*!
  Convert primitive to conservative variables.
  
  \param [in]  v           array of primitive variables 
  \param [out] u           array of conservative variables 
  \param [in]  run_config  reference to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void PrimToCons (double *v, double *u, RunConfig &run_config)
{
  double kinb2;

  #if EOS == IDEAL
  double gmm1 = run_config.gamma - 1.0;
  #else
  #endif

  u[RHO] = v[RHO];
  u[MX1] = v[RHO]*v[VX1];
  u[MX2] = v[RHO]*v[VX2];
  u[MX3] = v[RHO]*v[VX3];

  u[BX1] = v[BX1];
  u[BX2] = v[BX2];
  u[BX3] = v[BX3];

  kinb2   = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  kinb2   = v[RHO]*kinb2  + v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  kinb2  *= 0.5;

  #if EOS == IDEAL
  u[ENG] = kinb2 + v[PRS]/gmm1;
  #endif
  #ifdef GLM_MHD
  u[PSI_GLM] = v[PSI_GLM];
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
  #if RADIATION
  u[FR2] *= metric.h2_inv;
  #endif
  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, u[MX2] *= metric.h2_inv);
  u[MX3] *= metric.h3_inv;

  IF (COVARIANT_BX2, u[BX2] *= metric.h2);
  u[BX3] *= metric.h3;
  #if RADIATION
  IF (COVARIANT_MX2, u[FR2] *= metric.h2_inv);
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
   u[MX2] *= metric.h2;
   u[BX2] *= metric.h2_inv;
  #if RADIATION
  u[FR2] *= metric.h2;
  #endif

  #elif GEOMETRY == SPHERICAL
    IF (COVARIANT_MX2, u[MX2] *= metric.h2);
    u[MX3] *= metric.h3;

    IF (COVARIANT_BX2, u[BX2] *= metric.h2_inv);
    u[BX3] *= metric.h3_inv;
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
  flux(i,BX2) *= metric.h2_inv;
  #if RADIATION
  flux(i,FR2) *= metric.h2;
  #endif

  #elif GEOMETRY == SPHERICAL
  IF (COVARIANT_MX2, flux(i,MX2) *= metric.h2);
  flux(i,MX3) *= metric.h3;

  IF (COVARIANT_BX2, flux(i,BX2) *= metric.h2_inv);
  flux(i,BX3) *= metric.h3_inv;
  #if RADIATION
  IF (COVARIANT_MX2, flux(i,FR2) *= metric.h2);
  flux(i,FR3) *= metric.h3;
  #endif
  #endif
}
#endif  // GEOMETRY != CARTESIAN 
