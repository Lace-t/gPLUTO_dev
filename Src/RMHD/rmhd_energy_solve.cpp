////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Inversion scheme for RMHD using total energy density.

  Try to recover gas pressure from conserved variables {D, m, E, B}
  using the algorithm outlined in Section A3 of Mignone \& McKinney (2007).
  Specifically, we solve Eq. (A4) or (A6) (depending on the value of
  \c ::RMHD_REDUCED_ENERGY) using a Newton-Raphson scheme.
  Here W = rho*h*lorentz^2, E, D, etc... have the same meaning
  as in the mentioned paper.

  \author A. Mignone (andrea.mignone@unito.it)\n
          S. Truzzi  (stefano.truzzi@unito.it)
  \date   Aug 21, 2025

  \b References
     - "Equation of sweep in relativistic magnetohydrodynamics: variable versus
        constant adiabatic index"\n
        Mignone \& Mc Kinney, MNRAS (2007) 378, 1118.

*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#pragma omp declare target
#pragma acc routine seq
inline int gpu_isnan(double x) {
    return isnan(x);
}
#pragma omp end declare target

#define MAX_ITER 50
////////////////////////////////////////////////////////////////////////
uint32_t RMHD_EnergySolve (double *u, double *v, RunConfig &run_config)
/*!
  
  Solve f(W) = 0, where f(W) is Equation (A4) or (A6).
  
  \return Error codes:
  - 0 = success
  - 1 ERR_NEGATIVE_ENG = negative energy ERR_NEGATIVE_ENG
  - 2 ERR_NEGATIVE_PRS = negative pressure
  - 3 ERR_ENG_SOLVE_INACCURATE = inaccurate solution (debug only)
  - 4 ERR_ENG_SOLVE_PRS_IS_NAN = NaN
  - 5 ERR_ENG_SOLVE_NO_SOLUTION = vel2 > 1 during iteration process (solution does not exist ?)*/
////////////////////////////////////////////////////////////////////////
{
  int  k;
  uint32_t error_RMHD = 0;
  int  done;
  double Y1, Y2, scrh;
  double rho;
  double W, W2, S2_W2, fW, dfW, dW;
  double dv2_dW, chi;
  double p, lor, lor2;
  double dp, dp_drho, dp_dchi;
  double dchi_dW, drho_dW;
  #if EOS == IDEAL
  double gamma = run_config.gamma;
  #elif EOS == TAUB
  double chi_p_rho;
  #endif
  #if PHYSICS == RMHD
  double acc = 1.e-11;
  #elif PHYSICS == ResRMHD
  double acc = 1.e-3;   // Needed for guess in IMEX scheme 
  #endif
  double one_m_v2, vel2, Wp;

  double D  = u[RHO];
  double E  = u[ENG];
  double S  = u[MX1]*u[BX1] + u[MX2]*u[BX2] + u[MX3]*u[BX3];
  double m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
  double B2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];
  double S2 = S*S;

  // 0. Define the input parameters of the
  //    parameter structure

  if (E < 0.0) return error_RMHD |= ERR_NEGATIVE_ENG;

  // 1. Provide initial guess by taking the
  //    positive root of Eq. (A27).

  #if RMHD_REDUCED_ENERGY == YES
  Y1 = -4.0*(E + D - B2);
  Y2 = m2 - 2.0*(E + D)*B2 + B2*B2;
  #else
  Y1 = -4.0*(E - B2);
  Y2 = m2 - 2.0*E*B2 + B2*B2;
  #endif
  chi = Y1*Y1 - 12.0*Y2;
  chi = MAX(0.0,chi);
  W = ( - Y1 + sqrt(chi))/6.0;
  W = MAX(D, W);

  #if RMHD_REDUCED_ENERGY == YES
  Wp = W - D;
  #endif

  // 2. Start main iterative loop using Newton-Raphson.

  done = 0; p = -1.0;
  for (k = 1; k < MAX_ITER; k++) {

    #if RMHD_REDUCED_ENERGY == YES
    W = Wp + D;
    #endif

    W2    = W*W;
    S2_W2 = S2/W2;
    Y1    = 1.0/(W + B2);
    Y2    = Y1*Y1;

    vel2 = S2_W2*Y1*(Y1*W + 1.0) + m2*Y2;       // Eq (A3) 
    one_m_v2 = 1.0 - vel2;
    lor2 = 1.0/one_m_v2;
    lor  = sqrt(lor2);
    if (vel2 > 1.0){
      return error_RMHD |= ERR_ENG_SOLVE_NO_SOLUTION;
    }
    #if RMHD_REDUCED_ENERGY == YES
    chi = Wp/lor2 - D*vel2/(lor + 1.0);
    #else
    chi = (W - D*lor)*one_m_v2;
    #endif

    dv2_dW  = -2.0*Y2*(3.0*S2_W2 + Y1*(S2_W2*B2*B2/W + m2)); // Eq (A16) 

    // if chi < 0 we let it converge anyway 
    // kinematical terms 

    dchi_dW =  one_m_v2 - 0.5*lor*(D + 2.0*chi*lor)*dv2_dW;
    drho_dW = -0.5*D*lor*dv2_dW;

    #if EOS == IDEAL

    dp_dchi = (gamma - 1.0)/gamma;    // Eq. (A 18) 
    dp_drho = 0.0;
    p       = chi*dp_dchi;

    #elif EOS == TAUB
    rho = D/lor;

    chi_p_rho = chi + rho;
    scrh = sqrt(9.0*chi*chi + 18.0*rho*chi + 25.0*rho*rho);
    p    = 2.0*chi*(chi_p_rho + rho)/(5.0*chi_p_rho + scrh);  // Eq (A22) 

    scrh    = 1.0/(5.0*chi_p_rho - 8.0*p);
    dp_dchi = (2.0*chi_p_rho - 5.0*p)*scrh;    // Eq (A20) 
    dp_drho = (2.0*chi - 5.0*p)*scrh;          // Eq (A21) 

    #endif

    if (done) break;

    dp  = dp_dchi*dchi_dW + dp_drho*drho_dW;
    #if RMHD_REDUCED_ENERGY == YES
    fW  = Wp + 0.5*(B2 + (B2*m2 - S2)*Y2) - (E + p);   // Eq. (A25) 
    dfW = 1.0 - dp - (B2*m2 - S2)*Y2*Y1;               // Eq. (A8)  
    dW  = fW/dfW;
    Wp -= dW;
    if (fabs(dW) < acc*Wp || fabs(fW) < acc) done = 1;
    #else
    fW  = W + 0.5*(B2 + (B2*m2 - S2)*Y2) - (E + p);  // Eq. (A25) 
    dfW = 1.0 - dp - (B2*m2 - S2)*Y2*Y1;             // Eq. (A8) 
    dW  = fW/dfW;
    W -= dW;
    if (fabs(dW) < acc*W || fabs(fW) < acc) done = 1;
    #endif

  }

  if (p < 0.0)       return error_RMHD |= ERR_NEGATIVE_PRS;
  if (k == MAX_ITER) return error_RMHD |= ERR_TOO_MANY_ITERATIONS;
  if (gpu_isnan(p))  return error_RMHD |= ERR_NAN;

  // 3. Assign output values

  #if RMHD_REDUCED_ENERGY == YES
  W = Wp + D;
  #endif

  v[RHO] = rho = u[RHO]/lor;
  v[PRS] = p;

  Y1   = 1.0/(W + B2);
  scrh = S/W;
  v[VX1] = Y1*(u[MX1] + scrh*u[BX1]);
  v[VX2] = Y1*(u[MX2] + scrh*u[BX2]);
  v[VX3] = Y1*(u[MX3] + scrh*u[BX3]);

  // 4. Recompute entropy consistently

  #if ENTROPY_SWITCH
  {
    #if EOS == IDEAL
    u[ENTR] = p*lor/pow(rho,gamma-1.0);
    #elif EOS == TAUB
    double th = p/rho;
    u[ENTR] = p*lor/pow(rho,2.0/3.0)*(1.5*th + sqrt(2.25*th*th + 1.0));
    #endif
  }
  #endif

  
  return(error_RMHD);  // -- normal exit -- 
}
#undef MAX_ITER

