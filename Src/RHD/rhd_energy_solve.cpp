////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Inversion scheme for RHD using total energy density.

  Try to recover primitive variables by solving for the
  gas pressure using total energy density.
  We employ the algorithm outlined in Section 2 of Mignone, Plewa Bodo (2005).
  Specifically, we solve Eq. (12) using a Newton-Raphson root finder.


  \b References
     - "The Piecewise Parabolic Method for Multidimensional Relativistic
        Fluid Dynamics"\n
        Mignone, Plewa \& Bodo, ApJS (2005) 160, 199.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Feb 11, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define MAX_ITER 20
////////////////////////////////////////////////////////////////////////
uint32_t RHD_EnergySolve (double *u, double *v, RunConfig &run_config)
/*!
  
  \param [in,out]  u      array of conservative variables (entropy will
  be redefined)
  \param [out]     v      array of primitive variables
  
  \return Error codes:
  - 0 = success
  - 1 = solution does not exist ERR_ENG_SOLVE_NO_SOLUTION
  - 2 = negative pressure ERR_NEGATIVE_PRS
  - 3 = inaccurate solution (debug only) ERR_ENG_SOLVE_INACCURATE_SOL
  - 4 = NaN ERR_ENG_SOLVE_PRS_IS_NAN
*/
////////////////////////////////////////////////////////////////////////
{
  int      iter;
  uint32_t err_RHD = 0;
  double   p, D_1, alpha, alpha2, lor2, lor, m, tol=1.e-11;
  double   tau, theta, h, dh_dp, dh_dtau;
  double   yp, dyp, dp, scrh;
  double   D, E, m2, Q;

  D  = u[RHO];
  E  = u[ENG];
  m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];

  Q = E - sqrt(m2 + D*D);

  if (Q < 0.0) {
    return err_RHD |= ERR_ENG_SOLVE_NO_SOLUTION; // Equation does not admit a solution
  }

  #if EOS == IDEAL
  double gmmr = run_config.gamma/(run_config.gamma - 1.0);
  #endif
  m    = sqrt(m2);
  p    = m - E;
  p    = MAX(p, 0.0);
  D_1  = 1.0/D;

  double eps2 = 1.e-12;  // Maximum 1/gamma^2
  double pmin = sqrt(m2/(1.0 - eps2)) - E;

  // --------------------------------------------
  // 1. Solve f(p) = 0 by Newton's method
  // --------------------------------------------

  p = MAX(p, pmin);

  for (iter = 0; iter < MAX_ITER; iter++) {

    alpha  = E + p;
    alpha2 = alpha*alpha;
    lor2   = 1.0 - m2/alpha2;
    lor2   = 1.0/lor2;

    if (lor2 < 1.0){
      return err_RHD |= ERR_ENG_SOLVE_NO_SOLUTION;
    }
    lor  = sqrt(lor2);

    tau   = lor*D_1;
    theta = p*tau;

    #if EOS == IDEAL
    h       = 1.0 + gmmr*theta;
    dh_dp   = gmmr*tau;
    dh_dtau = gmmr*p;
    #elif EOS == TAUB
    h       = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
    scrh    = (5.0*h - 8.0*theta)/(2.0*h - 5.0*theta);
    dh_dp   = tau*scrh;
    dh_dtau = p*scrh;
    #endif

    yp  = D*h*lor - E - p;
    dyp = D*lor*dh_dp - m2*lor2*lor/(alpha2*alpha)*(lor*dh_dtau + D*h) - 1.0;
    dp  = yp/dyp;
    p  -= dp;

    if (p < pmin) p = pmin;
    if (fabs (dp) < tol*E) break;
  }

  // 2. Check if solution is consistent

  if (p < 0.0)          return err_RHD |= ERR_NEGATIVE_PRS;
  if (isnan(p))         return err_RHD |= ERR_NAN;
  if (iter == MAX_ITER) return err_RHD |= ERR_TOO_MANY_ITERATIONS;
  if (   fabs(yp/(E+p)) > 1.e-4 
      || p < (m-E))     return err_RHD |= ERR_ENG_SOLVE_INACCURATE;

  // 3. Solution has been found, update to
  //    converged value of p.

  alpha  = E + p;
  alpha2 = alpha*alpha;
  lor2   = alpha2/(alpha2 - m2);
  lor    = sqrt(lor2);
  tau    = lor*D_1;
  theta  = p*tau;

  #if EOS == IDEAL
  h = 1.0 + gmmr*theta;
  #elif EOS == TAUB
  h = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
  #endif

  v[RHO] = 1.0/tau;
  v[PRS] = p;

  scrh   = 1.0/(u[ENG] + p);  // = 1 / W 

  v[VX1] = u[MX1]*scrh;
  v[VX2] = u[MX2]*scrh;
  v[VX3] = u[MX3]*scrh;

  // 4. Recompute entropy consistently

  #if ENTROPY_SWITCH
    double rho = v[RHO];
    #if EOS == IDEAL
    double gamma = run_config.gamma;
    u[ENTR] = p*lor/pow(rho,gamma-1);
    #elif EOS == TAUB
    theta = p/rho;
    u[ENTR] = p*lor/pow(rho,2.0/3.0)*(1.5*theta + sqrt(2.25*theta*theta + 1.0));
    #endif
  #endif


  return 0; // -- success -- 
}

#undef MAX_ITER
