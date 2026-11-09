////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Relativistic shock tube problems.

  Shock tube problems are widely used for code benchmarking since they
  directly probe the solver's ability in resolving elementary waves and their
  propagation.
  In addition, numerical results can be compared with analytical solutions.

  Here we provide a simple suite of shock tube for the relativistic hydro
  module, by following Mignone \& Bodo (2005).
  \f[
     \left(\rho,\, v_x,\, p\right) = \left\{\begin{array}{ll}
       \left(\rho,\, v_x,\, p\right)_L & \quad\mathrm{for}\quad x < 0.5
       \\[0.5em]
       \left(\rho,\, v_x,\, p\right)_R & \quad\mathrm{otherwise}
      \end{array}\right.
  \f]
  The first four configurations correspond to Problem 1-4 of Mignone
  \& Bodo (2005):

  <CENTER>
  Test | rho  | vx  | p   |   Gamma
  -----|------|-----|-----|--------------
   1L  |  1   | 0.9 | 1   |    4/3
   1R  |  1   | 0.0 | 10  |    -
   2L  |  1   |-0.6 | 10  |    5/3
   2R  |  10  | 0.5 | 20  |    -
   3L  |  10  |-0.6 | 40/3|    5/3
   3R  |  1   | 0.0 | 0   |    -
   4L  |  1   | 0.0 | 1.e3|    5/3
   4R  |  1   | 0.0 |1.e-2|    -
  </CENTER>

  Conf #05 and #06 are equivalent to the first two,
  using the entropy equation.
  Results are shown in the four figures below:

  \author A. Mignone (mignone@to.infn.it)
  \date   May 03, 2020

  \b References
     - "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics"
        Mignone \& Bodo, MNRAS (2005) 364, 126
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*
*/ 
////////////////////////////////////////////////////////////////////////
{
  run_config.gamma = run_config.inputParam[GAMMA_EOS];
  if (x < 0.5){
    v[RHO] = run_config.inputParam[DN_L];
    v[VX1] = run_config.inputParam[VX_L];
    v[VX2] = run_config.inputParam[VY_L];
    v[PRS] = run_config.inputParam[PR_L];
  }else{
    v[RHO] = run_config.inputParam[DN_R];
    v[VX1] = run_config.inputParam[VX_R];
    v[VX2] = run_config.inputParam[VY_R];
    v[PRS] = run_config.inputParam[PR_R];
  }

  #if ENTROPY_SWITCH
  v[ENTR] = Entropy(v,run_config);
  #endif

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*!
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
}
