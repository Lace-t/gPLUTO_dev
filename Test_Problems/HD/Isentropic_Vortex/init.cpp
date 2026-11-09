////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Isentropic vortex problem.

  The isentropic vortex is a smooth exact solution of the 2D Euler
  equations.
  It consists of a single vortex centered at \c (5,5) in pressure
  equilibrium:
  \f[
    \Big(\delta v_x\,, \delta v_y\Big)
    = - (y_c,\,x_c)\frac{\epsilon}{2\pi}
                   \exp\left(\frac{1-r_c^2}{2}\right) \,,\quad
    T = \frac{p}{\rho}
      = 1 - \frac{(\Gamma-1)\epsilon^2}{8\Gamma\pi^2}\exp(1-r_c^2)\,,\quad
    s = \frac{p}{\rho^\Gamma}
  \f]
  where \f$ (x_c,\,y_c) = (x-5,\,y-5),\, r_c=\sqrt{x_c^2+y_c^2},\,\epsilon=5\f$
  while \f$\Gamma = 1.4\f$ is the adiabatic index.
  A constant entropy is used with value \c s=1.
  The vortex shifts along the main diagonal of the computational domain
  with uniform velocity \c (1,1) and returns in its original position
  after \c t=10.
  The domain is assumed to be periodic and its size must be taken large
  enough to ensure there is no interaction between off-domain vortices.

  The vortex problem is often used as computational benchmark to test
  the accuracy and dissipation of a numerical scheme in reproducing
  the vortex structure after several revolutions.

  There're no input parameters for this problem.

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CARTESIAN  |2  |LINEAR     |RK3        |hll
   #02 |HD  |CARTESIAN  |2  |LINEAR     |RK2        |roe
   #03 |HD  |CARTESIAN  |2  |PARABOLIC  |RK3        |tvdlf
   #04 |HD  |CARTESIAN  |2  |MP5        |RK3        |hll
   #05 |HD  |CARTESIAN  |2  |MP5        |RK4        |hll
   #06 |HD  |CARTESIAN  |2  |MP5        |RK4        |hll
   #07 |HD  |CARTESIAN  |2  |WENOZ      |RK4        |hll
   #08 |HD  |CARTESIAN  |2  |LINEAR     |RK2        |two_shock
  <\CENTER>

  \author A. Mignone (mignone@ph.unito.it)
          V. Berta   (vittoria.berta@unito.it)

  \date   Oct 3, 2014

  \b References
     - C.W. Shu, "High-order Finite Difference and Finite Volume WENO
       Schemes and Discontinuous Galerkin Methods for CFD",
       ICASE Report No.2001-11,  NASA/CR-2001-210865
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
 *
 *
*/
////////////////////////////////////////////////////////////////////////
{
  double s, T, r2, xt, yt, eps;
  double k0;

  run_config.gamma = 1.4;

  xt = x1 - 5.0;
  yt = x2 - 5.0;
  r2 = xt*xt + yt*yt;

  eps = 5.0;
  s   = 1.0;
  T   = 1.0 - (run_config.gamma - 1.0)*eps*eps/(8.0*run_config.gamma*CONST_PI*CONST_PI)*exp(1.0 - r2);

  k0      = eps/(2.0*CONST_PI)*exp(0.5*(1.0 - r2));
  v[RHO] = pow(T/s, 1.0/(run_config.gamma - 1.0));
  v[VX1] = 1.0 - k0*yt;
  v[VX2] = 1.0 + k0*xt;
  v[VX3] = 0.0;
  v[PRS] = T*v[RHO];

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
/*!
*/
////////////////////////////////////////////////////////////////////////
{ 
}

////////////////////////////////////////////////////////////////////////
/*! 
*/
////////////////////////////////////////////////////////////////////////
#if MAPPED_GRID == YES
double MappedGrid (Grid *grid, int dir, double xi)
{
  if (dir == IDIR) return MappedGrid_Sinh(grid, dir, xi, 3);
  // if (dir == IDIR) {
  //   double q1  = 0.25, q2  = 0.75;
  //   double xq1 = 4.5, xq2 = 5.5;
  //   return MappedGrid_SUS(grid, dir, xi, q1, xq1, q2, xq2);
  // }
  if (dir == JDIR) return MappedGrid_Sin(grid, dir, xi);
  if (dir == KDIR) return MappedGrid_Uniform(grid, dir, xi);
  return -1.0;
}
#endif
