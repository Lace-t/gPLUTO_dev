////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD blast wave.

  In this example the initial conditions consists of a static medium with 
  uniform density  \f$ \rho = 1 \f$ while pressure and magnetic field 
  are given by 
  \f[
     p = \left\{\begin{array}{ll}
           p_{\rm in}   & \qquad\mathrm{for}\quad  r < r_0 \\[0.5em]
           1            & \qquad\mathrm{otherwise} \\ 
     \end{array}\right.\,;\qquad
     \vec{B} = B_0\hat{\bf x}
  \f]
  The values \f$p_{\rm in},\, p_{\rm out},\, B_0,\, \theta,\, \phi,\, r_0\f$ 
  are control parameters that can be changed from \c pluto.ini using, 
  respectively,

  -# <tt>inputParam[P_IN]</tt>
  -# <tt>inputParam[BMAG]</tt>

  The over-pressurized region drives a blast wave delimited by an outer
  fast forward shock propagating (nearly) radially while magnetic field 
  lines pile up behind the shock thus building a region of higher magnetic 
  pressure.
  In these regions the shock becomes magnetically dominated and only weakly
  compressive (\f$\delta\rho/\rho\sim 1.2\f$ in both cases).
  The inner structure is delimited by an oval-shaped slow shock adjacent to a
  contact discontinuity and the  two fronts tend to blend together as the
  propagation becomes perpendicular to the field lines.
  The magnetic energy increases behind the fast shock and decreases 
  downstream of the slow shock.
  The resulting explosion becomes highly anisotropic and magnetically confined.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Sep 10, 2024
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
{
  double *inputParam = run_config.inputParam;
  double r, theta, phi, B0;
  
  r = DIM_EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r);

  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = 1.0/run_config.gamma;
  B0     = sqrt(2.0*v[PRS]/inputParam[BETA]);

  if (r <= 0.1) v[PRS] = inputParam[P_IN]; // Change pressure internally   
 
  v[BX1] = B0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;
  
  v[AX1] = 0.0;
  v[AX2] =  v[BX3]*x1;
  v[AX3] = -v[BX2]*x1 + v[BX1]*x2;
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
/*!
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
