#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *us, double x, double y, double z, RunConfig &run_config)
/*
 * Assign initial conditions in terms of primitive variables
*/
////////////////////////////////////////////////////////////////////////
{
  run_config.unit_velocity = (1.e5);
  us[RHO]=  25./9.0;

#if DIMENSIONS == 2
  us[VX1]= -sin(2.0*M_PI*y);
  us[VX2]=  sin(2.0*M_PI*x);
  us[VX3]=  0.0;
#elif DIMENSIONS == 3 // Taken from Helzel, Rossmanith 2011 
  double eps = 0.2;
  double fact = (1.0 + eps*sin(2.0*M_PI*z));
  us[VX1]= -sin(2.0*M_PI*y)*fact;
  us[VX2]=  sin(2.0*M_PI*x)*fact;
  us[VX3]=  eps*sin(2.0*M_PI*z);
#endif


  #if PHYSICS == MHD
  us[BX1]= -sin(2.0*M_PI*y);
  us[BX2]=  sin(4.0*M_PI*x);
  us[BX3]=  0.0;
  #endif

  us[PRS]=  5.0/3.0;

  #ifdef GLM_MHD
  us[PSI_GLM] = 0.0;
  #endif
  #if NTRACER == 1
  us[TR1]=  y > 0.5 ? 1.0:0.0;
  #endif

  us[AX1] = 0.0;
  us[AX2] = 0.0;
  us[AX3] = (cos(2.0*M_PI*y) + 0.5*cos(4.0*M_PI*x))/(2.0*M_PI);

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
  Perform runtime data analysis.
  
  \param [in] d the PLUTO Data structure
  \param [in] grid   pointer to array of Grid structures
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

#if BODY_FORCE != NO
////////////////////////////////////////////////////////////////////////
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
////////////////////////////////////////////////////////////////////////
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

////////////////////////////////////////////////////////////////////////
double BodyForcePotential(double x1, double x2, double x3)
////////////////////////////////////////////////////////////////////////
{
  return 0.0;
}
#endif

