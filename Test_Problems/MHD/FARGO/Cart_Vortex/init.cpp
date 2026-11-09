#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
 * Assign initial conditions in terms of primitive variables
*/
////////////////////////////////////////////////////////////////////////
{
  double rnd;
  double x, y, z = 0.0, r2; 
  double kp, mu;
  double dvx, dvy, dBx, dBy, arg;
  double *inputParam = run_config.inputParam;

#if EOS == IDEAL  
  run_config.gamma = 5./3.;
#endif

int64_t ix = (int64_t)(x1 * 1e6);
int64_t iy = (int64_t)(x2 * 1e6);
int64_t iz = (int64_t)(x3 * 1e6);
rnd = GridRandomNumber_new(ix, iy, iz, g_rank);

  DIM_EXPAND(x = x1;, y = x2;, z = x3;)

  kp  = inputParam[KPAR];
  mu  = inputParam[MUPAR];
  r2  = x*x + y*y + z*z;
  arg = (1.0 - r2);

  dvx = -y*kp*exp(0.5*arg);
  dvy =  x*kp*exp(0.5*arg);

  dBx = -y*mu*exp(0.5*arg);
  dBy =  x*mu*exp(0.5*arg);

  v[RHO]= 1.0;
  v[VX1]= dvx + inputParam[PERT_AMPL]*(rnd-0.5)*exp(-x1*x1/4.0);

  #ifdef FARGO
  v[VX2]     = dvy;
  v[FARGO_W] = 0.5*inputParam[MACH]*tanh(x1);
  #else
  v[VX2] = dvy + 0.5*inputParam[MACH]*tanh(x1);
  #endif
  v[VX3] = 0.0;

  #if DIMENSIONS == 3
  v[VX3] = 0.1*sin(2.0*CONST_PI*x3/10.0);
  #endif

  #if EOS == IDEAL
  v[PRS] = 1.0/run_config.gamma + exp(arg)*(mu*mu*(1.0 - r2 + z*z) - v[RHO]*kp*kp)*0.5;
  #endif

#if PHYSICS == MHD
  v[BX1] = dBx;
  v[BX2] = dBy;  
  v[BX3] = 0.0;

  #ifdef STAGGERED_MHD
  v[AX1] = v[AX2] = 0.0;
  v[AX3] = mu*exp(0.5*arg);
  #endif
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
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
////////////////////////////////////////////////////////////////////////
{
}
////////////////////////////////////////////////////////////////////////
double BodyForcePotential(const State& state, const RunConfig& run_config)
////////////////////////////////////////////////////////////////////////
{
  return 0.0;
}
#endif
