////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Two-dimensional Riemann problem.

  Sets the initial condition for the 2D Riemann problem
  described in Mignone et al. (2005).
  The computational domain is initially divided into four
  states and the outcoming wave pattern involves the formtation of
  2 shocks and two contact waves.

  \authors A. Mignone (mignone@ph.unito.it)
  \date   July 09, 2014

  \b Reference:
     - "The Piecewise Parabolic Method for Multidimensional Relativistic
        Fluid Dynamics", Mignone, Plewa & Bodo, ApJS (2005)
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*
*/ 
////////////////////////////////////////////////////////////////////////
{
  v[VX1] = v[VX2] = 0.0;
  if (x > 0.0 && y > 0.0){
    v[RHO] = 5.477875e-3;
    v[PRS] = 2.762987e-3;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
  }else if(x < 0.0 && y > 0.0){
    v[RHO] = 0.1;
    v[PRS] = 1.0;
    v[VX1] = 0.99;
    v[VX2] = 0.0;
  }else if(x < 0.0 && y < 0.0){
    v[RHO] = 0.5;
    v[PRS] = 1.0;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
  }else if(x > 0.0 && y < 0.0){
    v[RHO] = 0.1;
    v[PRS] = 1.0;
    v[VX1] = 0.0;
    v[VX2] = 0.99;
  }

  #if ENTROPY_SWITCH != NO
  v[ENTR] = Entropy(v, run_config);
  #else
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
