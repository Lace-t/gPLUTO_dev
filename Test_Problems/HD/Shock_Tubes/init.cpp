////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Hydro shock tube problems

  Shock tube problems are popular benchmark used to test
  shock-capturing schemes.
  It is a one-dimensional problem with initial condition given
  by a discontinuity separating two constant states:
  \f[
     {\bf V} = \left\{ \begin{array}{ll}
     \left(\rho,\, v_x,\, v_y,\, v_z,\, p\right)_L & 
    \qquad\mathrm{for}\quad x < 0.5
     \\[0.5em]
   \left(\rho,\, v_x,\, v_y,\, v_z,\, p\right)_R & 
      \qquad\mathrm{for}\quad x > 0.5
   \end{array}\right.
  \f]

  Left and right states are provided at runtime using the parameters 
  defined in pluto.ini. 

  Below the pre-defined configurations:

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CARTESIAN  |1  |LINEAR     |RK2        |roe
   #02 |HD  |CARTESIAN  |1  |MP5        |RK3        |two_shock
   #03 |HD  |CARTESIAN  |1  |PARABOLIC  |RK3        |roe
   #04 |HD  |CARTESIAN  |1  |LimO3      |RK3        |hllc
   #05 |HD  |CARTESIAN  |1  |WENOZ      |RK3        |hll
   #06 |HD  |CARTESIAN  |1  |PARABOLIC  |RK3        |hllc
   #07 |HD  |CARTESIAN  |1  |LINEAR     |RK2        |roe
   #08 |HD  |CARTESIAN  |1  |LINEAR     |RK2        |roe
   #09 |HD  |CARTESIAN  |1  |WENO3      |RK3        |tvdlf
  <\CENTER>

  \author A. Mignone (andrea.mignone@unito.it)

  \date   Apr 15, 2025

  \b References
     - Sod, G. A. (1978).
       "A Survey of Several Finite Difference Methods for Systems of
        Nonlinear Hyperbolic Conservation Laws".
        JCP (1978) 27:1-31
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*
*/
////////////////////////////////////////////////////////////////////////
{
  double *inputParam = run_config.inputParam;

  run_config.gamma = 1.4;

  if (x < 0.5){
    v[RHO] = inputParam[RHO_LEFT];
    v[VX1] = inputParam[VX1_LEFT];
    v[VX2] = inputParam[VX2_LEFT];
    v[VX3] = inputParam[VX3_LEFT];
    #if EOS == IDEAL
    v[PRS] = inputParam[PRS_LEFT];
    #endif
  }else{
    v[RHO] = inputParam[RHO_RIGHT];
    v[VX1] = inputParam[VX1_RIGHT];
    v[VX2] = inputParam[VX2_RIGHT];
    v[VX3] = inputParam[VX3_RIGHT];
    #if EOS == IDEAL
    v[PRS] = inputParam[PRS_RIGHT];
    #endif
  }

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*
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