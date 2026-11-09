////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Two-Dimensional Riemann problem.

  Sets the initial condition for a 2D Riemann problem in terms of 4 four
  constant states at each corner of the computational domain.
  For example, the top right corner, \f$x > 0\f$, \f$y > 0\f$, is denoted
  with \c PP and its initial values are:

  \f[
    \rho = \rho_{\rm PP}
    \,,\quad\quad
    P = P_{\rm PP}
    \,,\quad\quad
    v_x = v_{x\,{\rm PP}}
    \,,\quad\quad
    v_y = v_{y\,{\rm PP}}\,.
  \f]

  Similar conditions hold for the bottom right corner 
  (\f$x > 0\f$, \f$y < 0\f$; \c PM), the top left 
  (\f$x < 0\f$, \f$y > 0\f$; \c MP), and the bottom left 
  (\f$x < 0\f$, \f$y < 0\f$; \c MM).
  The 4 flow quantities of each corner are defined by the 16 parameters
  that are read from pluto.ini.
  - <tt>DN_PP, PR_PP, VX_PP, VY_PP</tt>
  - <tt>DN_MP, PR_MP, VX_MP, VY_MP</tt>
  - <tt>DN_PM, PR_PM, VX_PM, VY_PM</tt>
  - <tt>DN_MM, PR_MM, VX_MM, VY_MM</tt>

  - Configurations #01, #02 and #05 involve the interaction between four contact
    discontinuities (see fig. below).
  - Configurations #03 and #04 (with AMR) involve the interaction of shocks.

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CARTESIAN  |2  |LINEAR     |RK3        |hll
   #02 |HD  |CARTESIAN  |2  |PARABOLIC  |RK3        |roe
   #03 |HD  |CARTESIAN  |2  |LINEAR     |RK2        |hllc
  <\CENTER>

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 15, 2026

  \b References
     - C. W. Schulz-Rinne, "Classification of the Riemann problem for
        two-dimensional gas dynamics", SIAM J. Math. Anal., 24
        (1993), pp. 76-88.
     - C. W. Schulz-Rinne, J. P. Collins and H. M. Glaz,
       "Numerical solution of the Riemann problem for two-dimensional gas
        dynamics", SIAM J. Sci. Comp., 14 (1993), pp. 1394-1414.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  double *inputParam = run_config.inputParam;
  run_config.gamma = 1.4;

  if (x > 0.0 && y > 0.0){
    v[RHO] = inputParam[DN_PP];
    #if EOS == IDEAL
    v[PRS] = inputParam[PR_PP];
    #endif
    v[VX1] = inputParam[VX_PP];
    v[VX2] = inputParam[VY_PP];
  }else if(x < 0.0 && y > 0.0){
    v[RHO] = inputParam[DN_MP];
    #if EOS == IDEAL
    v[PRS] = inputParam[PR_MP];
    #endif
    v[VX1] = inputParam[VX_MP];
    v[VX2] = inputParam[VY_MP];
  }else if(x < 0.0 && y < 0.0){
    v[RHO] = inputParam[DN_MM];
    #if EOS == IDEAL
    v[PRS] = inputParam[PR_MM];
    #endif
    v[VX1] = inputParam[VX_MM];
    v[VX2] = inputParam[VY_MM];
  }else if(x > 0.0 && y < 0.0){
    v[RHO] = inputParam[DN_PM];
    #if EOS == IDEAL
    v[PRS] = inputParam[PR_PM];
    #endif
    v[VX1] = inputParam[VX_PM];
    v[VX2] = inputParam[VY_PM];
  }
  v[VX3] = 0.0;

// -- Initialize Tracers -- 

#if NTRACER == 2
  v[TR1] = (x*y) > 0.0;
  v[TR2] = (x*y) < 0.0;
#endif
 
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
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
/*!
  Assign user-defined boundary conditions.
  
  \param [in,out] data  pointer to the PLUTO data structure
  \param [in]     box   pointer to a RBox structure containing the 
                        boundary box corresponding to side
  \param [in]     side  side of computational domain (special value side == 0 
                        is used to control INTERNAL_BOUNDARY)
  \param [in]     grid  pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{ 
}
