////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Double Mach reflection test problem.

  Sets the initial condition for a planar shock front making an angle of
  \f$ \pi/3 \f$ with a reflecting wall:
  \f[
     \left(\rho,v_x,v_y,p\right)
     =
    \left\{\begin{array}{ll}
    (1.4,0,0,1) & \qquad {\rm for} \quad x >x_s(0)  \\[0.5em]
    (8,8,25,-8.25,116.5) & \qquad {\rm otherwise}
    \end{array}\right.
  \f]
  where \f$x_s(t) = 10t/\sin\alpha + 1/6 + y/\tan\alpha\f$ is the shock
  position.
  The ideal equation of state with \f$ \Gamma = 1.4\f$ is used.
  The wedge is represented by a reflecting boundary starting at x=1/6
  along the lower y-boundary.
  As the shock reflects off the lower wall, a complex flow structure
  develops with two curved reflected shocks propagating at directions
  almost orthogonal to each other and a tangential discontinuity
  separating them.
  At the wall, a pressure gradient sets up a denser fluid jet
  propagating along the wall. Kelvin-Helmholtz instability
  patterns may be identified with the rolls developing at
  the slip line.
  This feature is very sensitive to numerical diffusion and it
  becomes visible at high resolution and/or low dissipative schemes.

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CARTESIAN  |2  |PARABOLIC  |RK3        |roe
   #02 |HD  |CARTESIAN  |2  |PARABOLIC  |RK3        |roe
   #03 |HD  |CARTESIAN  |2  |WENOZ      |RK3        |roe
   #04 |HD  |CARTESIAN  |2  |WENOZ      |RK3        |two_shock
   #05 |HD  |CARTESIAN  |2  |MP5        |RK3        |hllc
  <\CENTER>

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 15, 2026

  \b References
    - "The Numerical Simulation of Two-Dimensional Fluid Flow with Strong Shocks"\n
       Woodward, P.R., & Colella, P., JCP (1984) 54, 115.
    - "The PLUTO Code for computational astrophysics"\n
       Mignone et al., ApJS(2007) 170,228
    - "Comparison of some FLux Corrected Transport and Total variatiion
       diminishing numerical schemes for hydrodynamics and magnetohydrodynamic
       problems" \n
       Toth, G., Odstrcil, D., JCP (1996) 128, 82.

*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
*/
////////////////////////////////////////////////////////////////////////
{
  double alpha, xs;

  run_config.gamma = 1.4;

  alpha   = 1.0/3.0*CONST_PI;
  xs      = 1.0/6.0 + x2/tan(alpha);

  if (x1 > xs){
    v[RHO] = 1.4;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[PRS] = 1.0;
  }else{
    v[RHO] = 8.0;
    v[VX1] =  8.25*sin(alpha);
    v[VX2] = -8.25*cos(alpha);
    v[PRS] = 116.5;
  }
  v[VX3] = 0.0;
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
  Assign user-defined boundary conditions:
  - left side (x beg):   constant shocked values
  - bottom side (y beg): constant shocked values for
  <tt> x < 1/6 </tt> and reflective boundary
  otherwise.
  - top side (y end):   time-dependent boundary:
  for \f$ x < x_s(t) = 10 t/\sin\alpha + 1/6 + 1.0/\tan\alpha \f$ we use
  fixed (post-shock) values. Unperturbed values otherwise.
*/
////////////////////////////////////////////////////////////////////////
{
  int     i, j, k;
  double  x1;
  double  alpha,xs;
  double  time = data->Dts->time;

  alpha = 60./180.*CONST_PI;

  if (side == X1_BEG){

    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k) {
      data->Vc(i,j,k,RHO) =   8.0;
      data->Vc(i,j,k,VX1) =   8.25*sin(alpha);
      data->Vc(i,j,k,VX2) = - 8.25*cos(alpha);
      data->Vc(i,j,k,VX3) =   0.0;
      data->Vc(i,j,k,PRS) = 116.5;
    }

  } else if (side == X2_BEG){

    int jbeg = grid->lbeg[JDIR];
    #pragma acc parallel loop collapse(3) present(data, grid, box)
    BOX_LOOP(*box,i,j,k) {
      x1 = grid->xC[IDIR][i];
      if (x1 < 1.0/6.0){
        data->Vc(i,j,k,RHO) =   8.0;
        data->Vc(i,j,k,VX1) =   8.25*sin(alpha);
        data->Vc(i,j,k,VX2) = - 8.25*cos(alpha);
        data->Vc(i,j,k,VX3) =   0.0;
        data->Vc(i,j,k,PRS) = 116.5;
      }else{                                   // reflective boundary 
        data->Vc(i,j,k,RHO) =   data->Vc(i,(2*jbeg-j-1),k,RHO);
        data->Vc(i,j,k,VX1) =   data->Vc(i,(2*jbeg-j-1),k,VX1);
        data->Vc(i,j,k,VX2) = - data->Vc(i,(2*jbeg-j-1),k,VX2);
        data->Vc(i,j,k,VX3) =   0.0;
        data->Vc(i,j,k,PRS) =  data->Vc(i,(2*jbeg-j-1),k,PRS);
      }
    }

  } else if (side == X2_END) {

    #pragma acc parallel loop collapse(3) present(data, grid, box)
    BOX_LOOP(*box,i,j,k){
      x1 = grid->xC[IDIR][i];
      xs = 10.0*time/sin(alpha) + 1.0/6.0 + 1.0/tan(alpha);
      if (x1 < xs){
        data->Vc(i,j,k,RHO) =   8.0;
        data->Vc(i,j,k,VX1) =   8.25*sin(alpha);
        data->Vc(i,j,k,VX2) = - 8.25*cos(alpha);
        data->Vc(i,j,k,VX3) =   0.0;
        data->Vc(i,j,k,PRS) = 116.5;
      }else{                                   // reflective boundary 
        data->Vc(i,j,k,RHO) = 1.4;
        data->Vc(i,j,k,VX1) = 0.0;
        data->Vc(i,j,k,VX2) = 0.0;
        data->Vc(i,j,k,VX3) = 0.0;
        data->Vc(i,j,k,PRS) = 1.0;
      }
    }
  }
}
