////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Resistive RMHD rotor test problem in 2 dimensions.

  The 2D rotor problem [dZBL03][Mig12] consists of a rapidly spinning disk
  embedded in a uniform background medium treaded by a costant magnetic field.
  The initial conditin reads as
  \f[
     \left(\rho, p, v_x, v_y, B_x\right) = \left\{\begin{array}{ll}
      \left(10,  1,-\omega y,\,\omega x,\, 1\right) & \quad\mathrm{for}\quad r < 0.1
      \\[0.5em]
      \left(1,\, 1,\, 0,\, 0,\, 1\right) & \quad\mathrm{otherwise}
    \end{array}\right.
  \f]
  where \f$\omega = 8.5\f$ is the angular frequency of rotation.
  The electric field is initialized to its ideal value while remaining
  variables are set to zero.
  The computational domain is the unit square and outflow boundary conditions
  are imposed everywhere.

  A list of tested configurations is given in the following table:

  <CENTER>
  Conf.| DIM       | divB          | eta   
  -----|-----------|---------------|-------
  #01  | 2         | GLM           | 1.e-6
  #03  | 2         | GLM           | 1.e-1
  #04  | 2         | CT/ARITHMETIC | 1.e-6
  #06  | 3         | CT/MAXWELL    | 1.e-6
  #07  | 3         | CT/MAXWELL    | 1.e-6
  #08  | 2 (ho)    | CT/MAXWELL    | 1.e-6
  #10  | 2 (ho)    | CT/MAXWELL    | 1.e-4
  </CENTER>

  The final solution is computed at t = 0.3.
  Conf. #06 and #07 are 2D solutions computed in the XY plane
  and YZ plane. Used to check invariance under rotations.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Aug 07, 2024

  \b Reference: 
     - [BdZ13] "A fully covariant mean-field dynamo closure for numerical 3 + 1
                resistive GRMHD"
        Bucciantini \& Del Zanna,  MNRAS (2013) 428, 71
        
     - [MAR17] "An HLLC Riemann Solver for Resistive RMHD"
       Miranda-Aranguren, Aloy, Rembiasz, MNRAS (2017) ?? - ?? 
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifndef ROTATE 
  #define ROTATE 0
#endif

static void VectorCyclicShift(double *, int);

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  double r, r0, Bx, omega;
  double x[3] = {x1, x2, x3};

  run_config.gamma = (4.0/3.0);

  VectorCyclicShift(x, -ROTATE);  // Get coordinates in unrotated frame 
  x1 = x[0]; x2 = x[1]; x3 = x[2];

  omega   = 8.5;
  Bx      = 1.0;
  r0      = 0.1;

  r = x1*x1 + x2*x2;
  r = sqrt(r);

  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = 1.0;

  v[BX1] = Bx;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  if (r <= r0) {
    v[RHO] = 10.0;
    v[VX1] = -omega*x2;
    v[VX2] =  omega*x1;
  }

  #if ASSIGN_VECTOR_POTENTIAL == YES
  v[AX1] = v[AX2] = 0.0;
  v[AX3] = v[BX1]*x2;
  #endif
  
  #if PHYSICS == ResRMHD
  v[EX1] = -(v[VX2]*v[BX3] - v[VX3]*v[BX2]);
  v[EX2] = -(v[VX3]*v[BX1] - v[VX1]*v[BX3]);
  v[EX3] = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]);
  #ifdef CRG
  v[CRG] = 0.0;
  #endif
  #endif

  VectorCyclicShift(v+AX1, ROTATE);
  VectorCyclicShift(v+BX1, ROTATE);
  VectorCyclicShift(v+EX1, ROTATE);
  VectorCyclicShift(v+VX1, ROTATE);
}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*! 
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
void VectorCyclicShift(double *v, int shift)
/*!  
  Perform a cyclic perumutation of a 3D vector. 
  On output, each vector components has been shifted, e.g.
  
  v = {A,B,C} --> {B,C,A}   if shift = 1, -2
  v = {A,B,C} --> {C,A,B}   if shift = 2, -1
*/
////////////////////////////////////////////////////////////////////////
{
  int i;
  double tmp[3];

  if (shift < 0) shift = 3 - abs(shift)%3;

  for (i = 0; i < 3; i++) tmp[i] = v[i];
  for (i = 0; i < 3; i++) {
    v[(i+shift)%3] = tmp[i];
  }
  return;

  // printf ("----------- Shift = %d  -------------------------------\n", shift);
  // for (i = 0; i < 3; i++) {
  //   printf ("Original[%d] = %f, Rotated[%d] = %f\n",
  //            i,tmp[i], i, v[i]);
  // }
}

