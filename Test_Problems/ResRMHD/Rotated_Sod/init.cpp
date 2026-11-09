////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Rotated shock tube problem.


  Conf #01: Intended to replace MHD/Rotated_Sod/#01

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Oct 16, 2018

  \b Reference:
     - [MT10] "A second-order unsplit Godunov scheme for cell-centered MHD: The CTU-GLM scheme"
        Mignone \& Tzeferacos, JCP (2010) 229, 2217
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "rotate.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*
*/ 
////////////////////////////////////////////////////////////////////////
{
  double *inputParam = run_config.inputParam;
  int    jx = inputParam[INT_JX];
  int    jy = inputParam[INT_JY];
  int    kx = inputParam[INT_KX];
  int    kz = inputParam[INT_KZ];
  double xcoord[3];
  
  run_config.gamma = inputParam[GAMMA_EOS];

  // 0. Initialize Rotation

  RotateSet(jx, jy, kx, kz, run_config);

  // 1. Rotate coordinate to 1D frame.
  //    The center of rotation is the point (0.5,0) of the original frame.

  xcoord[0] = x; xcoord[1]= y; xcoord[2] = z;
  RotateVector(xcoord, 1);
  x = xcoord[0]; y = xcoord[1]; z = xcoord[2];

  // 2. Assign initial conditions.
  //    We use a small negative value (-1.e-12) so that the row of zones 
  //    immediately above the midplane line is discontinuous at  x = 0;

  if ( x < -1.e-12){
    v[RHO] = inputParam[RHO_LEFT];
    v[VX1] = inputParam[VX1_LEFT];
    v[VX2] = inputParam[VX2_LEFT];
    v[VX3] = inputParam[VX3_LEFT];
    v[BX1] = inputParam[BX1_CONST];
    v[BX2] = inputParam[BX2_LEFT];
    v[BX3] = inputParam[BX3_LEFT];
    v[PRS] = inputParam[PRS_LEFT];
  }else{
    v[RHO] = inputParam[RHO_RIGHT];
    v[VX1] = inputParam[VX1_RIGHT];
    v[VX2] = inputParam[VX2_RIGHT];
    v[VX3] = inputParam[VX3_RIGHT];
    v[BX1] = inputParam[BX1_CONST];
    v[BX2] = inputParam[BX2_RIGHT];
    v[BX3] = inputParam[BX3_RIGHT];
    v[PRS] = inputParam[PRS_RIGHT];
  }

  #if DIVIDE_BY_4PI == TRUE
  v[BX1] /= sqrt(4.0*CONST_PI);
  v[BX2] /= sqrt(4.0*CONST_PI);
  v[BX3] /= sqrt(4.0*CONST_PI);
  #endif

  #if PHYSICS == ResRMHD
  v[EX1] = -(v[VX2]*v[BX3] - v[VX3]*v[BX2]);
  v[EX2] = -(v[VX3]*v[BX1] - v[VX1]*v[BX3]);
  v[EX3] = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]);

  #ifdef CRG
  v[CRG] = 0.0;
  #endif
  #endif

  v[AX1] = 0.0;
  v[AX2] = x*v[BX3];
  v[AX3] = y*v[BX1] - x*v[BX2];

  // 3. Rotate vectors from 1D frame to 2D frame
  //    [counter-clockwise]

  RotateVector (v + AX1, -1);
  RotateVector (v + VX1, -1);
  RotateVector (v + BX1, -1);

  #if PHYSICS == ResRMHD
  RotateVector (v + EX1, -1);
  #endif

  #if ENTROPY_SWITCH != NO
  v[ENTR] = Entropy(v, run_config);
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
/*!
  
*/
////////////////////////////////////////////////////////////////////////
{
}
