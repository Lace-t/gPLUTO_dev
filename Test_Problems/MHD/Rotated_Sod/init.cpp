////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD rotated shock tube problem.

  2D-3D rotated version of standard shock-tube problem.
  The initial condition is expressed in the original 
  (non-rotated) 1D frame perpendicular to the plane of 
  discontinuity.

  For each case, we rotate the initial condition by 
  an angle \f$\alpha\f$ around the \f$z\f$-axis
  and then by an angle \f$\gamma\f$ around the \f$y\f$-axis.
 
  The coordinate transformation used for 2D/3D rotation is explained 
  in [MT10] (S 4.2.2) and leads to a planar symmetry condition whereby,
  for any quantity \f$q(\vec{x})\f$, one must have 
  \f$q(\vec{x}+\vec{s}) = q(\vec{x})\f$, 
  where \f$\vec{s} = (n_x\Delta x,\, n_y\Delta y,\, n_z\Delta z)\f$ 
  corresponds to an integer shift of cells.
 
  The runtime parameters that are read from \c pluto.ini are:
  - <tt>run_config.inputParam[RHO_LEFT]</tt>: sets the left density, value: 1.08;
  - <tt>run_config.inputParam[VX_LEFT]</tt>: sets the left velocity in the x direction, value: 1.2;
  - <tt>run_config.inputParam[VY_LEFT]</tt>: sets the left velocity in the y direction, value: 0.01;
  - <tt>run_config.inputParam[VZ_LEFT]</tt>: sets the left velocity in the z direction, value: 0.5;
  - <tt>run_config.inputParam[BY_LEFT]</tt>: sets the left magnetic field in the y direction, value: 3.6;
  - <tt>run_config.inputParam[BZ_LEFT]</tt>: sets the left magnetic field in the z direction, value: 2.0;
  - <tt>run_config.inputParam[PR_LEFT]</tt>: sets the left pressure, value: 0.95;
  - <tt>run_config.inputParam[RHO_RIGHT]</tt>: sets the right density, value: 1.0;
  - <tt>run_config.inputParam[VX_RIGHT]</tt>: sets the right velocity in the x direction, value: 0.0;
  - <tt>run_config.inputParam[VY_RIGHT]</tt>: sets the right velocity in the y direction, value: 0.0;
  - <tt>run_config.inputParam[VZ_RIGHT]</tt>: sets the right velocity in the z direction, value: 0.0;
  - <tt>run_config.inputParam[BY_RIGHT]</tt>: sets the right magnetic field in the y direction, value: 4.0;
  - <tt>run_config.inputParam[BZ_RIGHT]</tt>: sets the right magnetic field in the z direction, value: 2.0;
  - <tt>run_config.inputParam[PR_RIGHT]</tt>: sets the right pressure, value: 1.0;
  - <tt>run_config.inputParam[BX_CONST]</tt>: sets the constant magnetic field in the x direction, value: 2.0;
  - <tt>run_config.inputParam[INT_JX]</tt>: sets the integer shift in the x direction, value: 1;
  - <tt>run_config.inputParam[INT_JY]</tt>: sets the integer shift in the y direction, value: -2;
  - <tt>run_config.inputParam[INT_KX]</tt>: sets the integer shift in the k direction, value: 1;
  - <tt>run_config.inputParam[INT_KZ]</tt>: sets the integer shift in the z direction, value: 4;

  A list of the available configurations is given in the following table:

  <CENTER>
  Conf.  | GEOMETRY  |DIM| T.STEPPING |RECONSTRUCTION| divB     |  RS  |  HO  
  -------|-----------|---|------------|--------------|----------|------|------
   #01   |CARTESIAN  | 2 |  RK2       |  LINEAR      | GLM      | roe  |  NO  
   #02(*)|CARTESIAN  | 2 |  RK2       |  LINEAR      | GLMex    | roe  |  NO 
   #03   |CARTESIAN  | 2 |  RK2       |  LINEAR      |CT Arithm | tvdlf|  NO  
   #04   |CARTESIAN  | 2 |  RK3       |  LINEAR      |CT Contact| hllc |  NO  
   #05   |CARTESIAN  | 3 |  RK3       |  MP5         |CT gforce |gforce|  NO  
   #06   |CARTESIAN  | 3 |  RK2       |  LINEAR      |CT hlld   | hlld |  NO  
   #07   |CARTESIAN  | 3 |  RK3       |  WENO3       | GLM      | roe  |  NO  
   #08   |CARTESIAN  | 3 |  RK3       |  WENOZ       |CT hll    | hllc |  NO
   #09   |CARTESIAN  | 3 |  HANCOK    |  LINEAR      | GLM      |      |  YES 
   #10   |CARTESIAN  | 2 |  RK2       |  LINEAR      | CT       |      |  YES  
  </CENTER>
  (*) To be implemented 

  \author A. Mignone (andrea.mignone@unito.it) \n
          V. Berta   (vittoria.berta@unito.it)

  \date   Jul 23, 2024

  \b Reference: 
     - [MT10] "A second-order unsplit Godunov scheme 
        for cell-centered MHD: The CTU-GLM scheme"
        Mignone \& Tzeferacos, JCP (2010) 229, 2217
     - [BM24] "A 4th-order accurate finite volume method
        for ideal classical and special relativistic MHD 
        based on pointwise reconstructions"
        Berta, Mignone, Bugli & Mattia JCP (2024) 499, 112701
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "rotate.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #ifndef GAMMA_EOS
  run_config.gamma=1.6666666666666667;
  #endif
  
//   static int first_call = 1;
  int    jx = run_config.inputParam[INT_JX];
  int    jy = run_config.inputParam[INT_JY];
  int    kx = run_config.inputParam[INT_KX];
  int    kz = run_config.inputParam[INT_KZ];
  double xm = 0.5*(run_config.domEnd[IDIR] + run_config.domBeg[IDIR]); // Initial discontinuity 
                                                                       // location             
//   double ym = 0.5*(run_config.domEnd[JDIR] + run_config.domBeg[JDIR]); /* Initial discontinuity */
//                                                                        /* location             */
  double xcoord[3];

  // 0. Initialize Rotation

  RotateSet(jx, jy, kx, kz, run_config);

  // 1. Rotate coordinate to 1D frame.
  //    The center of rotation is the point (0.5,0) of the original frame.

  xcoord[0] = x-xm; xcoord[1] = y; xcoord[2] = z;
  RotateVector(xcoord, 1);
  x = xcoord[0]; y = xcoord[1]; z = xcoord[2];
  
  // 2. Assign initial conditions.
  //    We use a small negative value (-1.e-12) so that the row of zones 
  //    immediately above the midplane line is discontinuous at  x = 0;

  if (x <= 1.e-6){  
     v[RHO] = run_config.inputParam[RHO_LEFT];
     v[VX1] = run_config.inputParam[VX_LEFT];
     v[VX2] = run_config.inputParam[VY_LEFT];
     v[VX3] = run_config.inputParam[VZ_LEFT];
     v[BX1] = run_config.inputParam[BX_CONST];
     v[BX2] = run_config.inputParam[BY_LEFT];
     v[BX3] = run_config.inputParam[BZ_LEFT];
     v[PRS] = run_config.inputParam[PR_LEFT];
  }else{
     v[RHO] = run_config.inputParam[RHO_RIGHT];
     v[VX1] = run_config.inputParam[VX_RIGHT];
     v[VX2] = run_config.inputParam[VY_RIGHT];
     v[VX3] = run_config.inputParam[VZ_RIGHT];
     v[BX1] = run_config.inputParam[BX_CONST];
     v[BX2] = run_config.inputParam[BY_RIGHT];
     v[BX3] = run_config.inputParam[BZ_RIGHT];
     v[PRS] = run_config.inputParam[PR_RIGHT];
  }

  #if DIVIDE_BY_4PI == YES
  v[BX1] /= sqrt(4.0*CONST_PI);
  v[BX2] /= sqrt(4.0*CONST_PI);
  v[BX3] /= sqrt(4.0*CONST_PI);
  #endif

  #if PHYSICS == ResRMHD
  v[EX1] = -(v[VX2]*v[BX3] - v[VX3]*v[BX2]); 
  v[EX2] = -(v[VX3]*v[BX1] - v[VX1]*v[BX3]); 
  v[EX3] = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]); 

  v[CRG] = 0.0;
  #endif

  v[AX1] = 0.0;
  v[AX2] = x*v[BX3];
  v[AX3] = y*v[BX1] - x*v[BX2];

  #ifdef PSI_GLM 
  v[PSI_GLM] = 0.0;
  #endif

  // 3. Rotate vectors from 1D frame to 2D  or 3D frame
  //    [counter-clockwise]

  RotateVector (v + AX1, -1);
  RotateVector (v + VX1, -1);
  RotateVector (v + BX1, -1);

  #if PHYSICS == ResRMHD
  RotateVector (v + EX1, -1);
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
