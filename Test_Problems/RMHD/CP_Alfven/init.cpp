////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief RMHD Circularly polarized Alfven waves

  Setup the initial conditions for the large-amplitude circularly
  polarized Alfven wave test, as in section 4.1 of Del Zanna (2007).

  Note: in 2D the solution is rotated around the $z$-axis by
  assuming one  wavelength in both directions.
  This implies that \f$ k_x = 2\pi/L_x\f$ and \f$ k_y = 2\pi/L_y\f$.
  The transformation from the 1D (unrotated) system with primed
  coordinates to the actual computational system is given by
  \f[
      \mathbf{V}  = R\mathbf{V}' \qquad{\rm where}\quad
      R = \left(\begin{array}{ccc}
              \cos\alpha  & -\sin\alpha   &  0 \\[0.5em]
              \sin\alpha  &  \cos\alpha   &  0 \\[0.5em]
                     0  &           0   &  1 \\[0.5em]
          \end{array}\right)
  \f]
  while the invers transormation is
  \f[
      \mathbf{V}'  = R^{-1}\mathbf{V} \qquad{\rm where}\quad
      R^{-1} = \left(\begin{array}{ccc}
                \cos\alpha  &  \sin\alpha   &  0 \\[0.5em]
               -\sin\alpha  &  \cos\alpha   &  0 \\[0.5em]
                       0  &           0   &  1 \\[0.5em]
          \end{array}\right)
  \f]
  Note that the wave phase
  \f[
      \phi = \mathbf{k}\cdot\mathbf{x} = \mathbf{k}'\cdot\mathbf{x}'
  \f]
  is invariant under rotations.

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM| RECONST |TIME_STEP|DIVB| Solver
  -----|----|-----------|---|---------|---------|----|-------
   #01 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hlld
   #02 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     |GLM | tvdlf
   #03 |MHD |CARTESIAN  | 2 | MP5     | RK4 [HO]| CT | gforce
   #04 |MHD |CARTESIAN  | 2 | MP5     | RK4 [HO]| CT | gforce
   #05 |MHD |CARTESIAN  | 2 | WENOZ   | RK4 [HO]| CT | gforce
  </CENTER>

  \author A. Mignone (andrea.mignone@unito.it) \n
          M. Bugli   (matteo.bugli@unito.it) \n
          V. Berta   (vittoria.berta@unito.it)
  \date   Apr 16, 2026

  \b References:
     - Del Zanna et al, A&A (2007) 473, 11-30
     - Berta, Mignone, Bugli & Mattia, JCP (2024) 499, 112701
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "rotate.hpp"

// static int first_call = 1;
// #pragma acc declare create(first_call)

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*!
  The Init() function can be used to assign initial conditions as
  as a function of spatial position.
  
  \param [out] v          pointer to a vector of primitive variables
  \param [in] x1          coordinate point in the 1st dimension
  \param [in] x2          coordinate point in the 2nd dimension
  \param [in] x3          coordinate point in the 3rd dimension
  \param [in] run_config  reference to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
{
  double Lx = run_config.domEnd[IDIR] - run_config.domBeg[IDIR];
  double Ly = run_config.domEnd[JDIR] - run_config.domBeg[JDIR];
  double Lz = run_config.domEnd[KDIR] - run_config.domBeg[KDIR];
  double kx = 2.0*CONST_PI/Lx;  // Choose wavevector so that we have one  
  double ky = 2.0*CONST_PI/Ly;  // wavelenght in each direction.          
  double kz = 2.0*CONST_PI/Lz;   

  // ------------------------------------------------------
  // 0. Define rotation shift vectors from the 
  //    domain size.
  //    Assume jshift_y = kshift_z = 1 and square 
  //    cells (dx=dy=dz).
  // ------------------------------------------------------

  double jshift_y = 1;
  double jshift_x = round(-jshift_y/Ly);
  double kshift_z = 1;
  double kshift_x = round(-kshift_z/Lz);

  #if DIMENSIONS == 1
  jshift_y = 0; jshift_x = 0.0;
  kshift_z = 0; kshift_x = 0.0;
  ky = kz = 0.0;
  #endif

  // ------------------------------------------------------
  // 1. Define solution in the unrotated (1D) frame
  // ------------------------------------------------------

  run_config.gamma = 4.0/3.0;
  v[RHO] = 1.0;
  v[PRS] = 1.0;

  double B0 = 1.0;
  double eta = 1.0;

  // Alfven speed 
  double w  = v[RHO] + run_config.gamma*v[PRS]/(run_config.gamma - 1.0);
  double vA  = B0*B0/(w + B0*B0*(1.0 + eta*eta));
  vA /= 0.5*(1.0 + sqrt(1.0 - 4.0*eta*eta*vA*vA));
  vA  = sqrt(vA);

  // Wave phase (invariant under rotations) 
  double phi = DIM_EXPAND(kx*x, + ky*y, + kz*z);

  v[BX1] = B0;
  v[BX2] = eta*B0*cos(phi);
  v[BX3] = eta*B0*sin(phi);

  v[VX1] = 0.0;
  v[VX2] = -vA*v[BX2]/B0;
  v[VX3] = -vA*v[BX3]/B0;

  // ------------------------------------------------------
  // 2. Rotate vectors
  // ------------------------------------------------------

  #if DIMENSIONS > 1
  RotateSet(jshift_x, jshift_y, kshift_x, kshift_z, run_config);
  RotateVector(v + VX1,-1);
  RotateVector(v + BX1,-1);
  #endif

  // ------------------------------------------------------
  // 3. Compute vector potential
  // ------------------------------------------------------

  double kvec[3]={kx, ky, kz};
  double xvec[3]={x, y, z};
  
  RotateVector(kvec,1);
  RotateVector(xvec,1);

  v[AX1] = 0.0;
  v[AX2] = -B0*cos(phi)/kvec[0];
  v[AX3] = -B0*sin(phi)/kvec[0] + B0*xvec[1];  // Assume rho = 1 
  RotateVector(v + AX1,-1);

  // if (first_call == 1){
  //   print ("o vA    = %18.15e\n",vA);
  //   print ("o omega = %18.15e\n",kmod*vA);
  //   print ("o T     = %18.15e\n",2*CONST_PI/(kmod*vA));
  //   first_call = 0;
  // }

  #if PHYSICS == ResRMHD
  v[EX1] = -(v[VX2]*v[BX3] - v[VX3]*v[BX2]);
  v[EX2] = -(v[VX3]*v[BX1] - v[VX1]*v[BX3]);
  v[EX3] = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]);
  #endif

  #if ENTROPY_SWITCH != NO
  v[ENTR] = Entropy(v,run_config);
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
  -
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
                        boundary box of the corresponding to side
  \param [in]     side  side of computational domain (special value side == 0 
                        is used to control INTERNAL_BOUNDARY)
  \param [in]     grid  pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
}

