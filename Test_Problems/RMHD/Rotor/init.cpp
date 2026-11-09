////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief RMHD rotor test problem in 2 or 3 dimensions.

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
  where \f$\omega = 9.95\f$ is the angular frequency of rotation.
  The computational domain is the unit square and outflow boundary conditions
  are imposed everywhere.

  As the disk rotates, strong torsional Alfven waves form and propagate
  outward carrying angular momentum from the disk to the ambient.
  The emerging flow structure is enclosed by a circular fast forward shock
  traveling into the surrounding medium.
  An inward fast shock bounds the innermost oval region where density has
  been depleted to lower values.
  The presence of the magnetic field slows down the rotor, and the maximum
  Lorentz factor decreases from the nominal value of 10 to 2.2 (approx).

  A list of tested configurations is given in the following table:

  <CENTER>
  Conf.| GEOMETRY  | DIM |T. STEPPING|RECONSTRUCTION| divB |  HO
  -----|-----------|-----|-----------|--------------|------|-----
  #01  | CARTESIAN |  2  |  RK3      |  LINEAR      | CT   |  NO
  #02  | CARTESIAN |  2  |  RK2      |  LINEAR      | CT   |  NO
  #03  | CARTESIAN |  2  |  RK2      |  LINEAR      | 8W   |  NO
  #04  | CARTESIAN |  2  |  RK2      |  LINEAR      | CT   |  NO
  #05  | CARTESIAN |  3  |  RK2      |  LINEAR      | GLM  |  NO
  #07  | CARTESIAN |  3  |  RK2      |  LINEAR      | CT   |  NO
  #08  | CARTESIAN |  3  |  RK3      |  LimO3       | GLM  |  NO
  #10  | CARTESIAN |  3  |  RK2      |  LINEAR      | CT   |  NO
  #11  | CARTESIAN |  2  |  RK3      |  WENO3       | CT   |  NO
  #12  | CARTESIAN |  2  |  RK3      |  WENO3       | GLM  | YES
  #13  |CYLINDRICAL|  2  |  RK2      |  LINEAR      | 8W   |  NO
  #14  |CYLINDRICAL|  2  |  RK3      |  PARABOLIC   | CT   |  NO
  </CENTER>

  The final solution at t = 0.4 on a grid of 400x400 grid points is shown below.

  The three dimensional version extends the current configuration to a
  spinning sphere and it is described in [MUB09].

  \b References:
     - [dZBL03] "An efficient shock-capturing central-type scheme for 
       multidimensional relativistic flows"
        Del Zanna, Bucciantini, Londrillo, A&A (2003) 400, 397
     - [MUB09] "A five-wave Harten-Lax-van Leer Riemann solver for relativistic
        magnetohydrodynamics"
        Mignone, Ugliano \& Bodo, MNRAS (2009) 393, 1141.
     - [Mig12] "The PLUTO Code for Adaptive Mesh Computations 
        in Astrophysical Fluid Dynamics" Mignone et al., 
        ApJS (2012) 198, 7 (see Sect. 6.3)

  \author A. Mignone (andrea.mignone@unito.it)

  \date   Feb 24, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
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
  double r, r0, Bx, omega;

  omega = 0.995;
  Bx    = 1.0;
  r0    = 0.1;

  #if GEOMETRY == CARTESIAN

  r = DIM_EXPAND(x1*x1, + x2*x2, + x3*x3);
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
    v[VX1] = -omega*x2/r0;
    v[VX2] =  omega*x1/r0;
  }

  v[AX1] = v[AX2] = 0.0;
  v[AX3] = v[BX1]*x2;

  #elif GEOMETRY == CYLINDRICAL

  r = x1;

  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[BX1] =  Bx*cos(x2);
  v[BX2] = -Bx*sin(x2);
  v[BX3] = 0.0;
  v[PRS] = 1.0;

  if (r <= r0) {
    v[RHO] = 10.0;
    v[VX2] = omega*r/r0;
  }

  v[AX1] = v[AX2] = 0.0;
  v[AX3] = - r*sin(x2)*Bx;
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
  int     i, j, k, nv;
  int     ibeg = grid->lbeg[IDIR];
  double  *r   = grid->xC[IDIR];

  if (side == X1_BEG){
    #pragma acc parallel loop collapse(3) present(data,box)
    BOX_LOOP(*box,i,j,k){
      double slp = r[i]/r[ibeg];
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(ibeg,j,k,nv);

      data->Vc(i,j,k,VX1) *= slp;
      data->Vc(i,j,k,VX2) *= slp;
    }
    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(data,box)
    BOXI_LOOP(*box,i,j,k) data->Vs[BX2s](i,j,k) = data->Vs[BX2s](ibeg,j,k);
    #endif
  }
}
