////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Relativistic shock-cloud interaction.

  Sets initial condition for the relativistic shock-cloud interaction
  problem in 2D or 3D. 
  A high-density spherical cloud with \f$ \rho = 10 \f$ is initially
  located at \f$ (0.8,0) \f$ adjacent to a shock wave located at
  \f$ x=0 \f$.
  The post-shock region (\f$ x < 0 \f$) is at rest with
  \f[
     \begin{array}{lcll}
      (\rho,\, v_z,\, B_z,\, p) &=& (42.5942,\,0,\, -2.12971,\, 127.9483)
      & \qquad\mathrm{Ideal}\quad\mathrm{EOS, [MB06]}
     \\[0.5em]
      (\rho,\, v_z,\, B_z,\, p) &=& (39.5052,\,0,\, 1.9753,\, 129.72386)
      & \qquad\mathrm{TM}\quad\mathrm{EOS, [Mig12]}
    \end{array}
  \f]
  while the upstream region is characterized by
  \f[
      (\rho,\, v_z,\, B_z,\, p) = (1, \, -\sqrt{0.99}, \,0.5,\, 10^{-3})
  \f]
  This test has no input parameters

  <CENTER>
  Conf.|PHYS| GEOMETRY |DIM|RECONSTR.|TIME_STEP|DIVB|Solver|
  -----|----|----------|---|---------|---------|----|------|
   #03 |RMHD|CARTESIAN | 2 |LINEAR   | RK2     | NO |hll   |
   #04 |RMHD|CARTESIAN | 2 |PARABOLIC| RK3     | NO |tvdlf |
  </CENTER>

  \b References
     - [MB06]: "An HLLC Riemann solver for relativistic flows -- 
       II. Magnetohydrodynamics", Mignone & Bodo MNRAS(2006) 368,1040
     - [Mig12]: "The PLUTO Code for Adaptive Mesh Computations in Astrophysical
       Fluid dynamics", Mignone et al., ApJS (2012), 198:7
       
  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Apr 07, 2025
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
  double x,y,z,r;
  
  #if EOS == IDEAL
  run_config.gamma = 4.0/3.0;
  #endif

  x = x1; y = x2; z = x3;

  if (x < 0.6){

    #if EOS == IDEAL
    v[RHO] =   42.5942;
    v[PRS] =  127.9483;
    v[BX3] =   -2.12971;
    #elif EOS == TAUB
    v[RHO] =   39.5052;
    v[PRS] =  129.72386;
    v[BX3] =   1.9753;
    #endif
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[BX1] = 0.0;
    v[BX2] = 0.0;

  }else{

    v[RHO] = 1.0;
    v[PRS] = 1.e-3;
    v[VX1] = -sqrt(1.0 - 1.0/100.);
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.5;
  }

  r = DIM_EXPAND( (x-0.8)*(x-0.8), + y*y, + z*z);
  r = sqrt(r);
  if (r < 0.15) v[RHO] = 10.0;

  #if ENTROPY_SWITCH == YES
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



