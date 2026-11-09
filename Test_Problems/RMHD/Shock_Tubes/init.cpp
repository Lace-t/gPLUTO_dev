#include "gpluto.hpp"
////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Collection of relativistic MHD shock-tube problems.

  This directory contains several 1D shock-tube configurations
  commonly used for benchmarking of numerical methods.
  There are 16  input parameters that can be used to specify entirely
  the left and right states with respect to the discontinuity initially
  placed at \c x=1/2.
  For the left state one has:

  - <tt>inputParam[RHO_LEFT]</tt>: density for the left state
  - <tt>inputParam[VX_LEFT]</tt>: normal (x) velocity for the left state
  - <tt>inputParam[VY_LEFT]</tt>: tangential (y) velocity for the left state
  - <tt>inputParam[VZ_LEFT]</tt>: tangential (z) velocity for the left state
  - <tt>inputParam[BY_LEFT]</tt>: tangential (y) magnetic field for the left state
  - <tt>inputParam[BZ_LEFT]</tt>: tangential (z) magnetic field for the left state
  - <tt>inputParam[PR_LEFT]</tt>: pressure in the left state

  Similar parameters are used to specify the right states.
  In addition,

  - <tt>inputParam[BX_CONST]</tt>: normal component of magnetic field (this
    cannot have jump)
  - <tt>inputParam[GAMMA_EOS]</tt>: the specific heat ratio

  The IDEAL equation of state is used.
  The available configurations are taken from

  <CENTER>
  Conf.|PHYS| GEOMETRY |DIM|   REC   |TIME_STEP|DIVB|Solver|Reference
  -----|----|----------|---|---------|---------|----|------|---------
   #01 |RMHD|CARTESIAN | 1 |LINEAR   | RK3     | NO |hll   |1st prob in [MB06]
   #02 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hllc  |1st prob in [MB06]
   #03 |RMHD|CARTESIAN | 1 |LINEAR   | RK3     | NO |hll   |2nd prob in [MB06]
   #04 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hllc  |2nd prob in [MB06]
   #05 |RMHD|CARTESIAN | 1 |LINEAR   | RK3     | NO |hllc  |3rd prob in [MB06]
   #06 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hllc  |3rd prob in [MB06]
   #07 |RMHD|CARTESIAN | 1 |LINEAR   | RK3     | NO |hllc  |4th prob in [MB06]
   #08 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hllc  |4th prob in [MB06]
   #09 |RMHD|CARTESIAN | 1 |LINEAR   | RK3     | NO |hll   |2nd prob in [MUB09]
   #10 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hlld  |2nd prob in [MUB09]
   #11 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hlld  |4th prob in [MUB09]
   #12 |RMHD|CARTESIAN | 1 |PARABOLIC| RK3     | NO |tvdlf |4th prob in [MUB09]
   #13 |RMHD|CARTESIAN | 1 |LINEAR   | RK2     | NO |hllc  |[MMB05]
   #14 |RMHD|CARTESIAN | 1 |LimO3    | RK3     | NO |hllc  |[MMB05]
  </CENTER>

  \b References
     - [MB06]: "An HLLC Riemann solver for relativistic flows --
       II. Magnetohydrodynamics", Mignone & Bodo MNRAS(2006) 368,1040
     - [MUB09]: "A five-wave Harten-Lax-van Leer Riemann solver for
       relativistic magnetohydrodynamics", Mignone et al., MNRAS(2009) 393,1141
     - [MMB05]: "Relativistic MHD simulations of jets with toroidal
       magnetic fields", Mignone, Massaglia & Bodo, SSRv (2005), 121,21
     - Balsara, ApJS (2001), 132, 83
     - Del Zanna, Bucciantini, Londrillo, A&A (2003)

  \authors A. Mignone (andrea.mignone@unito.it)\n

  \date    Feb 27, 2025
*/
////////////////////////////////////////////////////////////////////////

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
  double *inputParam = run_config.inputParam;
  
  run_config.gamma = run_config.inputParam[GAMMA_EOS];

  if (x1 < 0.5){
    v[RHO] = inputParam[RHO_LEFT];
    v[VX1] = inputParam[VX_LEFT];
    v[VX2] = inputParam[VY_LEFT];
    v[VX3] = inputParam[VZ_LEFT];
    v[BX1] = inputParam[BX_CONST];
    v[BX2] = inputParam[BY_LEFT];
    v[BX3] = inputParam[BZ_LEFT];
    v[PRS] = inputParam[PR_LEFT];
  }else{
    v[RHO] = inputParam[RHO_RIGHT];
    v[VX1] = inputParam[VX_RIGHT];
    v[VX2] = inputParam[VY_RIGHT];
    v[VX3] = inputParam[VZ_RIGHT];
    v[BX1] = inputParam[BX_CONST];
    v[BX2] = inputParam[BY_RIGHT];
    v[BX3] = inputParam[BZ_RIGHT];
    v[PRS] = inputParam[PR_RIGHT];
  }

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*!
  [Not Used]
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*!
  [Not Used]
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
