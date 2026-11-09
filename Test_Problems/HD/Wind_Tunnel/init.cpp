////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HD wind tunnel test with internal obstacle.

  This test models a two-dimensional hydrodynamic (HD) wind tunnel
  where a uniform supersonic flow interacts with a solid obstacle.

  The computational domain is initialized with constant density and
  pressure, and a uniform flow moving along the x-direction. An internal
  rectangular region acts as a solid obstruction, producing shocks and
  flow deflection.

  Initial conditions:
  - Density: \f$ \rho = 1.4 \f$
  - Pressure: \f$ P = 1.0 \f$
  - Velocity:
    - \f$ v_x = 3.0 \f$ (supersonic inflow)
    - \f$ v_y = v_z = 0 \f$

  Obstacle:
  A rectangular region defined by:
  \f[
    x \ge 0.6,\quad y \le 0.2
  \f]
  where the velocity is set to zero and cells are flagged as
  \c FLAG_INTERNAL_BOUNDARY, mimicking a solid body.

  Boundary conditions:
  - X1_BEG (left boundary): fixed inflow with constant primitive variables
  - Internal region: enforced solid boundary via user-defined flagging

  Physics:
  - Equation of state: ideal gas
  - Adiabatic index: \f$ \gamma = 1.4 \f$

  Purpose:
  This setup is a standard benchmark to test shock formation, reflection,
  and flow interaction with solid obstacles in compressible hydrodynamics.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!

<CENTER>
Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
-----|----|-----------|---|-----------|-----------|-------
 #01 |HD  |CARTESIAN  |2  |LINEAR     |RK2        |two_shock
<\CENTER>

*/
////////////////////////////////////////////////////////////////////////
void Init (double *us, double x1, double x2, double x3, RunConfig &run_config)
{
  run_config.gamma = 1.4;

  us[VX1] = (x1 > 0.6 && x2 < 0.2 ? 0.0:3.0);
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[TRC] = 0.0;

  us[RHO] = 1.4;
  us[PRS] = 1.0;

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
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
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
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double *x = grid->xC[IDIR];
  double *y = grid->xC[JDIR];

  if (side == 0){
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){
      if (y[j] <= 0.2 && x[i] >= 0.6) {
        data->flag(i,j,k) |= FLAG_INTERNAL_BOUNDARY;   
        data->Vc(i,j,k,RHO) = 1.4;
        data->Vc(i,j,k,VX1) = 0.0;
        data->Vc(i,j,k,VX2) = 0.0;
        data->Vc(i,j,k,VX3) = 0.0;
        data->Vc(i,j,k,PRS) = 1.0;
      }
    }
  }
  
  if (side == X1_BEG){
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){
      data->Vc(i,j,k,RHO) = 1.4;
      data->Vc(i,j,k,VX1) = 3.0;
      data->Vc(i,j,k,VX2) = 0.0;
      data->Vc(i,j,k,VX3) = 0.0;
      data->Vc(i,j,k,PRS) = 1.0;
    }
  } 
}
