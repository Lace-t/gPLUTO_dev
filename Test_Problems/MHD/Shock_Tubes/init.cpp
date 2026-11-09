////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD shock tube problem.

  Shock tube problems are popular benchmark used to test
  shock-capturing schemes.
  It is a one-dimensional problem with initial condition given
  by a discontinuity separating two constant states.

  Left and right states are provided at runtime using the parameters 
  defined in pluto.ini. 

  Below the pre-defined configurations:

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |DIVB| Solver
  -----|----|-----------|---|-----------|-----------|----|-------
   #01 |MHD |CARTESIAN  |1  |LINEAR     |RK2        |  - | hlld
   #02 |MHD |CARTESIAN  |1  |LINEAR     |RK2        |  - | roe
   #03 |MHD |CARTESIAN  |1  |LINEAR     |RK2        |  - | roe
   #04 |MHD |CARTESIAN  |1  |PARABOLIC  |RK3        |  - | hll
   #05 |MHD |CARTESIAN  |1  |LINEAR     |RK3        |  - | hll
   #06 |MHD |CARTESIAN  |1  |LINEAR     |RK3        |  - | hlld
   #07 |MHD |CARTESIAN  |1  |PARABOLIC  |RK3        |  - | roe
   #08 |MHD |CARTESIAN  |1  |WENO3      |RK3        |  - | gforce
  </CENTER>

  \author A. Mignone (mignone@ph.unito.it)
  \date   17 April 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define LEFTGOING   YES

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
//
 // * Assign initial conditions in terms of primitive variables
// 
////////////////////////////////////////////////////////////////////////
{
#ifdef GAMMA_EOS
run_config.gamma=run_config.inputParam[GAMMA_EOS];
#endif

#if LEFTGOING == YES
  if (x1 < 0.5){

    v[RHO] = run_config.inputParam[RHO_LEFT];
    v[VX1] = run_config.inputParam[VX1_LEFT];
    v[VX2] = run_config.inputParam[VX2_LEFT];
    v[VX3] = run_config.inputParam[VX3_LEFT];
    v[BX1] = run_config.inputParam[BX1_CONST];
    v[BX2] = run_config.inputParam[BX2_LEFT];
    v[BX3] = run_config.inputParam[BX3_LEFT];
    #if EOS == IDEAL
     v[PRS] = run_config.inputParam[PRS_LEFT];
    #endif

  }else{

    v[RHO] = run_config.inputParam[RHO_RIGHT];
    v[VX1] = run_config.inputParam[VX1_RIGHT];
    v[VX2] = run_config.inputParam[VX2_RIGHT];
    v[VX3] = run_config.inputParam[VX3_RIGHT];
    v[BX1] = run_config.inputParam[BX1_CONST];
    v[BX2] = run_config.inputParam[BX2_RIGHT];
    v[BX3] = run_config.inputParam[BX3_RIGHT];
    #if EOS == IDEAL
     v[PRS] = run_config.inputParam[PRS_RIGHT];
    #endif
  }

#else

  if (x1 > 0.5){

    v[RHO] =  run_config.inputParam[RHO_LEFT];
    v[VX1] = -run_config.inputParam[VX1_LEFT];
    v[VX2] =  run_config.inputParam[VX2_LEFT];
    v[VX3] =  run_config.inputParam[VX3_LEFT];
    v[BX1] = -run_config.inputParam[BX1_CONST];
    v[BX2] =  run_config.inputParam[BX2_LEFT];
    v[BX3] =  run_config.inputParam[BX3_LEFT];
    #if EOS == IDEAL
     v[PRS] = run_config.inputParam[PRS_LEFT];
    #endif

  }else{

    v[RHO] =  run_config.inputParam[RHO_RIGHT];
    v[VX1] = -run_config.inputParam[VX1_RIGHT];
    v[VX2] =  run_config.inputParam[VX2_RIGHT];
    v[VX3] =  run_config.inputParam[VX3_RIGHT];
    v[BX1] = -run_config.inputParam[BX1_CONST];
    v[BX2] =  run_config.inputParam[BX2_RIGHT];
    v[BX3] =  run_config.inputParam[BX3_RIGHT];
    #if EOS == IDEAL
    v[PRS] = run_config.inputParam[PRS_RIGHT];
    #endif
  }
#endif

  if (run_config.inputParam[DIVIDE_BY_4PI] > 0.5){
    v[BX1] /= sqrt(4.0*CONST_PI);
    v[BX2] /= sqrt(4.0*CONST_PI);
    v[BX3] /= sqrt(4.0*CONST_PI);
  }

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
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
/*
 *
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
