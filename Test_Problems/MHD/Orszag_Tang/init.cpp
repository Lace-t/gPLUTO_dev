////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Orszag-Tang MHD vortex.

  The Orszag Tang vortex system describes a doubly periodic fluid 
  configuration leading to two-dimensional supersonic MHD turbulence.
  Although an analytical solution is not known, its simple and reproducible 
  set of initial conditions has made it a widespread benchmark for 
   inter-scheme comparison.

  The computational domain is the unit periodic box \f$[0,1]^D\f$ where 
  \c D is the number of spatial dimensions.
  In 2D, the initial condition is given by
  \f[
     \mathbf{v} = \left(-2\pi\sin y,\, 2\pi\sin x, 0\right) \,,\qquad
     \mathbf{B} = \left(-2\pi\sin y,\, 2\pi\sin 2x, 0\right) \,,\qquad
     \rho = 25/9\,,\qquad
     p    = 5/3
  \f]
  
  This test problem does not have any input parameter.
  Below the pre-defined configurations:

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM|   REC     |TIME_STEP  |DIVB| Solver
  -----|----|-----------|---|-----------|-----------|----|-------
   #01 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | 8W | tvdlf
   #02 |MHD |CARTESIAN  | 2 | PARABOLIC | RK3       | 8W | roe
   #03 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       |GLM | hlld
   #04 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | 8W | roe
   #05 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | tvdlf
   #06 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | roe
   #07 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | hlld
   #08 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | hlld
   #09 |MHD |CARTESIAN  | 2 | LimO3     | RK3       | CT | gforce
   #10 |MHD |CARTESIAN  | 2 | LINEAR    | RK4 [HO]  | CT | roe
   #11 |MHD |CARTESIAN  | 3 | WENOZ     | RK3       | CT | hlld
  </CENTER>

  \author A. Mignone (andrea.mignone@unito.it)

  \date   April 15, 2026

  \b References
     - "Comparison of some Flux Corrected Transport and TVD numerical 
        schemes for hydrodynamic and magnetohydrodynamic problems",
        Toth & Odstrcil, JCP (1996) 128, 82
     - "High-order conservative finite difference GLM-MHD schemes for 
        cell-centered MHD", Mignone, Tzeferacos & Bodo, JCP (2010) 229, 5896.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifndef PLANE_ROTATION 
 #define PLANE_ROTATION NO
#endif

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
//
   // Assign initial conditions in terms of primitive variables
// 
////////////////////////////////////////////////////////////////////////
{
  double x = x1, y = x2;
  double z = x3;

  // Shift coordinates when rotating plane

  #if PLANE_ROTATION == NO

  #elif PLANE_ROTATION == 1
  x = x2; y = x3; z = x1;   // yz plane 
  #elif PLANE_ROTATION == 2
  x = x3; y = x1; z = x2;   // zx plane 
  #else   
    #error Invalid PLANE_ROTATION
  #endif

  // Assign initial conditions using xy plane vectors

  v[RHO]=  25./9.0;

  #if (DIMENSIONS == 2) || (PLANE_ROTATION != 0)
  v[VX1]= -sin(2.0*M_PI*y);
  v[VX2]=  sin(2.0*M_PI*x);
  v[VX3]=  0.0;
  #elif (DIMENSIONS == 3) && (PLANE_ROTATION == 0) 

  // Taken from Helzel, Rossmanith 2011 

  double eps = 0.2;
  double fact = (1.0 + eps*sin(2.0*M_PI*z));
  v[VX1]= -sin(2.0*M_PI*y)*fact;
  v[VX2]=  sin(2.0*M_PI*x)*fact;
  v[VX3]=  eps*sin(2.0*M_PI*z);
  #endif

  #if PHYSICS == MHD
  v[BX1]= -sin(2.0*M_PI*y);
  v[BX2]=  sin(4.0*M_PI*x);
  v[BX3]=  0.0;
  #endif

  v[PRS]=  5.0/3.0;

  #ifdef GLM_MHD
  v[PSI_GLM] = 0.0;
  #endif
  #if NTRACER == 1
  v[TR1]=  y > 0.5 ? 1.0:0.0;
  #endif

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = (cos(2.0*M_PI*y) + 0.5*cos(4.0*M_PI*x))/(2.0*M_PI);

  // Shift vectors
  
  #if PLANE_ROTATION == 1  // yz plane 
  { 
    int nv;
    double v0[32];
    for (nv = 0; nv < 32; nv++) v0[nv] = v[nv];

    v[VX1] = v0[VX3]; v[VX2] = v0[VX1]; v[VX3] = v0[VX2];
    v[BX1] = v0[BX3]; v[BX2] = v0[BX1]; v[BX3] = v0[BX2];
    v[AX1] = v0[AX3]; v[AX2] = v0[AX1]; v[AX3] = v0[AX2];
  }
  #elif PLANE_ROTATION == 2  // zx plane
  {
    double v0[NVAR];
    NVAR_LOOP(nv) v0[nv] = v[nv];

    v[VX1] = v0[VX1]; v[VX2] = v0[VX3]; v[VX3] = v0[VX1];

    int nv;
    double v0[32];
    for (nv = 0; nv < 32; nv++) v0[nv] = v[nv];

    v[VX1] = v0[VX2]; v[VX2] = v0[VX3]; v[VX3] = v0[VX1];
    v[BX1] = v0[BX2]; v[BX2] = v0[BX3]; v[BX3] = v0[BX1];
    v[AX1] = v0[AX2]; v[AX2] = v0[AX3]; v[AX3] = v0[AX1];
  }
  #endif
  
  // Compute entropy

  #if ENTROPY_SWITCH != NO
  v[ENTR] = Entropy(v,run_config);
  #else
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
  Perform runtime data analysis.
  
  \param [in] d the PLUTO Data structure
  \param [in] grid   pointer to array of Grid structures
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
                        boundary box corresponding to side
  \param [in]     side  side of computational domain (special value side == 0 
                        is used to control INTERNAL_BOUNDARY)
  \param [in]     grid  pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
}

#if BODY_FORCE != NO
////////////////////////////////////////////////////////////////////////
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
////////////////////////////////////////////////////////////////////////
{ 
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

////////////////////////////////////////////////////////////////////////
double BodyForcePotential(const State& state, const RunConfig& run_config)
////////////////////////////////////////////////////////////////////////
{
  return 0.0;
}
#endif
