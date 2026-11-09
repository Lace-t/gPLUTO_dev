////////////////////////////////////////////////////////////////////////
/*! 
  \file
  \brief Disk-Planet interaction problem.

  Simulate the interaction of a planet embedded in a disk as described
  in section 3.4 of Mignone et al., A&A (2012) 545, A152.
  This test is a nice benchmark for the FARGO module and the
  \c ROTATING_FRAME switch.
  For testing-purposes  no viscosity is used here.
  The initial condition consists of a locally isothermal configuration
  with temperature profile \f$\propto T^{-1}\f$  yielding a disk vertical
  height to radius of \c 0.05.
  The gravitational potential due to the presence of the star and the
  planet is defined in BodyForcePotential() function.

  The conventions used throught the implementation are the following:
 
  - \c r  = spherical radius
  - \c R  = cylindrical radius
  - \c z  = cylindrical height
  - \c th = meridional angle
  
  The test can be carried out in polar (2D or 3D) or spherical (3D)
  coordinates and the following parameters determine the initial configuration:
   
  -# <tt>inputParam[Mstar]</tt>: controls the star mass (in solar masses)
  -# <tt>inputParam[Mdisk]</tt>: controls the disk mass (in solar masses)
  -# <tt>inputParam[Mplanet]</tt>: sets the planet mass (in earth masses)
  -# <tt>inputParam[Viscosity]</tt>: sets the amount of viscosity


  Computation can be carried in the rotating or in the observer's frame
  of reference (\c ROTATING_FRAME to \c YES or \c NO, respectively).
  In particular:

  - Configurations #01 and #02 are in 2D polar coordinates without and with
    FARGO, in the rotating frame.
  - Configurations #03, #04 and #05 are in spherical 3D coordinates with
    and without FARGO in the rotating frame
  - Configurations #06 and #07 are in 2D polar coordinates but in the
    observer's frame
  - Configuration #08 employs static AMR  (grid levels are spatial
    dependent but not dependent on time) in the rotating frame.

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CYLINDRICAL|2  |LINEAR     |RK2        |roe
   #02 |HD  |CYLINDRICAL|2  |LINEAR     |RK2        |roe
   #03 |HD  |SPHERICAL  |3  |LINEAR     |RK2        |hllc
   #04 |HD  |SPHERICAL  |3  |LINEAR     |RK2        |hllc
   #06 |HD  |CYLINDRICAL|2  |LINEAR     |RK2        |roe
  <\CENTER>


  \author A. Mignone (andrea.mignone@unito.it)
  \date   Feb 5, 2025

  \b References:
     - "A Conservative orbital advection scheme for simulations
        of magnetized shear flows with the PLUTO Code"
        Mignone et al, A&A (2012)
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define MIN_DENSITY 1e-8

#if DIMENSIONS == 3
static void NormalizeDensity (const Data *, Grid *, RunConfig &);
#endif

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
*/
////////////////////////////////////////////////////////////////////////
{
  double r, th, R, H, OmegaK;
  double z;
  double scrh, OmegaZ;
  double *inputParam = run_config.inputParam;

  run_config.gamma = 1.01;

  // ----------------------------------------------
  // Define c.g.s units 
  // ----------------------------------------------

  double ul = 5.2*CONST_au;
  double ud = CONST_Msun/(ul*ul*ul);
  double uv = (CONST_G*inputParam[Mstar]*CONST_Msun/ul)/(2.*CONST_PI); 

  uv = sqrt(uv);
  run_config.unit_length   = ul;
  run_config.unit_velocity = uv;
  run_config.unit_density  = ud;

  // ----------------------------------------------
  // Rotating frame speed
  // ----------------------------------------------

  OmegaZ = 0.0;
  #if ROTATING_FRAME == YES
  scrh    = inputParam[Mplanet]/inputParam[Mstar]*CONST_Mearth/CONST_Msun;
  OmegaZ  = sqrt(1.0 + scrh);
  OmegaZ *= 2.0*CONST_PI;
  #endif
  run_config.OmegaZ = OmegaZ;

#if GEOMETRY == CYLINDRICAL
  R  = x1;
  #if DIMENSIONS == 2
  z  = 0.0;
  r  = R;
  th = 0.5*CONST_PI;
  #else
  z  = x3;
  r  = sqrt(R*R + z*z);
  th = atan2(R,z);
  #endif
#elif GEOMETRY == SPHERICAL
  r  = x1;
  th = x2;
  R  = r*sin(th);
  z  = r*cos(th);
#endif
  
  H      = 0.05*R;
  OmegaK = 2.0*CONST_PI/(R*sqrt(R));
  
  scrh    = (0.5*CONST_PI - th)*r/H;
  v[RHO] = 1.0/(R*sqrt(R))*exp(-0.5*scrh*scrh);
  v[VX1] = v[VX2] = v[VX3] = 0.0;

  v[iVPHI] = R*(OmegaK - OmegaZ);
  #ifdef FARGO
  v[FARGO_W] = R*(OmegaK - OmegaZ);
  v[iVPHI]   = 0.0;
  #endif  
  #if EOS == IDEAL
  double cs     = H*OmegaK;
  v[PRS] = v[RHO]*cs*cs;
  #elif EOS == ISOTHERMAL
  run_config.isoSoundSpeed = 0.1*CONST_PI;
  #endif
}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*
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
/* 
*/
////////////////////////////////////////////////////////////////////////
{
  int   i, j, k, nv;
  int   ibeg = grid->lbeg[IDIR];
  int   iend = grid->lend[IDIR];
  double *x1, R;
  double OmegaK, OmegaZ = data->run_config->OmegaZ;
  
  x1 = grid->xC[IDIR];
  #if GEOMETRY == SPHERICAL
  double *x2;
  x2 = grid->xC[JDIR];
  #endif

  #if DIMENSIONS == 3
  static int do_once = 1;
  if (side == 0){
    if (do_once){
      NormalizeDensity(data, grid, *(data->run_config));
      do_once = 0;
    }
  }
  #endif

  if (side == X1_BEG){
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box, i,j,k){
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(2*ibeg - i - 1,j,k,nv);
      data->Vc(i,j,k,VX1) *= -1.0;
      #if GEOMETRY == CYLINDRICAL
      R = x1[i];
      #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      
      data->Vc(i,j,k,iVPHI) = R*(OmegaK - OmegaZ);
      #ifdef FARGO
      data->Vc(i,j,k,iVPHI) = 0.0;
      #endif      
    }
  }

  if (side == X1_END){
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box, i,j,k){
      NVAR_LOOP(nv)  data->Vc(i,j,k,nv) = data->Vc(iend,j,k,nv);
      #if GEOMETRY == CYLINDRICAL
      R = x1[i];
      #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
      data->Vc(i,j,k,VX1) = 0.0;
      data->Vc(i,j,k,VX2) = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      data->Vc(i,j,k,iVPHI) = R*(OmegaK - OmegaZ);
      #ifdef FARGO
      data->Vc(i,j,k,iVPHI) = 0.0;
      #endif      
    }
  }
}

#if DIMENSIONS == 3
////////////////////////////////////////////////////////////////////////
void NormalizeDensity (const Data *data, Grid *grid, RunConfig &run_config)
/*
 * Normalize density and pressure as   rho -> K*rho, where
 *
 *   K = M/(\sum rho*dV)
*/
////////////////////////////////////////////////////////////////////////
{
  int   i, j, k;
  double mc;
  double ul = run_config.unit_length;
  double ud = run_config.unit_density;

  RBox box;

  RBoxDefine (grid->lbeg, grid->lend, &box);

  mc  = 0.5*run_config.inputParam[Mdisk]*CONST_Msun;
  mc /= ud*ul*ul*ul;
  #pragma acc parallel loop collapse(3) present(data) copyin(box)
  BOX_LOOP(box, i,j,k){
    data->Vc(i,j,k,RHO) *= mc;
    #if EOS == IDEAL
    data->Vc(i,j,k,PRS) *= mc;
    #endif
  }
}
#endif

#if (BODY_FORCE & VECTOR)
////////////////////////////////////////////////////////////////////////
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
/*
 *
 *
*/
////////////////////////////////////////////////////////////////////////
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double BodyForcePotential(const State& state, const RunConfig& run_config)
/*
 *
 *
*/
////////////////////////////////////////////////////////////////////////
{
  const double x1 = state.x1;
  const double x2 = state.x2;

  double d, R, r, z, x, y, phiplanet, rsm;
  double xp, yp, phi;
  const double *param = run_config.inputParam;

#if GEOMETRY == CYLINDRICAL

  R  = x1;
  #if DIMENSIONS == 2
  z  = 0.0;
  r  = R;
  // th = 0.5*CONST_PI;
  #else
  const double x3 = state.x3;
  z  = x3;
  r  = sqrt(R*R + z*z);
  // th = atan2(R,z);
  #endif
  
  x  = R*cos(x2);
  y  = R*sin(x2);

#elif (GEOMETRY == SPHERICAL)
  const double x3 = state.x3;
  r  = x1;
  // double th = x2;
  R = r*sin(x2);
  z = r*cos(x2);
  x = r*sin(x2)*cos(x3);
  y = r*sin(x2)*sin(x3);
#endif

  // ---------------------------------------------
  // Planet position
  // ---------------------------------------------

#if ROTATING_FRAME == NO
  double OmegaZ;
  double t = run_config.time;
//  if (g_stepNumber == 2) t += g_dt; // !!!!! THIS IS NOT POSSIBLE IN GPLUTO
  OmegaZ  = sqrt(1.0 + param[Mplanet]/param[Mstar]*CONST_Mearth/CONST_Msun);
  OmegaZ *= 2.0*CONST_PI;

  xp = cos(OmegaZ*t);
  yp = sin(OmegaZ*t);
#else
  xp = 1.0/sqrt(2.0);  // initial planet position 
  yp = 1.0/sqrt(2.0); 
#endif

  d = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z);
  rsm = 0.03*R;
  if (d > rsm) phiplanet = param[Mplanet]/d;
  else phiplanet = param[Mplanet]/d*(pow(d/rsm,4.)-2.*pow(d/rsm,3.)+2.*d/rsm);
  
  phi  = - 4.0*CONST_PI*CONST_PI/param[Mstar];
  phi *= (param[Mstar]/r + phiplanet*CONST_Mearth/CONST_Msun);

  return phi;
}

#endif

