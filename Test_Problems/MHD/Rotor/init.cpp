////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD Rotor test problem.

  The rotor problem consists of a rapidly spinning cylinder embedded in 
  a static background medium with uniform density and pressure and
  a constant magnetic field along the x direction 
  \f$B_x = 5/\sqrt{4\pi}\f$.
  The cylinder rotates uniformly with constant angular velocity 
  \f$ \omega \f$ and has larger density:
  \f[
     (\rho, v_\phi) = \left\{\begin{array}{ll}
     (10, \omega r) & \qquad {\rm for}\quad r < r_0 
     \\[0.5em]
     (1+9f,f\omega r_0) & \qquad {\rm for}\quad r_0 \le r \le r_1
     \\[0.5em]
     (1,0) & \qquad {\rm otherwise}
   \end{array}\right.
  \f]
  Here \f$ r_0=0.1\f$ and \f$ = (r_1-r)/(r_1-r_0)\f$ is a taper function.
  The ideal equation of state with \f$\Gamma = 1.4\f$ is used.
  As the disk rotates, strong torsional Alfven waves form and propagate 
  outward carrying angular momentum from the disk to the ambient.

  A list of tested configurations is given in the following table:

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM|   REC   |TIME_STEP|DIVB| Solver |BCK FIELD
  -----|----|-----------|---|---------|---------|----|--------|---------
   #01 |MHD |CARTESIAN  | 2 |LINEAR   |RK2      | CT | hlld   |  -
   #02 |MHD |CYLINDRICAL| 2 |LINEAR   |RK2      | CT | hlld   |  -
   #03 |MHD |CYLINDRICAL| 2 |LINEAR   |RK2      | 8W | roe    |  -
   #04 |MHD |CYLINDRICAL| 2 |LINEAR   |RK2      | CT | hll    | YES
   #05 |MHD |CYLINDRICAL| 2 |PARABOLIC|RK3      | GLM| tvdlf  |  -
   #08 |MHD |CARTESIAN  | 2 |LINEAR   |RK2      | CT | hllc   |  -
   #09 |MHD |CYLINDRICAL| 2 |LINEAR   |RK2      | CT | hlld   | YES
   #11 |MHD |CYLINDRICAL| 3 |LINEAR   |RK2      | GLM| hlld   | YES
   #12 |MHD |CARTESIAN  | 3 |WENO3    |RK3      | GLM| hlld   |  -
  </CENTER>

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 3, 2025

  \b Reference: 
     - "On the divergence-free condition in Godunov-type schemes 
        for ideal MHD: the upwind constrained transport method"
        P. Londrillo, L. Del Zanna,  JCP (2004), 195, 17
     - "PLUTO: A numerical code for computational astrophysics"
        Mignone et al., ApJS (2007) 170, 228 (see Sect. 5.4)
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*
*/
////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
{
  double r, r0, r1, Bx, f, omega;

  #ifndef GAMMA_EOS
  run_config.gamma=1.4;
  #endif 
  
#if HAVE_ENERGY
  v[PRS]    = 1.0;
#endif

  omega = 20.0;
  Bx    = 5.0/sqrt(4.0*CONST_PI);
  r0    = 0.1;
  r1    = 0.115;

  #if GEOMETRY == CARTESIAN
   r = sqrt(x1*x1 + x2*x2);

   v[BX1] = Bx;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   f = (r1 - r)/(r1 - r0);
   if (r <= r0) {
     v[RHO] = 10.0;
     v[VX1] = -omega*x2;
     v[VX2] =  omega*x1;
   }else if (r < r1) {
     v[RHO] = 1.0 + 9.0*f;
     v[VX1] = -f*omega*x2*r0/r;
     v[VX2] =  f*omega*x1*r0/r;
   } else {
     v[RHO] = 1.0;
     v[VX1] = 0.0;
     v[VX2] = 0.0;
   }
   v[VX3] = 0.0;

   v[AX1] = v[AX2] = 0.0;
   v[AX3] = Bx*x2;
  #elif GEOMETRY == CYLINDRICAL
   r = x1;

   v[BX1] =  Bx*cos(x2);
   v[BX2] = -Bx*sin(x2);
   v[BX3] = 0.0;
   v[VX1] = 0.0;

  f = (r1 - r)/(r1 - r0);
   if (r <= r0) {
     v[RHO] = 10.0;
     v[VX2] = omega*r;
   }else if (r < r1) {
     v[RHO] = 1.0 + 9.0*f;
     v[VX2] = f*omega*r0;
   } else {
     v[RHO] = 1.0;
     v[VX2] = 0.0;
   }
   v[VX3] = 0.0;
  
   v[AX1] = v[AX2] = 0.0;
   v[AX3] = - r*sin(x2)*Bx;
  #endif

  #if BACKGROUND_FIELD == YES  // With background splitting, the initial 
                                  // condition should contain deviation only 
   v[BX1] = v[BX2] = v[BX3] =
   v[AX1] = v[AX2] = v[AX3] = 0.0;
  #endif

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
//
  // Assign initial condition by looping over the computational domain.
  // Called after the usual Init() function to assign initial conditions
  // on primitive variables.
  // Value assigned here will overwrite those prescribed during Init().
// 
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

#if PHYSICS == MHD
////////////////////////////////////////////////////////////////////////
/*
  Define the component of a static, curl-free background 
  magnetic field.
  
  \param [out] state        reference to a State structure containing a 
  valid pointer to Bgf
  \param [in]  run_config   reference to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
void BackgroundField (State &state, RunConfig &run_config)
{
  const double Bx = 5.0/sqrt(4.0*CONST_PI);

  #if GEOMETRY == CARTESIAN
  state.Bgf[IDIR] = Bx;
  state.Bgf[JDIR] = 0.0;
  state.Bgf[KDIR] = 0.0;
  #elif GEOMETRY == CYLINDRICAL
  const double x2 = state.x2;
  state.Bgf[IDIR] =   Bx*cos(x2);
  state.Bgf[JDIR] = - Bx*sin(x2);
  state.Bgf[KDIR] = 0.0;
  #endif

}
#endif

////////////////////////////////////////////////////////////////////////
/*
  Provide inner radial boundary condition in polar geometry.
  Zero gradient is prescribed on density, pressure and magnetic field.
  For the velocity, zero gradient is imposed on v/r (v = vr, vphi).
*/
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double  *r, slp;
  int  ibeg = grid->lbeg[IDIR];

  if (side == X1_BEG){
    r = grid->xC[IDIR];
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){
      slp = r[i]/r[ibeg];
      data->Vc(i,j,k,RHO) = data->Vc(ibeg,j,k,RHO);
      data->Vc(i,j,k,VX1) = slp*data->Vc(ibeg,j,k,VX1);
      data->Vc(i,j,k,VX2) = slp*data->Vc(ibeg,j,k,VX2);
      data->Vc(i,j,k,VX3) = 0.0;
      #if HAVE_ENERGY
      data->Vc(i,j,k,PRS) = data->Vc(ibeg,j,k,PRS);
      #endif
      data->Vc(i,j,k,BX1) = data->Vc(ibeg,j,k,BX1);
      data->Vc(i,j,k,BX2) = data->Vc(ibeg,j,k,BX2);
      data->Vc(i,j,k,BX3) = 0.0;
      #ifdef GLM_MHD 
      data->Vc(i,j,k,PSI_GLM) = 0.0;
      #endif
    }
    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(data, box)
    BOXJ_LOOP(*box, i,j,k){
      data->Vs[BX2s](i,j,k) = data->Vs[BX2s](ibeg,j,k);
    } 
    #endif
  }
}


////////////////////////////////////////////////////////////////////////
/*! 
*/
////////////////////////////////////////////////////////////////////////
#if MAPPED_GRID == YES
double MappedGrid (Grid *grid, int dir, double xi)
{
  if (dir == IDIR) return MappedGrid_Sinh(grid, dir, xi, 3);
  if (dir == JDIR) return MappedGrid_Sin(grid, dir, xi);
  if (dir == KDIR) return MappedGrid_Uniform(grid, dir, xi);
  return -1.0;
}
#endif
