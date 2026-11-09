////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Magnetic field diffusion in 2D and 3D.

  Sets the initial conditions for a magnetic field diffusion problem
  in 2 or 3 dimensions.
  This is a useful test to check the ability of the code to solve 
  standard diffusion problems.
  The magnetic field has initially a Gaussian profile, and an anisotropic
  resistivity is possible.
  This problem has an analytical solution given, in 2D, by
  \f[
    B_x(y,t) = \exp(-y^2/4\eta_zt)/\sqrt{t}
    \quad\quad\quad
    B_y(x,t) = \exp(-x^2/4\eta_zt)/\sqrt{t}
    \quad\quad\quad
    B_z(x,y,t) = \exp(-x^2/4\eta_yt)\exp(-y^2/4\eta_xt)/t
  \f]
  and in 3D by
  \f[ 
    \begin{array}{lcl}
    B_x(y,z,t) &=& \exp(-y^2/4\eta_zt)\exp(-z^2/4\eta_yt)/t \\[0.5em]
    B_y(x,z,t) &=& \exp(-x^2/4\eta_zt)\exp(-z^2/4\eta_xt)/t \\[0.5em]
    B_z(x,y,t) &=& \exp(-x^2/4\eta_yt)\exp(-y^2/4\eta_xt)/t
    \end{array}
  \f]
  The initial condition is simply set using the previous profiles 
  with \c t=1.
  In order to solve only the parabolic term in the induction equation
  (\f$ \nabla\times\vec{J}=\nabla\times(\eta\nabla\times\vec{B})\f$) 
  we give to the fluid a very large intertia using a high density 
  value.
  Moreover, to avoid any fluid motion, the velocity is reset to zero 
  at each time step by using the \c INTERNAL_BOUNDARY.

  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[ETAX]</tt>: sets the resistivity along the \f$x\f$ direction;
  - <tt>g_inputParam[ETAY]</tt>: sets the resistivity along the \f$y\f$ direction;
  - <tt>g_inputParam[ETAZ]</tt>: sets the resistivity along the \f$z\f$ direction;

  The configurations use both \c EXPLICIT and \c STS time integrators and
  different geometries are explored:

  <CENTER>
  Conf.|GEOMETRY    | DIM | T.STEPPING | divB | RESISTIVITY
  -----|----------  |-----|------------|------| -----------
   #01 |CARTESIAN   |  3  | RK2        |  8W  | EXPLICIT
   #02 |CARTESIAN   |  3  | HANCOCK    |  GLM | STS       WiP
   #03 |CARTESIAN   |  3  | RK3        |  CT  | EXPLICIT
   #04 |CARTESIAN   |  3  | HANCOCK    |  GLM | EXPLICIT  WiP
   #05 |CYLINDRICAL |  3  | RK2        |  8W  | EXPLICIT
   #06 |CYLINDRICAL |  3  | RK2        |  8W  | STS       WiP
   #07 |SPHERICAL   |  3  | RK2        |  8W  | EXPLICIT
   #08 |SPHERICAL   |  3  | RK2        |  8W  | STS       WiP
   #09 |SPHERICAL   |  3  | RK3        |  CT  | EXPLICIT
   #10 |CARTESIAN   |  2  | HANCOCK    |  CT  | EXPLICIT  WiP
   #11 |CARTESIAN   |  3  | HANCOCK    |  CT  | EXPLICIT  WiP
   #12 |CARTESIAN   |  2  | HANCOCK    |  CT  | STS       WiP
   #13 |SPHERICAL   |  3  | RK2        |  CT  | STS       WiP
   #14 |CARTESIAN   |  3  | HANCOCK    |  GLM | RKL       WiP
   #15 |SPHERICAL   |  3  | RK2        |  GLM | RKL       WiP
  </CENTER>

  ALL THE MISSING CONFIGURATIONS WILL BE ADDED LATER

  \authors A. Mignone (andrea.mignone@unito.it) \n
           T. Matsakos (titos@oddjob.uchicago.edu) \n
           D. Crocco
  \date    December 02, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

void BoundValues (double *v, double x1, double x2, double x3, double t, RunConfig &run_config);

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  BoundValues(v, x1, x2, x3, 1.0, run_config);
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
// 
// 
////////////////////////////////////////////////////////////////////////
{
}

#if PHYSICS == MHD
////////////////////////////////////////////////////////////////////////
/*!
  Define the component of a static, curl-free background 
  magnetic field.
  
  \param [out] state        reference to a State structure containing a 
  valid pointer to Bgf
  \param [in]  run_config   reference to RunConfig structure
*/
////////////////////////////////////////////////////////////////////////
void BackgroundField (State &state, RunConfig &run_config)
{
  state.Bgf[IDIR] = 0.0;
  state.Bgf[JDIR] = 0.0;
  state.Bgf[KDIR] = 0.0;

}
#endif

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
//
 // *                                                                    
////////////////////////////////////////////////////////////////////////
{
  #pragma acc update host(data, grid, box)  // Ensure that the device has the latest data before starting the computation 
  int  i, j, k, nv;
  double *x1, *x2, *x3;
  double t, vb[NVAR];

  if (side == 0){  
    #if GEOMETRY == CARTESIAN
    #pragma acc parallel loop collapse(3) present(data, grid, box)
    BOX_LOOP((*box),i,j,k){
      data->Vc(i,j,k,VX1) = 0.0;
      data->Vc(i,j,k,VX2) = 0.0;
      data->Vc(i,j,k,VX3) = 0.0;
      #if ENTROPY_SWITCH == YES
      data->flag(i,j,k) |= FLAG_ENTROPY;
      #endif
    }

    PrimToCons3D(data,0,grid,box);

    #endif
  }
  #pragma acc update host(data, grid, box)

  t = data->run_config->time + 1.0;

  x1 = grid->xC[IDIR]; 
  x2 = grid->xC[JDIR]; 
  x3 = grid->xC[KDIR]; 

  #pragma acc update host(data, grid, box)  // Ensure that the device has the latest data before starting the computation 
  if (side == X1_BEG || side == X2_BEG || side == X3_BEG ||
      side == X1_END || side == X2_END || side == X3_END){  // -- All Boundaries -- 
    #pragma acc parallel loop collapse(3) present(data, grid, box,data->run_config) private (vb[:NVAR])
    BOX_LOOP((*box),i,j,k){
      BoundValues (vb, x1[i], x2[j], x3[k], t, *(data->run_config));
      #pragma acc loop seq
      for (nv = 0; nv < NVAR; nv++) data->Vc(i,j,k,nv) = vb[nv];
    }

    #ifdef STAGGERED_MHD
    double x1p = grid->xL[IDIR];
    double x2p = grid->xL[JDIR];
    double x3p = grid->xL[KDIR];

    if (side == X2_BEG || side == X2_END || side == X3_BEG || side == X3_END){
      #pragma acc parallel loop collapse(3) present(data, grid, box,data->run_config) private (vb[:NVAR])
      BOXI_LOOP((*box),i,j,k){
        BoundValues (vb, x1p[i], x2[j], x3[k], t, *(data->run_config));
        data->Vs[BX1s](i,j,k) = vb[BX1];
      }
    }
    if (side == X1_BEG || side == X1_END || side == X3_BEG || side == X3_END){
      #pragma acc parallel loop collapse(3) present(data, grid, box,data->run_config) private (vb[:NVAR])     
      BOXJ_LOOP((*box),i,j,k){
        BoundValues (vb, x1[i], x2p[j], x3[k], t, *(data->run_config));
        data->Vs[BX2s](i,j,k) = vb[BX2];
      }
    }
    if (side == X1_BEG || side == X1_END || side == X2_BEG || side == X2_END){
    #if (DIMENSIONS == 3)
      #pragma acc parallel loop collapse(3) present(data, grid, box,data->run_config) private (vb[:NVAR])
      BOXK_LOOP((*box),i,j,k){
        BoundValues (vb, x1[i], x2[j], x3p[k], t, *(data->run_config));
        data->Vs[BX3s](i,j,k) = vb[BX3];
      }
    #endif
    }
    #endif
  }

  #pragma acc update host(data, grid, box)
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void BoundValues (double *v, double x1, double x2, double x3, double t, RunConfig &run_config)
//
 // * Time-dependent exact solution
// 
////////////////////////////////////////////////////////////////////////
{
  //#pragma acc update host(run_config,run_config.inputParam,v,x1,x2,x3,t)  /* Ensure that the device has the latest data before starting the computation */

  double Bx, By, Bz;
  double x, y, z;
  double eta_x, eta_y, eta_z;

  double *inputParam = run_config.inputParam;

  eta_x = inputParam[ETAX];
  eta_y = inputParam[ETAY];
  eta_z = inputParam[ETAZ];

// -- find Cartesian coordinates from (x1,x2,x3) -- 

  #if GEOMETRY == CARTESIAN
   x = x1; y = x2; z = x3;
  #elif GEOMETRY == CYLINDRICAL
   x = x1*cos(x2) - 5.0; 
   y = x1*sin(x2) - 5.0; 
   z = x3 - 5.0;
  #elif GEOMETRY == SPHERICAL
   x = x1*cos(x3)*sin(x2) - 5.0;
   y = x1*sin(x3)*sin(x2) - 5.0;
   z = x1*cos(x2) - 5.0;
  #else
   #error geometry not valid
  #endif

  v[RHO] = 1.0e9;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if EOS != ISOTHERMAL
   v[PRS] = 1.0;
  #endif
  #ifdef GLM_MHD
   v[PSI_GLM] = 0.0;
  #endif
 
  #if DIMENSIONS == 2 
   double st = sqrt(t);

   Bx = exp(-0.25*y*y/(eta_z*t))/st;
   By = exp(-0.25*x*x/(eta_z*t))/st;
   Bz = exp(-0.25*x*x/(eta_y*t))*exp(-0.25*y*y/(eta_x*t))/t;
  #elif DIMENSIONS == 3
   Bx = exp(-0.25*y*y/(eta_z*t))*exp(-0.25*z*z/(eta_y*t))/t;
   By = exp(-0.25*x*x/(eta_z*t))*exp(-0.25*z*z/(eta_x*t))/t;
   Bz = exp(-0.25*x*x/(eta_y*t))*exp(-0.25*y*y/(eta_x*t))/t;
  #endif

  // -----------------------------------------------------------
  // Transform back to original coordinate system
  // -----------------------------------------------------------
  
  #if GEOMETRY == CARTESIAN
   v[BX1] = Bx;
   v[BX2] = By;
   v[BX3] = Bz;
  #elif GEOMETRY == CYLINDRICAL
   v[BX1] =  cos(x2)*Bx + sin(x2)*By;
   v[BX2] = -sin(x2)*Bx + cos(x2)*By;
   v[BX3] =  Bz;
  #elif GEOMETRY == SPHERICAL
   v[BX1] =  cos(x3)*sin(x2)*Bx + sin(x3)*sin(x2)*By + cos(x2)*Bz;
   v[BX2] =  cos(x3)*cos(x2)*Bx + sin(x3)*cos(x2)*By - sin(x2)*Bz;
   v[BX3] = -sin(x3)*Bx + cos(x3)*By;
  #else
   QuitPLUTO("! BoundValues: GEOMETRY not defined\n", 1);
  #endif  

  //#pragma acc update device(run_config,v)  /* Ensure that the device has the latest data before starting the computation */
}

#if BODY_FORCE != NO
////////////////////////////////////////////////////////////////////////
/*!
  Prescribe the acceleration vector as a function of the coordinates
  and the vector of primitive variables *v.
  
  \param [in] v  pointer to a cell-centered vector of primitive 
  variables
  \param [out] g acceleration vector
  \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
  \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
  \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
////////////////////////////////////////////////////////////////////////
/*!
  Return the gravitational potential as function of the coordinates.
  
  \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
  \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
  \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
  
  \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double BodyForcePotential(const State& state, const RunConfig& run_config)
{
  return 0.0;
}
#endif
