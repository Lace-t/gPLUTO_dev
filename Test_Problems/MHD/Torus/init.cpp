////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Magnetized accretion torus.

  Magnetized accretion torus in 2D cylindricsal / spherical coordinates 
  in a radially stratified atmosphere. 
  The accretion torus is determined by the equilibrium condition
  \f[
     -\frac{1}{R} + \frac{1}{2(1-a)}\frac{L_k^2}{r^{2-2a}} 
     +(n+1)\frac{p}{\rho} = const
  \f]
  where \f$p = K\rho^\gamma \f$, \c R is the spherical radius, \c r 
  is the cylindrical radius.
  Since \f$ n = 1/(\gamma-1) \f$ then \f$ 1+n = \gamma/(\gamma-1) \f$
  we can also write using \f$ a = 0 \f$ (constant angular momentum 
  distribution, Kuwabara Eq [8]):
  \f[
     -\frac{1}{R} + \frac{1}{2}\frac{L_k^2}{r^2} 
     + \frac{\gamma}{\gamma-1}K\rho^{\gamma-1} = c
     = -\frac{1}{r_{\min}} + \frac{1}{2}\frac{L_k^2}{r_{\min}^2} 
  \f]
  where the constant on the right hand side  is determined at the 
  zero pressure surface.
  The previous equation is used to compute \c K at the point where
  the density is maximum and then again to obtain the density 
  distribution:
  \f[
     K = \frac{\gamma-1}{\gamma}\left[
          c + \frac{1}{r_{\max}} - \frac{1}{2}\frac{L_k^2}{r_{\max}^2}
         \right]\frac{1}{\rho_{\max}^{\gamma-1}}\,;\qquad
     \rho =  \left[\frac{\gamma-1}{K\gamma}\left(c + \frac{1}{R} 
           - \frac{1}{2}\frac{L_k^2}{r^2}\right)\right]^{1/(\gamma-1)}
  \f]
  Thus, specifying \f$ r_{\min} \f$, \f$ r_{\max} \f$ and \f$ \rho_{\max}\f$
  determine the torus structure completely.
  The torus azimuthal velocity is compute from \f$v_\phi = L_k/r\f$.

  The atmosphere is assumed to be isothermal and it is given by the
  condition of hydrostatic balance:
  \f[
    \rho_a(R) = \eta\rho_{\max}
                \exp\left[\left(\frac{1}{R}-\frac{1}{r_{\min}}\right)
                          \frac{GM}{a^2}\right]
  \f]
  where \f$\eta\f$ is the density contrast between the (maximum) torus
  density and the atmosphere density at \f$(r_{\min},0)\f$.

  The torus surface is defined by the condition \f$p_t > p_a\f$.

  The dimensionless form of the equation employs on the following units:

  - \f$ [\rho] = \rho_{\max}\f$: reference density;
  - \f$ [L]    = R_{*} \f$: star radius;
  - \f$ [v^2] = GM/L \f$: reference (squared) velocity.

  The user defined parameters of this problem, as they appear 
  in pluto.ini, allows for control in the Torus shape and the 
  contrast of physical values: 

  -# <tt>run_config.inputParam[RMIN]</tt>: minimum cylindrical radius for the torus
                                  (inner rim)
  -# <tt>run_config.inputParam[RMAX]</tt>: radius of the Torus where pressure is maximum
  -# <tt>run_config.inputParam[RHO_CUT]</tt>:  minimum density to define the last
                                      contour of the magnetic vec. pot. 
  -# <tt>run_config.inputParam[BETA]</tt>: plasma beta 
  -# <tt>run_config.inputParam[ETA]</tt>   density contrast between atmosphere and Torus
  -# <tt>run_config.inputParam[SCALE_HEIGHT]</tt>:  atmosphere scale height \f$ H=GM/a^2\f$
 
  The magnetic field can be specified to be inside the torus or as a
  large-scale dipole field.
  In the first case, we give the vector potential:
  \f[
    A_\phi = B_0(\rho_t - \rho_{\rm cut})
  \f]
  In the second case, a dipole field is specified in the function
  ::DipoleField().
  The user-defined constant \c USE_DIPOLE can be used to select between the
  two configurations.
  
  \author A. Mignone (andrea.mignone@unito.it) \n
          T. Matsakos
  \date   Mar 23, 2025

  \b References
     - "Global Magnetohydrodynamical Simulations of Accretion Tori", 
       Hawley, J.F., ApJ (2000) 528 462.
     - "The dynamical stability of differentially rotating discs with
        constant specific angular momentum"
        Papaloizou, J.C.B, Pringle, J.E, MNRAS (1984) 208 721
     - "The Acceleration Mechanism of Resistive Magnetohydrodynamic Jets
        Launched from Accretion Disks." 
        Kuwabara, T., Shibata, K., Kudoh, T., Matsumoto, R. ApJ (2005) 621 921   
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if USE_DIPOLE == YES
static void DipoleField(double x1, double x2, double x3, 
                        double *Bx1, double *Bx2, double *A, RunConfig &run_config);
#endif

// static double kappa; // Static variables are not supported in acc routine
#define KAPPA  31   // Last element of run_config->inputParam

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
//
 // * Assign initial conditions in terms of primitive variables
// 
////////////////////////////////////////////////////////////////////////
{
  double rmin, rmax;
  double r, R, A_phi, H;
  double a,c, l_k, l_k2;
  double Vo;
  double prt, pra;
  double rhot, rhoa;
  double kappa;
  
  double gmmr = run_config.gamma/(run_config.gamma - 1.0);

  run_config.small_prs = (1.e-8);

  #if GEOMETRY == CYLINDRICAL
  r  = x1;   // cylindrical radius 
  // z  = x3;   
  R  = sqrt(x1*x1 + x3*x3);           // spherical radius 
  #elif GEOMETRY == SPHERICAL
  r  = x1*sin(x2);   // cylindrical radius 
  // z  = x1*cos(x2);   
  R  = x1;           // spherical radius 
  #endif

  rmax    = run_config.inputParam[RMAX];
  rmin    = run_config.inputParam[RMIN];
  H       = run_config.inputParam[SCALE_HEIGHT];  

  l_k     = sqrt(rmax);
  l_k2    = l_k*l_k;
  Vo      = l_k;

  // ----------------------------------------
  // Solve for kappa and c assuming pra=0
  // ----------------------------------------  
 
  c     =    - 1.0/rmin + 0.5*l_k2/(rmin*rmin);
  kappa = (c + 1.0/rmax - 0.5*l_k2/(rmax*rmax))/gmmr;

  run_config.inputParam[KAPPA] = kappa;

  // ----------------------------
  // Torus density + pressure
  // ----------------------------  

  a    = c + 1.0/R - 0.5*l_k2/(r*r);
  a    = MAX(a, 0.0);
  rhot = pow(a/(gmmr*kappa), 1.0/(run_config.gamma-1));
  prt  = kappa*pow(rhot,run_config.gamma);

  // -------------------------------------
  // Atmospheric density + pressure 
  // (force equilibrium Fg and \nabla P)
  // -------------------------------------

  rhoa = run_config.inputParam[ETA]*exp((1./R - 1./rmin)/H); 
  pra  = rhoa*H;

  if (prt > pra && r > 2.0) {    // Torus 
    v[RHO]   = rhot;
    v[PRS]   = prt;
    v[VX1]   = v[VX2] = v[VX3]= 0.0;
    v[iVPHI] = Vo/r;
    v[TRC]   = 1.0;
  }else{                         // Atmosphere 
    v[RHO]   = rhoa;
    v[PRS]   = pra;
    v[VX1]   = v[VX2] = v[VX3]= 0.0;
    v[iVPHI] = 0.0;
    v[TRC]   = 0.0;
  }

  #if PHYSICS == MHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  A_phi  = 0.0;
  #if USE_DIPOLE == NO
  double beta    = run_config.inputParam[BETA];
  double Bo = sqrt(2.0*kappa/beta);
  double rhocut  = run_config.inputParam[RHO_CUT];  
  if (rhot > rhocut && r > 2.0) A_phi = Bo*(rhot - rhocut);
  #endif
// A_phi  = 0.0;  
// v[BX3] = 0.01;
  
  #if GEOMETRY == SPHERICAL 
  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = A_phi;
  #elif GEOMETRY == CYLINDRICAL
  v[AX1] = 0.0;
  v[AX2] = A_phi;
  v[AX3] = 0.0;
  #endif

  #endif

  #if (USE_DIPOLE == YES) && (BACKGROUND_FIELD == NO)
  #if GEOMETRY == SPHERICAL
  DipoleField (x1,x2,x3,v+BX1, v+BX2, v+AX3);
  #elif GEOMETRY == CYLINDRICAL
  DipoleField (x1,x2,x3,v+BX1, v+BX3, v+AX2);
  #endif
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
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  int i, j, k, nv;
  int ibeg = grid->lbeg[IDIR];

  if (side == X1_BEG){
    #pragma acc parallel loop collapse(3) present(box, grid)
    BOX_LOOP(*box, i,j,k){
      #pragma acc loop seq
      for (nv = 0; nv <  NVAR; nv++) data->Vc(i,j,k,nv) = data->Vc(ibeg,j,k,nv);

      data->Vc(i,j,k,VX1) = MIN(data->Vc(i,j,k,VX1),0.0); 
    }
      // }
//     }else if (box->vpos == X2FACE){
//       #ifdef STAGGERED_MHD
// //       BOX_LOOP(*box,i,j,k) data->Vs[BX2s][k][j][i] = data->Vs[BX2s][k][j][IBEG];
//       #endif
//     }
  }
}

#if (BODY_FORCE & VECTOR)
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  double  R;
  const double x1 = state.x1;

  #if GEOMETRY == CYLINDRICAL
  const double x3 = state.x3;
  R = sqrt(x1*x1 + x3*x3);
  g[IDIR] = -1.0/(R*R*R)*x1;
  g[JDIR] =  0.0; 
  g[KDIR] = -1.0/(R*R*R)*x3;
  #elif GEOMETRY == SPHERICAL
  R = x1;
  g[IDIR] = -1.0/(R*R);
  g[JDIR] =  0.0;
  g[KDIR] =  0.0; 
  #endif

}
#endif

#if (BODY_FORCE & POTENTIAL)
////////////////////////////////////////////////////////////////////////
double BodyForcePotential(const State& state, const RunConfig& run_config)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  return 0.0;
}
#endif

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

  #if USE_DIPOLE == YES
  double A;
  DipoleField (state.x1,state.x2,state.x3, state.Bgf, state.Bgf+1, &A, run_config);
  #else
  #endif
}
#endif

#if USE_DIPOLE == YES
////////////////////////////////////////////////////////////////////////
void DipoleField(double x1, double x2, double x3, 
                 double *Bx1, double *Bx2, double *A, RunConfig &run_config)
/*!
  Assign the polodial component of magnetic fields and vector
  potential for a dipole.
*/
////////////////////////////////////////////////////////////////////////
{
  double r, r2, r3, Bmag;
  double kappa = run_config.inputParam[KAPPA];

  Bmag = sqrt(2.0*kappa/run_config.inputParam[BETA]);

  #if GEOMETRY == CYLINDRICAL
   double R  = x1; double z = x2; r = sqrt(x1*x1 + x2*x2);
   r2 = r*r;
   r3 = r2*r;
   *Bx1 = 3.0*Bmag*z*R/(r2*r3);
   *Bx2 =    -Bmag*(R*R - 3.0*z*z)/(r2*r3);
   *A   = Bmag*R/r3;
  #elif GEOMETRY == SPHERICAL
   r  = x1;
   r2 = r*r;
   r3 = r2*r;
   *Bx1 = 2.0*Bmag*cos(x2)/r3;
   *Bx2 =     Bmag*sin(x2)/r3;
   *A   = Bmag*sin(x2)/r2;
  #endif

}
#endif