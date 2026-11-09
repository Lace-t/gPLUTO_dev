////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Axisymmetric jet propagation.
 
  \b Description.\n
  The jet problem is set up in axisymmetric cylindrical coordinates 
  \f$ (R,z) \f$.
  We label ambient and jet values with the suffix "a" and "j", respectively.
  The ambient medium is at rest and is characterized by uniform density 
  and pressure, \f$ \rho_a \f$ and \f$ p_a \f$.
  The beam enters from the lower-z boundary from a circular nozzle of 
  radius \f$ R_j \f$ and carries a constant poloidal field \f$ B_z \f$
  and a (radially-varying) toroidal component \f$ B_\phi(R) \f$.
  Flow variables are prescribed as 
  \f[ \left\{\begin{array}{lcl}
     \rho(R)   &=& \rho_j  \\[0.5em]
     v_z(R)    &=& v_j     \\[0.5em]
     B_\phi(R) &=& \DS \left\{\begin{array}{ll}
                        -B_m R/a & \quad{\rm for}\quad R<a \\[0.5em]
                        -B_m a/R & \quad{\rm for}\quad R>a 
                       \end{array}\right. \\[0.5em]
     B_z(R)    &=& B_{z0}\quad\mathrm{(const)}
  \end{array}\right. \f]
  Here \c a is the magnetization radius while \f$ v_R = B_R = 0 \f$.
  These profiles are similar to the ones used by Tesileanu et al. (2008).
  The pressure distribution is found by solving the radial momentum 
  balance between thermal, centrifugal and magnetic forces:
  \f[
    \frac{dp}{dR} = \frac{\rho v_\phi^2}{R} -
                    \frac{1}{2}\left[\frac{1}{R^2}\frac{d(R^2B^2_\phi)}{dR} 
                                    +\frac{dB_z^2}{dR}\right]
  \f]
  Neglecting rotation and assuming \c Bz to be constant the solution of 
  radial momentum balance becomes
  \f[
    p(R) = p_a + B_m^2\left[1-\min\left(\frac{R^2}{a^2},1\right)\right]
  \f]
  and the jet on-axis pressure increases for increasing toroidal field:
  \f[
    p(R=0) \equiv p_j = p_a + B_m^2
  \f]
  where \f$p_a\f$ is the ambient pressure.

  \b Normalization.\n
  Since the MHD equations are scale-invariant we have the freedom to 
  specify a reference length, density and velocity. Here we choose
  - Length: jet radius \f$R_j=1\f$;
  - Density: ambient density \f$\rho_a=1\f$; 
  - Velocity:
    - adiabtic setup (no cooling):
      ambient sound speed: \f$c_a = \Gamma p_a/\rho_a = 1\f$
      (from which it follows \f$p_a = 1/\Gamma\f$).
    - radiative setup (with cooling): 1 Km/s (set by \c unit_velocity).
      In this case, the ambient pressure is computed from the ambient 
      temperature <tt> Ta = 2500 K </tt>.
       
  In this way the number of parameters is reduced to 4:
  -# <tt>inputParam[ETA]</tt>: density contrast 
     \f$ \eta = \rho_j/\rho_a
         \qquad\Longrightarrow\qquad \rho_j = \eta
     \f$
  -# <tt>inputParam[JET_VEL]</tt>:  Jet velocity  \f$ v_j \f$.
     This is also the Mach number of the beam with respect to the ambient
     in the adiabatic setup;
  -# <tt>inpurParam[SIGMA_Z]</tt>: poloidal magnetization;
  -# <tt>inpurParam[SIGMA_PHI]</tt>: toroidal magnetization;
  Magnetization are defined as follows:
  \f[ \left\{\begin{array}{ll}
       \sigma_z    &= \DS \frac{B_z^2}{2p_a}      \\[0.5em]
       \sigma_\phi &= \DS \frac{<B_\phi^2>}{2p_a} 
     \end{array}\right.
   \qquad\Longrightarrow\qquad
     \left\{\begin{array}{ll}
        B_z^2 &= \DS 2\sigma_zp_a       \\[0.5em]
        B_m^2 &= \DS \frac{2\sigma_\phi}{a^2(1/2-2\log a)}p_a
      \end{array}\right.    
  \f]
  Here the average value of \f$B^2_\phi\f$ is simply found as
  \f[
     <B^2_\phi> = \frac{\int_0^1 B^2_\phi R\,dR}{\int_0^1 R\,dR}
                = B_m^2a^2\left(\frac{1}{2} - 2\log a\right)
  \f]
  
  The following MAPLE code can be used to verify the solution:
  \code
   restart;   # Solve radial equilibrium balance equation
   assume (a < 1, a > 0);
   B[phi] := -B[m]*R/a;
   ode    := diff(p(R),R) = -1/(2*R^2)*diff(R^2*B[phi]^2,R);
   dsolve({ode,p(a)=pa});   # only for R < a
  \endcode

  When cooling is enabled, two additional parameters controlling the
  amplitude and frequency of perturbation can be used:
  -# <tt>inputParam[PERT_AMPLITUDE]</tt>: perturbation amplitude (in 
     units of jet velocity);
  -# <tt>inputParam[PERT_PERIOD]</tt>: perturbation period (in years).

  The jet problem has been tested using a variety of configurations, namely:
  <CENTER>
  Conf.| GEOMETRY  | T. STP |RECONSTRUCTION|divB| EOS    | COOLING
  -----|-----------|--------|--------------| ---|--------|--------
   #01 |CYLINDRICAL| RK2    | LINEAR       | CT | IDEAL  |   NO
   #02 |CYLINDRICAL| RK3    | LINEAR       | CT | IDEAL  |   NO 
   #03 |CYLINDRICAL| RK3    | PARABOLIC    | CT | IDEAL  |   NO 
   #04 |CYLINDRICAL| RK2    | LINEAR       | NO | IDEAL  |   NO
   #05 |CYLINDRICAL| RK2    | LINEAR       | CT | IDEAL  |   NO 
   #07 |CYLINDRICAL| RK3    | PARABOLIC    | NO | IDEAL  |   SNEq    
   #10 |CARTESIAN  | RK3    | PARABOLIC    | CT | IDEAL  |   NO 
  </CENTER>

  The following image show the (log) density map at the end of 
  simulation for setup #01.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Nov 26, 2025

  \b References
     - "Simulating radiative astrophysical flows with the PLUTO code:
       a non-equilibrium, multi-species cooling function"
       Tesileanu, Mignone and Massaglia, A&A (2008) 488,429.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp" 

void   GetJetValues (double, double, double, double *, RunConfig &);
double Profile (double, double, RunConfig &);
double Radius (double, double, double);
void   JetPerturbation (double, double, double, double *, const RunConfig &);

#define MAG_RADIUS (USER_DEF_PARAMETERS+1)
#define T_AMBIENT  (USER_DEF_PARAMETERS+2)
#define P_AMBIENT  (USER_DEF_PARAMETERS+3)
#define MU_AMBIENT (USER_DEF_PARAMETERS+4)

// const double a = 0.8;      /* -- Magnetization radius -- */
// static double Ta, pa, mua;   /* -- Ambient pressure -- */
// #if ACC_INIT
// #pragma acc declare create(Ta, pa, mua) copyin(a)
// #endif

////////////////////////////////////////////////////////////////////////
/*!
  Define ambient values and units.
*/
////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
{
  double pa, a = 0.8;
  double *inputParam = run_config.inputParam;
  inputParam[MAG_RADIUS] = a;

  // Set c.g.s units & floor values

  run_config.unit_density  = CONST_amu*1.e3;
  run_config.unit_length   = 2.5e15;
  run_config.unit_velocity = 1.e5;
  #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 3)
  run_config.small_prs     = 1.e-3;
  #endif
 
  // Assign initial conditions

  v[RHO] = 1.0;
  v[VX1] = v[VX2] = v[VX3] = 0.0;
  v[PRS] = pa = 0.6;

  // Need cooling ? 

  #if COOLING == SNEq  // Used in conf #07
  {
    int nv;
    //static int first_call = 1;
    //static double veq[NVAR] = {-1.0};
    //int first_call = 1;

    double veq[NVAR] = {-1.0};
    double nH = 200.0;

    run_config.unit_density   = CONST_amu*200.0;
    run_config.unit_length    = 2.5e15;
    run_config.unit_velocity  = 1.e5;
    run_config.maxCoolingRate = 0.2;
    run_config.minCoolingTemp = 2500;

    double Ta = 2500.0; // Ambient temperature for atomic cooling 
    //if (first_call)
    CompEquil(nH, Ta, veq);
    NIONS_LOOP(nv) v[nv] = veq[nv];
    double mua    = MeanMolecularWeight(v);
    v[PRS] = pa = Ta*v[RHO]/(Kelvin(run_config)*mua);
    inputParam[T_AMBIENT] = Ta;
    inputParam[MU_AMBIENT] = mua;

    //first_call = 0;
  }
  #endif
  inputParam[P_AMBIENT] = pa;
  run_config.small_prs  = 1.e-3;  // pressure floor
  
  // Ambient magnetic field 

  #if PHYSICS == MHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = sqrt(2.0*run_config.inputParam[SIGMA_Z]*pa); 

  #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 3)
  v[AX1] = 0.0;
  v[AX2] = v[BX3]*x2;
  v[AX3] = 0.0;
  #elif GEOMETRY == CYLINDRICAL
  double R = Radius(x1,x2,x3);
  v[AX1] = 0.0;
  v[AX2] = 0.5*R*v[BX3];
  v[AX3] = 0.0;
  #endif

  #endif
  v[TRC] = 0.0;

}

////////////////////////////////////////////////////////////////////////
/*! 
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
  
*/
////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
{ 

}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
{ 

}

////////////////////////////////////////////////////////////////////////
/*! 
  Set the injection boundary condition at the lower z-boundary
  (\c X3-beg must be set to \c userdef in \c pluto.ini).
  For \f$ R \le 1 \f$ we set constant input values (given by the
  ::GetJetValues() function while for $ R > 1 $ the solution has
  equatorial symmetry with respect to the \c z=0 plane.
  To avoid numerical problems with a "top-hat" discontinuous jump,
  we smoothly merge the inlet and reflected value using a profile 
  function ::Profile().
*/
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
{
  int  i, j, k, nv;
  int  kbeg = grid->lbeg[KDIR];

  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];
  #ifdef STAGGERED_MHD
  double *x1s = grid->xL[IDIR];
  #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 3)
  double *x2s = grid->xL[JDIR];
  #endif
  #endif
  RunConfig &run_config = *(data->run_config);
  #if COOLING != NO  
  double time = data->Dts->time; 
  double *inputParam = run_config.inputParam;
  double unit_time   = run_config.unit_length/run_config.unit_velocity;
  #endif

  if (side == X3_BEG){  
	  
    // ------------------------------------------
    // 1. Assign cell-centered values
    // ------------------------------------------

    #pragma acc parallel loop collapse(3) present(box, data)
    BOX_LOOP(*box,i,j,k){
      double vjet[16], vout[16];
      double R = Radius(x1[i], x2[j], x3[k]);

      // 1a. Store jet values 

      GetJetValues (x1[i],x2[j],x3[k], vjet, run_config);

      // 1b. Add perturbation 

      #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 3)
      JetPerturbation (x1[i], x2[j], x3[k], vjet, run_config);
      #endif

      #if COOLING != NO  // Pulsed perturbation 
      double t     = time*unit_time;  // time in seconds 
      double omega = 2.0*CONST_PI/(86400.*365.*inputParam[PERT_PERIOD]);
      vjet[VX3] *= 1.0 + inputParam[PERT_AMPLITUDE]*sin(omega*t);
      #endif

      // 1c. Copy and reflect ambient medium values

      #pragma acc loop seq
      NVAR_LOOP(nv) vout[nv] = data->Vc(i,j,2*kbeg-k-1,nv);
      vout[VX3] *= -1.0;
      #if PHYSICS == MHD
      vout[BX1] *= -1.0;
      vout[BX2] *= -1.0;
      #endif

      // 1d. Blend jet and reflected ambient medium

      #pragma acc loop seq
      NVAR_LOOP(nv){
        double fp = Profile(R, nv, run_config);
        data->Vc(i,j,k,nv) = vout[nv] - (vout[nv] - vjet[nv])*fp;
      }
    }

    // ------------------------------------------
    // 2. Assign staggered magnetic fields
    // ------------------------------------------

    #ifdef STAGGERED_MHD

    // 2a. Assign Bx1 at i-1/2 faces

    #pragma acc parallel loop collapse(3) present(box, data)
    BOXI_LOOP(*box,i,j,k){
      double vjet[16], vout[16];
      double R = Radius(x1s[i],x2[j],x3[k]);
      GetJetValues (x1s[i],x2[j],x3[k], vjet, run_config);

      vout[BX1] = -data->Vs[BX1s](i,j,2*kbeg - k - 1);

      double fp = Profile(R, BX1, run_config);
      data->Vs[BX1s](i,j,k) = vout[BX1] - (vout[BX1] - vjet[BX1])*fp;
    }

    // 2b. Assign Bx2 at j-1/2 faces

    #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 3)
    #pragma acc parallel loop collapse(3) present(box, data)
    BOXJ_LOOP(*box,i,j,k){
      double vjet[16], vout[16];
      double R = Radius(x1[i],x2s[j],x3[k]);
      GetJetValues (x1[i],x2s[j],x3[k], vjet, run_config);

      vout[BX2] = -data->Vs[BX2s](i,j,2*kbeg - k - 1);

      double fp = Profile(R, BX2, run_config);
      data->Vs[BX2s](i,j,k) = vout[BX2] - (vout[BX2] - vjet[BX2])*fp;
    }
    #endif

    #endif
  }
}

////////////////////////////////////////////////////////////////////////
/*! 
  Define jet values as functions of the (cylindrical) radial
  coordinate.
*/
////////////////////////////////////////////////////////////////////////
void GetJetValues (double x1, double x2, double x3, double *vj, RunConfig &run_config)
{
  double  Bz, Bm = 0.0, xi;
  double *inputParam = run_config.inputParam;
  double  R  = Radius(x1, x2, x3);
  double pa  = inputParam[P_AMBIENT];
  double a   = inputParam[MAG_RADIUS];
  
  // static double veq[NVAR];

  xi  = R/a;

  // 1. Define jet values as a function of coordinates and parameters
  //    Vectors are defined in cylindrical coordinates

  vj[RHO] = 1.0 + (inputParam[ETA] - 1.0);
  Bz = sqrt(2.0*inputParam[SIGMA_Z]*pa);
  Bm = sqrt(2.0*inputParam[SIGMA_PHI]*pa/(a*a*(0.5 - 2.0*log(a))));
  
  #if PHYSICS == MHD
  vj[BX1] =  0.0;
  vj[BX2] = -Bm*(xi < 1.0 ? xi:1.0/xi);
  vj[BX3] =  Bz;

  vj[AX1] = 0.0;
  vj[AX2] = 0.5*R*Bz;
  vj[AX3] = 0.0;
  #endif

  vj[VX1] = 0.0;
  vj[VX2] = 0.0;
  vj[VX3] = inputParam[JET_VEL];
  
  vj[PRS] = pa + Bm*Bm*(1.0 - MIN(xi*xi,1.0));
  vj[TRC] = 1.0;

  // 2. Map vector to Cartesian coordinates

  #if GEOMETRY == CARTESIAN
  // VectorTransform  (vj + BX1, x1, x2, x3, CYLINDRICAL);
  // VectorTransform  (vj + VX1, x1, x2, x3, CYLINDRICAL);

  // Transform to Cartesian components
  
  double Bphi = vj[BX2];
  double sin_phi = x1/R;
  double cos_phi = x2/R;

  vj[BX1] =  Bphi*cos_phi;
  vj[BX2] = -Bphi*sin_phi;
  vj[BX3] =  Bz; 
  #endif

  // Need cooling ? 

  #if COOLING != NO
  {
    double Ta  = inputParam[T_AMBIENT];
    double veq[NVAR];
    //static int first_call = 1;
    //double static Tj,muj;

    #if COOLING == H2_COOL
      double Tj,muj,mua;
      Tj = Ta*(1.0 + Bm*Bm/pa)/inputParam[ETA]; // Guess Tj assuming muj = mua 
      //if (first_call){
      CompEquil(200.0*inputParam[ETA], Tj, veq);
      muj = MeanMolecularWeight(veq);
      mua = inputParam[MU_AMBIENT];
      Tj *= muj/mua;   // Improve guess 
      // }
    #else
      //if (first_call) 
      CompEquil(200.0*inputParam[ETA], Ta, veq);
    #endif
    int nv;
    NIONS_LOOP(nv) vj[nv] = veq[nv];
    //first_call = 0;
  }
  #endif  // COOLING != NO 

}

////////////////////////////////////////////////////////////////////////
void JetPerturbation (double x, double y, double z, double *v, 
                      const RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{       
  int m,l;
  double phi   = atan2(x,y);
  double om[8] = {0.03, 0.06, 0.12, 0.25, 0.5, 1.0, 2.0, 3.0};        // temporal frequencies 
  double b[8]  = {0.43, -0.71, 0.1, 0.23, -0.45, 0.17, 0.34, -0.01};  // Random phases
        
  // Add perturbation (as in Mignone etal (2013)
        
  double vr = 0.0;
  double time = run_config.time;
  #pragma acc loop collapse(2) seq reduction(+:vr)
  for (m = 0; m <= 2; m++){
  for (l = 0; l <= 7; l++){
    vr += cos(m*phi + om[l]*time + b[l]);
  }}
  vr *= 0.05*run_config.inputParam[JET_VEL]/24.0;
  
  v[VX2] += vr*cos(phi);
  v[VX1] += vr*sin(phi);

}

////////////////////////////////////////////////////////////////////////
double Profile (double R, double nv, RunConfig &run_config)
/*!
 * Smoothing function used to blend jet & ambient medium values
*/
////////////////////////////////////////////////////////////////////////
{
  double R4 = R*R*R*R, R8 = R4*R4;

  #if (GEOMETRY == CYLINDRICAL)
  if (run_config.inputParam[SIGMA_PHI] > 1.e-9) {    // Start with a smoother profile 
    if (run_config.time < 0.1)  return 1.0/cosh(R4); // with toroidal magnetic fields 
  }
  #else
  #endif

  return 1.0/cosh(R8);
}

////////////////////////////////////////////////////////////////////////
double Radius (double x1, double x2, double x3)
/*!
  Return the cylindrical radius as a function of coordinates
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == CYLINDRICAL
  double R = x1;
  #elif GEOMETRY == CARTESIAN 
  double R = sqrt(x1*x1 + x2*x2);
  #else
    #error Invalid Geometry
  #endif

  return R; 
}


