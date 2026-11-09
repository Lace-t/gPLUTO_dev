////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD blast wave.

  The MHD blast wave problem has been specifically designed to show 
  the scheme ability to handle strong shock waves propagating in 
  highly magnetized environments.
  Depending on the strength of the magnetic field, it can become a rather 
  arduous test leading to unphysical densities or pressures if the
  divergence-free condition is not adequately controlled and the numerical 
  scheme does not introduce proper dissipation across curved shock fronts.

  In this example the initial conditions consists of a static medium with 
  uniform density  \f$ \rho = 1 \f$ while pressure and magnetic field 
  are given by 
  \f[
     p = \left\{\begin{array}{ll}
           p_{\rm in}   & \qquad\mathrm{for}\quad  r < r_0 \\[1ex]
           p_{\rm out}  & \qquad\mathrm{otherwise}  
     \end{array}\right.\,;\qquad
     \mathbf{B} = B_0\left(  \sin\theta\cos\phi\hat{\bf x} 
                        + \sin\theta\sin\phi\hat{\bf y}
                        + \cos\theta\hat{\bf z}\right)
  \f]
  The values \f$p_{\rm in},\, p_{\rm out},\, B_0,\, \theta,\, \phi,\, r_0\f$ 
  are control parameters that can be changed from \c pluto.ini using, 
  respectively,

  -# <tt>run_config.inputParam[P_IN]</tt>
  -# <tt>run_config.inputParam[P_OUT]</tt>
  -# <tt>run_config.inputParam[BMAG]</tt>
  -# <tt>run_config.inputParam[THETA]</tt>
  -# <tt>run_config.inputParam[PHI]</tt>
  -# <tt>run_config.inputParam[RADIUS]</tt>

  The over-pressurized region drives a blast wave delimited by an outer
  fast forward shock propagating (nearly) radially while magnetic field 
  lines pile up behind the shock thus building a region of higher magnetic 
  pressure.
  In these regions the shock becomes magnetically dominated and only weakly
  compressive (\f$\delta\rho/\rho\sim 1.2\f$ in both cases).
  The inner structure is delimited by an oval-shaped slow shock adjacent to a
  contact discontinuity and the  two fronts tend to blend together as the
  propagation becomes perpendicular to the field lines.
  The magnetic energy increases behind the fast shock and decreases 
  downstream of the slow shock.
  The resulting explosion becomes highly anisotropic and magnetically confined.

  The available configurations are taken by collecting different setups
  available in literature:
  
  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM|   REC     |TIME_STEP  |DIVB| Solver | Ref
  -----|----|-----------|---|-----------|-----------|----|--------|--------------
   #01 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | tvdlf  |[BS99]
   #02 |MHD |CARTESIAN  | 3 | LINEAR    | RK2       | CT | roe    |[Z04]
   #03 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | hlld   |[BS99]
   #04 |MHD |CARTESIAN  | 3 | LINEAR    | RK2       | CT | roe    |[Z04]
   #05 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | tvdlf  |[MdZ21]
   #06 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | roe    |[MdZ21]
   #07 |MHD |CARTESIAN  | 3 | PARABOLIC | RK3       | CT | hlld   |[GS08],[MT10]
   #08 |MHD |CARTESIAN  | 2 | WENOZ     | RK4 [HO]  | CT | hlld   |[BM24]
   #09 |MHD |CARTESIAN  | 3 | WENOZ     | RK4 [HO]  | CT | gforce |[BM24]
   #10 |MHD |CARTESIAN  | 2 | WENOZ     | RK4 [HO]  | CT | roe    |[BM24]
   #11 |MHD |CARTESIAN  | 3 | WENOZ     | RK4 [HO]  | CT | hlld   |[BM24]
   #12 |MHD |CARTESIAN  | 3 | LINEAR    | RK2       |GLM | hlld   |[GS08],[MT10]
  </CENTER>

  (*)  Setups are in different coordinates and with different orientations 
       of magnetic field using constrained-transport MHD.
  (**) second version in sec. 4.7

  The snapshot below show the solution for configuration #11.

  This setup also works with the \c BACKGROUND_FIELD spliting.
  In this case the initial magnetic field is assigned in the 
  ::BackgroundField() function while the Init() function is used to 
  initialize the deviation to 0.

  \authors A. Mignone (andrea.mignone@unito.it)\n
           V. Berta   (vittoria.berta@unito.it)
  \date    Oct 31, 2025

  \b References: \n
     - [BS99]: "A Staggered Mesh Algorithm using High Order ...",
       Balsara \& Spicer, JCP (1999) 149, 270 (Sec 3.2)
     - [GS08]: "An unsplit Godunov method for ideal MHD via constrained 
        transport in three dimensions", Gardiner \& Stone, JCP (2008) 227, 4123
        (Sec 5.5)
     - [MT10] "A second-order unsplit Godunov scheme for cell-centered MHD:
       The CTU-GLM scheme", Mignone \& Tzeferacos, JCP (2010) 229, 2117
       (Sec 4.7)
     - [Z04]: "A central-constrained transport scheme for ideal 
       magnetohydrodynamics", Ziegler, JCP (2004) 196, 393 (Sec. 4.6)
     - [MdZ21]: "Systematic construction of upwind constrained transport 
       schemes for MHD", Mignone \& Del Zanna, JCP (2021) 424 109748
     - [BM24]: "A 4th-order accurate finite volume method for ..."", 
       Berta, V., Mignone, A., Bugli, M., \& Mattia, G., JCP (2024), 499, 112701.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
 * Assign initial conditions in terms of primitive variables
*/
////////////////////////////////////////////////////////////////////////
{
  double *inputParam = run_config.inputParam;
  double B0    = inputParam[BMAG];
  double theta = inputParam[THETA]*CONST_PI/180.0;
  double phi   = inputParam[PHI]*CONST_PI/180.0;
  double sth   = sin(theta);
  double cth   = cos(theta);
  double sphi  = sin(phi);
  double cphi  = cos(phi);

  // Transform to Cartesian components
  #if GEOMETRY == CARTESIAN
  double x = x1, y = x2, z = x3;
  #elif GEOMETRY == CYLINDRICAL
  double x = x1*cos(x2), y = x1*sin(x2), z = x3;
  #endif

  // Calculate squared distance in Cylindrical Coordinates
  
  #ifdef XOFFSET  // Used in 3D Cylindrical coordinates to test ring average
  double x0 = 0.1;
  double r2 = DIM_EXPAND((x-x0)*(x-x0), + y*y, + z*z);
  #else
  double r2 = DIM_EXPAND(x*x, + y*y, + z*z);
  #endif

  double r = sqrt(r2);

  run_config.gamma = run_config.inputParam[GAMMA_EOS];

  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[PRS] = inputParam[P_OUT];
  if (r <= inputParam[RADIUS]) v[PRS] = inputParam[P_IN];

#if PHYSICS == MHD

  // Reset out of plane component to zero when cth is almost zero

  #if DIMENSIONS == 2
  if (fabs(cth) < 1.e-12) cth = 0.0;
  if (fabs(sth) < 1.e-12) sth = 0.0;
  #endif

  v[BX1] = B0*sth*cphi;
  v[BX2] = B0*sth*sphi;
  v[BX3] = B0*cth;

  #if GEOMETRY == CARTESIAN
  v[AX1] =  0.0;
  v[AX2] =  v[BX3]*x;
  v[AX3] = -v[BX2]*x + v[BX1]*y;
  #elif GEOMETRY == CYLINDRICAL
  v[AX1] = 0.0;
  v[AX2] = 0.5*v[BX3]*x1;
  v[AX3] = 0.0;
  #endif


/*
restart;
Ax := 0:
Ay := B0*r*cos(phi):
Az := 0:

Ar   :=  Ax*cos(phi) + Ay*sin(phi):
Aphi := -Ax*sin(phi) + Ay*cos(phi);

Br   := 1/r*diff(Az,phi)   - diff(Aphi,z);
Bphi := diff(Ar,z)         - diff(Az,r);
Bz   := 1/r*diff(r*Aphi,r) - 1/r*diff(Ar,phi);

v[AX1] =  0.0;
v[AX2] =  v[BX3]*x;
v[AX3] = -v[BX2]*x + v[BX1]*y;
VectorTransform (v+AX1, x1, x2, x3, CARTESIAN);
VectorTransform (v+BX1, x1, x2, x3, CARTESIAN);
*/

  #if BACKGROUND_FIELD == YES
  v[BX1] = v[BX2] = v[BX3] =
  v[AX1] = v[AX2] = v[AX3] = 0.0;
  #endif
#endif   

}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
}

#if BACKGROUND_FIELD == YES
////////////////////////////////////////////////////////////////////////
void BackgroundField (State &state, RunConfig &run_config)
/*!
  Define the component of a static, curl-free background
  magnetic field.
*/
////////////////////////////////////////////////////////////////////////
{
  double sth,cth,sphi,cphi;

  double Bm    = run_config.inputParam[BMAG];
  double theta = run_config.inputParam[THETA]*CONST_PI/180.0;
  double phi   = run_config.inputParam[PHI]*CONST_PI/180.0;

  theta = theta*CONST_PI/180.0;
  phi   =   phi*CONST_PI/180.0;
  sth   = sin(theta);
  cth   = cos(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);

  state.Bgf[IDIR] = Bm*sth*cphi;  
  state.Bgf[JDIR] = Bm*sth*sphi; 
  state.Bgf[KDIR] = Bm*cth;

//VectorTransform (state.Bgf, state.x1, state.x2, state.x3, CARTESIAN);

}
#endif

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
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
