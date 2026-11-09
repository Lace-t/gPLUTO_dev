////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Relativistic magnetized blast wave.

  Set the initial condition for the relativistic magnetized blast wave
  problem in 2D or 3D.
  It consists of a highly pressurized region inside a circle (in 2D) or
  a sphere (in 3D) embeddd in a static uniform medium with lower pressure.
  The magnetic field is constant and threads the whole computational domain.
  \f[
     (\rho,\, p) = \left\{\begin{array}{lcl}
        (\rho_{\rm in},\, p_{\rm in})  & \quad\mathrm{or}\quad & r < r_c
        \\[0.5em]
        (\rho_{\rm out},\, p_{\rm out})  & \quad\mathrm{or}\quad & r \ge r_c
     \end{array}\right.
     \,,\qquad
      |\vec{B}| = B_0
   \f]

  \f[
     \vec{B}
  \f]
  In 3D, a linear smoothing is applied in the region \f$ r_c<r<1\f$.
  The input parameters used in this problem are:

  -# <tt>inputParam[PRS_IN]</tt>:  pressure inside the initial circular (2D)
                                    or spherical (3D) region.
  -# <tt>inputParam[PRS_OUT]</tt>: ambient pressure
  -# <tt>inputParam[RHO_OUT]</tt>: ambient density
  -# <tt>inputParam[BMAG]</tt>:    magnetic field intensity
  -# <tt>inputParam[THETA]</tt>:   angle between mag. field and z-axis
                                    (Cartesian only )
  -# <tt>inputParam[PHI]</tt>:     angle between mag. field and xy-plane
                                    (Cartesian only)
  -# <tt>inputParam[RADIUS]</tt>:  radius of the initial over-pressurized region.

  The domain, in each direction, spans over a grid of [-6.0,6.0] in all tests but
  #1 and #2, where a smaller grid of [-0.5,0.5] in all directions is used.
  Consequently, the <tt>inputParam[RADIUS]</tt> is set to 0.08 in the tests
  #1 and #2 and to 0.8 in all the other tests.

  Note that a given choice of parameters can be re-scaled by an arbitrary
  factor \f$\eta\f$  by letting
  \f$ \{\rho,\, p\} \to \eta^2\{\rho,\,p\},\, B \to \eta B \f$.

  Parameters for the different configurations are given in the following table:

  <CENTER>
  Conf.|rho, p (in) | rho, p (out) |  B  | Dim | theta | phi  | Ref
  -----|------------| -------------|-----|-----|-------|------|------------
  #01  | 1,  10^3   | 1.0, 0.01    | 4.0 |  2  | 90.0  | 0.0  | [dZBL2003], Sec 4.2
  #02  | 1,  10^3   | 1.0, 0.01    | 4.0 |  2  | 90.0  | 0.0  | [dZBL2003], Sec 4.2
  #03  | 0.01, 1    | 1.e-4, 3.e-5 | 1.0 |  3  | 0.0   | 90.0 | [Mig_etal2007], Sec 5.7
  #05  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  2  | 90.0  | 0.0  | [BS2011], Sec 4.6
  #06  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  2  | 90.0  | 0.0  | [BS2011], Sec 4.6
  #07  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  2  | 90.0  | 0.0  | [BS2011], Sec 4.6
  #08  | 0.01, 1    | 1.e-4, 5.e-4 | 0.1 |  2  | 90.0  | 0.0  | [Leis_etal2005], Sec B.2.1
  #09  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  2  | 90.0  | 45.0 | [BS2011], Sec 4.6
  #10  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  3  | 0.0   | 90.0 | [MM2022], Sec 4.4
  #11  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  2  | 90.0  | 45.0 | [BS2011], Sec 4.6
  #12  | 0.01, 1    | 1.e-4, 5.e-3 | 1.0 |  3  | 90.0  | 45.0 | [Ber_etal2024], Sec 5.7
  </CENTER>

  The corresponding algorithms are listed here:

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM| RECONST |TIME_STEP|DIVB| Solver
  -----|----|-----------|---|---------|---------|----|-------
   #01 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hllc
   #02 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hll
   #03 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hll
   #04 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hll
   #05 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hlld
   #06 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hlld
   #07 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | GLM| hlld
   #08 |MHD |CARTESIAN  | 2 | WENOZ   | RK3     | CT | hllc
   #09 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | hllc
   #10 |MHD |CARTESIAN  | 2 | LINEAR  | RK2     | CT | gforce
   #11 |MHD |CARTESIAN  | 2 | WENOZ   | RK4 [HO]| CT | gforce
   #12 |MHD |CARTESIAN  | 2 | WENOZ   | RK4 [HO]| CT | gforce
  </CENTER>

  Some of the strongly magnetized configurations can pass this test
  only by taking some precautions (e.g. correcting total energy with
  staggered magnetic field).

  \b References
     - [dZBL2003]      Del Zanna, Bucciantini, Londrillo. A&A (2003) 397
     - [Leis_etal2005] Leismann, et al., A&A (2005) 436, 503
     - [Mig_etal2007]  Mignone et al., ApJS (2007), 170, 228
     - [BS2011]        Beckwith & Stone, ApJS (2011), 193, 6
     - [MM2022]        Mattia & Mignone, MNRAS (2022), 510, 481
     - [Ber_etal2024]  Berta et al., JCP (2024), 499, 112701

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Apr, 15, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double r;
  double *inputParam = run_config.inputParam;
  double rc = inputParam[RADIUS];
  double dc = inputParam[RHO_IN];
  double pc = inputParam[PRS_IN];
  double de = inputParam[RHO_OUT];
  double pe = inputParam[PRS_OUT];
  double theta0 = inputParam[THETA]*CONST_PI/180.0;
  double phi0   = inputParam[PHI]*CONST_PI/180.0;
  double Bx     = inputParam[BMAG]*sin(theta0)*cos(phi0);
  double By     = inputParam[BMAG]*sin(theta0)*sin(phi0);
  double Bz     = inputParam[BMAG]*cos(theta0);

  run_config.gamma = 4.0/3.0;
  #if (GEOMETRY == CARTESIAN)
  r = sqrt(DIM_EXPAND(x1*x1, + x2*x2, + x3*x3));
  #elif (GEOMETRY == CYLINDRICAL)
  r = sqrt(DIM_EXPAND(x1*x1, + 0.0, + x3*x3));
  #elif GEOMETRY == SPHERICAL
  r = x1;
  #endif


  if (r <= rc) {
    v[RHO] = dc;
    v[PRS] = pc;
  #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 3)
  }else if (r > rc && r < 1.0){
    v[RHO] = de*(r - rc)/(1.0 - rc) + dc*(r - 1.0)/(rc - 1.0);
    v[PRS] = pe*(r - rc)/(1.0 - rc) + pc*(r - 1.0)/(rc - 1.0);
  #endif
  }else{
    v[RHO] = de;
    v[PRS] = pe;
  }

  v[VX1] = v[VX2] = v[VX3] = 0.0;

  #if PHYSICS == RMHD
  v[AX1] = v[AX2] = v[AX3] = 0.0;
  #if GEOMETRY == CARTESIAN
  v[BX1]  = Bx;
  v[BX2]  = By;
  v[BX3]  = Bz;

  v[AX1] =   0.0;
  v[AX2] =  v[BX3]*x1;
  v[AX3] = -v[BX2]*x1 + v[BX1]*x2;
  #elif GEOMETRY == CYLINDRICAL
  v[BX1]  = 0.0;
  v[BX2]  = 0.0;
  v[BX3]  = inputParam[BMAG];
  
  v[AX1] = v[AX3] = 0.0;
  v[AX2] = 0.5*v[BX3]*x1;
  #elif GEOMETRY == SPHERICAL
  v[BX1] =  Bx*sin(x2)*cos(x3) + By*sin(x2)*sin(x3) + Bz*cos(x2);
  v[BX2] = -Bx*cos(x2)*cos(x3) + By*cos(x2)*sin(x3) - Bz*sin(x2);
  v[BX3] = -Bx*sin(x3)         + By*cos(x3);
  #endif

// Force Bz to be zero when theta is very close to 90 deg 

  if (fabs(inputParam[THETA] - 90.0) < 1.e-5) v[BX3] = 0.0;
  #endif

  #if GEOMETRY != CARTESIAN
  run_config.small_prs = 1.e-6;
  #endif
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
