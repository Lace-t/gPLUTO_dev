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

  The different configurations are given in the following table

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
  \date    Aug 15, 2024\n
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double r, rc;
  double rho_in, rho_out;
  double prs_in, prs_out;
  double *inputParam = run_config.inputParam;

  run_config.gamma = (4.0/3.0);
  run_config.small_prs = 1.e-6;

  #if (GEOMETRY == CARTESIAN)
  r = sqrt(DIM_EXPAND(x1*x1, + x2*x2, + x3*x3));
  #elif (GEOMETRY == CYLINDRICAL)
  r = sqrt(DIM_EXPAND(x1*x1, + 0.0, + x3*x3));
  #elif GEOMETRY == SPHERICAL
  r = x1;
  #endif

  rc = inputParam[RADIUS];
  rho_in  = inputParam[RHO_IN];
  prs_in  = inputParam[PRS_IN];
  rho_out = inputParam[RHO_OUT];
  prs_out = inputParam[PRS_OUT];

  if (r <= rc) {
    v[RHO] = rho_in;
    v[PRS] = prs_in;
  #if (GEOMETRY == CARTESIAN)
  }else if (r > rc && r < 1.0){
    double e = exp(1.0);
    double ftap = (r - rc)/(1.0 - rc);
    ftap = (exp(ftap) - e)/(1.0 - e);   // Taper function 
    v[RHO] = rho_in*ftap + rho_out*(1.0 - ftap);
    v[PRS] = prs_in*ftap + prs_out*(1.0 - ftap);
  #endif
  }else{
    v[RHO] = rho_out;
    v[PRS] = prs_out;
  }

  v[VX1] = v[VX2] = v[VX3] = 0.0;
  v[AX1] = v[AX2] = v[AX3] = 0.0;
#if PHYSICS == ResRMHD
  v[EX1] = v[EX2] = v[EX3] = 0.0;
  #ifdef CRG
  v[CRG] = 0.0;
  #endif
#endif

  #if GEOMETRY == CARTESIAN
  v[BX1]  = inputParam[BMAG];
  v[BX2]  = 0.0;
  v[BX3]  = 0.0;

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
  v[BX1] =  inputParam[BMAG]*cos(x2);
  v[BX2] = -inputParam[BMAG]*sin(x2);
  v[BX3] = 0.0;
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
  [Not Used]
*/
////////////////////////////////////////////////////////////////////////
{ 
}
