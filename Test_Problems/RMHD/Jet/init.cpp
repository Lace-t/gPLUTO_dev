////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Relativistic magnetized jet in axisymmetric coordinates.

  The jet problem is set up in axisymmetric cylindrical coordinates
  \f$ (r,z) \f$ as described in section 4.2.3 of 
  Mignone, Ugliano \& Bodo (2009, [MUB09] hereafter).
  We label ambient and jet values with the suffix "a" and "j", respectively.
  The ambient medium is at rest and is characterized by uniform density
  and pressure, \f$ \rho_a \f$ and \f$ p_a \f$.
  The beam enters from the lower-z boundary from a circular nozzle of
  radius 1, carries a constant poloidal field \f$ B_z \f$
  and a (radially-varying) toroidal component \f$ B_\phi(R) \f$.
  Flow variables are prescribed as
  \f[ \left\{\begin{array}{lcl}
     \rho(R)   &=& \rho_j  \\[0.5em]
     v_z(R)    &=& v_j     \\[0.5em]
     B_\phi(R) &=& \displaystyle \left\{\begin{array}{ll}
             -\gamma_j  b_m R/a & \quad{\rm for}\quad R<a \\[0.5em]
             -\gamma_jf b_m a/R & \quad{\rm for}\quad R>a
             \end{array}\right. \\[0.5em]
     B_z(R)    &=& B_{z0}\quad\mathrm{(const)}
  \end{array}\right. 
  \f]
  Here \c a is the magnetization radius while \f$ v_R = B_R = 0 \f$.
  The pressure distribution is found by solving the radial momentum
  balance between thermal, centrifugal and magnetic forces. 
  Neglecting rotation and assuming \c Bz to be constant the solution 
  is given by
  \f[
    p(R) = p_j + b_m^2\left[1-\min\left(\frac{R^2}{a^2},1\right)\right]
  \f]
  Here \f$p_j=p_a\f$ is the jet/ambient pressure at \c R=1 and its
  value depends on the Mach number, see [MUB09].

  The parameters controlling the problem is
  -# <tt>inputParam[MACH]</tt>: the jet Mach number
     \f$ M = v_j/cs \f$ where \f$c_s=\sqrt{\Gamma p_j/(\rho h)}\f$ is 
     the sound speed and \f$\rho h = \rho + \Gamma p_j/(\Gamma-1)\f$.
     This is used to recover \f$p_j\f$;
  -# <tt>inputParam[LORENTZ]</tt>:  the jet Lorentz factor;
  -# <tt>inpurParam[RHOJ]</tt>: the jet density;
  -# <tt>inpurParam[SIGMA_POL]</tt>: magnetization strength for poloidal 
     magnetic field component (see [MUB09], Eq [65]);
  -# <tt>inpurParam[SIGMA_TOR]</tt>: magnetization strength for toroidal
     magnetic field component (see [MUB09], Eq [65]);

  The different configurations are:

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM|   REC     |TIME_STEP  |DIVB|Solver
  -----|----|-----------|---|-----------|-----------|----|------
   #01 |MHD |CYLINDRICAL| 2 | LINEAR    | RK2       | 8W |tvdlf 
   #02 |MHD |CYLINDRICAL| 2 | LINEAR    | RK2       | CT |hllc  
   #03 |MHD |CYLINDRICAL| 2 | LINEAR    | RK3       | CT |hlld  
   #04 |MHD |CYLINDRICAL| 2 | LINEAR    | RK2       | GLM|hlld  
   #05 |MHD |CYLINDRICAL| 2 | LimO3     | RK3       | CT |hlld  
   #06 |MHD |CYLINDRICAL| 2 | LimO3     | RK3       | CT |gforce
   #07 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | 8W |hll   
   #11 |MHD |SPHERICAL  | 2 | LINEAR    | RK2       | 8w |tvdlf 
   #12 |MHD |SPHERICAL  | 2 | LINEAR    | RK2       | CT |hllc  
   #13 |MHD |SPHERICAL  | 2 | LINEAR    | RK3       | CT |hllc  
   #14 |MHD |SPHERICAL  | 2 | LINEAR    | RK2       | GLM|hlld  
  </CENTER>

  - Configurations #01 and #02 employ a purely poloidal field and are similar 
    to section 4.3.3 of Mignone \& Bodo (2006);
  - Configurations #03 and #04 reproduce the configuration in section 4.2.3 of 
    Mignone, Ugliano \& Bodo (2009);
  - Configuration #05, #06 are more strongly magnetized.
  - Configurations #11 - #14 reproduce the previous 4 in spherical coordinates.

  \b References
     - "A five-wave Harten-Lax-van Leer Riemann solver for relativistic
        magnetohydrodynamics", Mignone, Ugliano \& Bodo, MNARS (200) 393, 1141
     - "An HLLC Riemann solver for relativistic flows - II. Magnetohydrodynamics"
        Mignone \& Bodo, MNRAS (2006) 368, 1040   

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 1, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

void   GetJetValues (double, double, double, double *, RunConfig &);
double pj;
#if ACC_INIT
#pragma acc declare create(pj)
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Assign initial conditions
*/
////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
{
  double vjet[32];

  run_config.gamma = 5./3.;
  GetJetValues (x1, x2, x3, vjet, run_config);

// -- set ambient values -- 

  v[RHO] = run_config.inputParam[RHOA];
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = pj;

#if PHYSICS == RMHD
  #if GEOMETRY == CARTESIAN
  v[BX1] = 0.0;
  v[BX2] = vjet[BX2];  
  v[BX3] = 0.0; 

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = vjet[AX3];
  #elif GEOMETRY == CYLINDRICAL
  v[BX1] = 0.0;
  v[BX2] = 0.0;  
  v[BX3] = vjet[BX3]; 

  v[AX1] = 0.0;
  v[AX2] = vjet[AX2];
  v[AX3] = 0.0;
  #elif GEOMETRY == SPHERICAL
  v[BX1] = vjet[BX1];
  v[BX2] = vjet[BX2];
  v[BX3] = 0.0; 

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = vjet[AX3];
  #endif
#endif

  #if NTRACER > 0
  v[TRC] = 0.0;
  #endif
  #ifdef PSI_GLM 
  v[PSI_GLM] = 0.0;
  #endif

}

////////////////////////////////////////////////////////////////////////
/*! 
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
  Assign user-defined inlet jet boundary conditions.
  
  With GLM and SPHERICAL geometry, we solve a reverse 2x2 Riemann 
  problem for psi[n] and psi[L] if we assume that Bx[L] is given and 
  that Bx[n] = Bx[L].
  Here's a MAPLE script to do this:
  
      restart;
      Bx[n] := Bx[L]:   # Normal field should match inner ghost (jet) value

      eq1 := Bx[n]   = (Bx[L]  + Bx[R])/2 - (psi[R] - psi[L])/(2*c):
      eq2 := psi[n] = (psi[L] + psi[R])/2 - (Bx[R] - Bx[L])*c/2:  

      solve({eq1,eq2},{psi[L], psi[n]});
*/
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
{
  int  i, j, k, nv;

  double *x1  = grid->xC[IDIR], *x1L = grid->xL[IDIR];
  double *x2  = grid->xC[JDIR], *x2L = grid->xL[JDIR];
  double *x3  = grid->xC[KDIR];

  double prof;

  RunConfig &run_config = *(data->run_config);

#if GEOMETRY == CARTESIAN
  int  jbeg = grid->lbeg[JDIR];
  if (side == X2_BEG){  // -- X2_BEG boundary -- 
    #pragma acc parallel loop collapse(3) present(box, data)
    BOX_LOOP(*box,i,j,k){
      
      double  vjet[16], vout[NVAR];
      GetJetValues (x1[i], x2[j], x3[k], vjet, run_config);  // Jet Values 
      #pragma acc loop seq
      NVAR_LOOP(nv) vout[nv] = data->Vc(i,2*jbeg-j-1,k,nv); // Ambient 

      vout[VX2] *= -1.0; 
      vout[BX1] *= -1.0;    
      vout[BX3] *= -1.0;
      #ifdef PSI_GLM
      vjet[PSI_GLM]  = data->Vc(i,jbeg,k,PSI_GLM) - data->Vc(i,jbeg,k,BX2);
      vout[PSI_GLM] *= -1.0;
      #endif

      prof = (fabs(x1[i]) <= 1.0); 
      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = vout[nv] - (vout[nv] - vjet[nv])*prof;
    }

    #ifdef STAGGERED_MHD
    double vjet_BX1 = 0.0;
    #pragma acc parallel loop collapse(3) present(box, data)
    BOXI_LOOP(*box,i,j,k){
      double vout_BX1 = -data->Vs[BX1s](i,2*jbeg-j-1,k);
      prof = (fabs(x1L[i]) <= 1.0);
      data->Vs[BX1](i,j,k) = vout_BX1 - (vout_BX1 - vjet_BX1)*prof;
    }
    #endif
  }
#endif

#if GEOMETRY == CYLINDRICAL
  int  kbeg = grid->lbeg[KDIR];
  if (side == X3_BEG){  // -- X2_BEG boundary -- 
    #pragma acc parallel loop collapse(3) present(box, data)
    BOX_LOOP(*box,i,j,k){
      
      double  vjet[16], vout[NVAR];
      GetJetValues (x1[i], x2[j], x3[k], vjet, run_config); // Jet Values 
      #pragma acc loop seq
      NVAR_LOOP(nv) vout[nv] = data->Vc(i,j,2*kbeg-k-1,nv); // Ambient 

      vout[VX3] *= -1.0; 
      vout[BX1] *= -1.0;
      vout[BX2] *= -1.0;
      #ifdef PSI_GLM
      vjet[PSI_GLM]  = data->Vc(i,j,kbeg,PSI_GLM) - data->Vc(i,j,kbeg,BX3);
      vout[PSI_GLM] *= -1.0;
      #endif

      prof = (fabs(x1[i]) <= 1.0); 
      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = vout[nv] - (vout[nv] - vjet[nv])*prof;
    }

    #ifdef STAGGERED_MHD
    double vjet_BX1 = 0.0; 
    #pragma acc parallel loop collapse(3) present(box, data)
    BOXI_LOOP(*box,i,j,k){
      double vout_BX1 = -data->Vs[BX1s](i,j,2*kbeg - k - 1);
      prof = (fabs(x1L[i]) <= 1.0);
      data->Vs[BX1s](i,j,k) = vout_BX1 - (vout_BX1 - vjet_BX1)*prof;
    }
   #endif
  }
#endif


#if GEOMETRY == SPHERICAL
  int  ibeg = grid->lbeg[IDIR];
  if (side == X1_BEG){  // -- X2_BEG boundary -- 
    #pragma acc parallel loop collapse(3) present(box, data) 
    BOX_LOOP(*box,i,j,k){
      double  vjet[16], vout[NVAR];
      GetJetValues (x1[i], x2[j], x3[k], vjet, run_config); // Jet Values 
      #pragma acc loop seq
      NVAR_LOOP(nv) vout[nv] = data->Vc(2*ibeg-i-1,j,k,nv); // Ambient 

      vout[VX1] *= -1.0; 
      vout[BX2] *= -1.0;
      vout[BX3] *= -1.0;
      #ifdef PSI_GLM 
      // vjet[PSI_GLM]  = data->Vc(ibeg,j,k,PSI_GLM) - data->Vc(ibeg,j,k,BX1);
      // vout[PSI_GLM] *= -1.0;
      // Reverse Riemann problem. See comment above. 
      vjet[PSI_GLM]  = vjet[BX1] + data->Vc(ibeg,j,k,PSI_GLM) - data->Vc(ibeg,j,k,BX1);
      vout[PSI_GLM] *= -1.0;
      #endif

      prof = (fabs(x2[j]) <= 1.0); 
      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = vout[nv] - (vout[nv] - vjet[nv])*prof;
    }

    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(box, data) 
    BOXJ_LOOP(*box,i,j,k){
      double  vjet[16];
      GetJetValues (x1[i], x2L[j], x3[k], vjet, run_config); // Jet Values 
      double vout_BX2 = -data->Vs[BX2s](2*ibeg-i-1,j,k);
      prof = (fabs(x2L[j]) <= 1.0);
      data->Vs[BX2s](i,j,k) = vout_BX2 - (vout_BX2 - vjet[BX2])*prof;
    }
    #endif
  }
#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Prescribe user-defined EMF 
*/
////////////////////////////////////////////////////////////////////////
#if CT_EMF_BOUNDARY == YES
void CT_EMF_Boundary(const Data *data, Grid *grid)
{
  int j,k;
  int ibeg = grid->lbeg[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR] + INCLUDE_JDIR;
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR] + INCLUDE_KDIR;
  EMF *emf = data->emf;

  if (grid->lbound[IDIR] != 0){
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
      emf->Ee[KDIR](ibeg,j,k) = 0.0;
    }}
  }
}
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Compute jet profiles.
*/
////////////////////////////////////////////////////////////////////////
void GetJetValues (double x1, double x2, double x3, 
                   double *vj, RunConfig &run_config)
{
  double bm, bphi, lor, a = 0.5;
  double r, x, scrh, sig_z, sig_phi;
  double vz, Bz, Bphi;
  const double gamma       = run_config.gamma;
  const double *inputParam = run_config.inputParam;

  #if GEOMETRY == CARTESIAN
  r = sqrt(DIM_EXPAND(x1*x1,  , + x3*x3));
  #elif GEOMETRY == CYLINDRICAL
  r = x1;
  #else
  r = x1*sin(x2);
  #endif

  lor     = inputParam[LORENTZ];
  sig_z   = inputParam[SIGMA_POL];
  sig_phi = inputParam[SIGMA_TOR];
  if (fabs(r) < 1.e-9) r = 1.e-9;

  x = r/a;
  vj[RHO] = inputParam[RHOJ];

  vz = sqrt(1.0 - 1.0/(lor*lor));  // 3-vel 

  scrh = inputParam[MACH]/vz;
  pj   = vj[RHO]/(gamma*scrh*scrh - gamma/(gamma - 1.0));

  bm  = 4.0*pj*sig_phi;
  bm /= a*a*(1.0 - 4.0*log(a) - 2.0*sig_phi);
  if (bm < 0.0){
    // QuitPLUTO("GetJetValues(): invalid choice of parameters\n",1);
  }
  bm   = sqrt(bm);
 
  scrh    = MIN(x*x, 1.0);
  vj[PRS] = pj + bm*bm*(1.0 - scrh);

  bphi  = bm*(fabs(x) < 1.0 ? x: 1.0/x);

  Bz   = sqrt(sig_z*(bm*bm*a*a + 2.0*pj)); 
  Bphi = lor*bphi;

  #if GEOMETRY == CARTESIAN
  vj[VX1] = 0.0;                        
  vj[VX2] = vz;
  vj[VX3] = 0.0;  

  vj[BX1] = 0.0;     
  vj[BX2] = Bz; 
  vj[BX3] = Bphi;

  vj[AX1] = 0.0;
  vj[AX2] = 0.0;
  vj[AX3] = 0.5*r*Bz;
  #endif

  #if GEOMETRY == CYLINDRICAL
  vj[VX1] = 0.0;                        
  vj[VX2] = 0.0;
  vj[VX3] = vz;  

  vj[BX1] = 0.0;     
  vj[BX2] = Bphi; 
  vj[BX3] = Bz;

  vj[AX1] = 0.0;
  vj[AX2] = 0.5*r*Bz;
  vj[AX3] = 0.0;
  #endif

  #if GEOMETRY == SPHERICAL
  vj[VX1] =  vz*cos(x2);
  vj[VX2] = -vz*sin(x2);
  vj[VX3] = 0.0;

  vj[BX1] =  Bz*cos(x2);
  vj[BX2] = -Bz*sin(x2);
  vj[BX3] =  Bphi;

  vj[AX1] = 0.0;
  vj[AX2] = 0.0;
  vj[AX3] = 0.5*r*sin(x2)*Bz;
  #endif

// -- set tracer value -- 

  #if NTRACER > 0
  vj[TRC] = 1.0;
  #endif
  #ifdef PSI_GLM 
  vj[PSI_GLM] = 0.0;
  #endif

}

