////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Integration methods for cooling netowrk ODEs.

  Solve the system of ordinary differential equations describing
  optically thin cooling and ionization network (if any).

  \b References
     - "Simulating radiative astrophysical flows with the PLUTO code:
        a non-equilibrium, multi-species cooling function" \n
       Tesileanu, Mignone and Massaglia, A&A (2008) 488, 429

  \authors A. Mignone (andrea.mignone@unito.it)\n
           S. Truzzi \n
           O. Tesileanu
  \date    Jan 20, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if 0 // not implemented yet
  #if COOLING == MINEq
  #include "cooling_defs.hpp"
  #endif
#endif // not implemented yet end

//void Numerical_Jacobian (double *, double J[][(NIONS+1)], RunConfig &);  

////////////////////////////////////////////////////////////////////////
/*!
  Explicit 2nd order integrator with 1st-order embedded solution.
  Attempt to integrate the cooling network with a single step
  if the error is within the given tolerance.
  The error is computed from the difference between 1st and 2nd order
  solutions.
  
  \param [in]    v0    initial condition array
  \param [in]    k1    right hand side, R(v0)
  \param [out]   v2nd  final solution array
  \param [in]    dt    the time step we want to achieve
  \param [in]    tol   tolerance used to compute the error
  
  \return The (scaled) error between 1st and 2nd-order solutions
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
uint32_t SolveODE_RKF12 (double *v0, double *k1, double *v2nd, 
                        double dt, double tol, RunConfig& run_conf)
{ 
  int  nv;
  double err;
  double v1[NVAR_COOLING], v1st[NVAR_COOLING], vscal[NVAR_COOLING];
  double k2[NVAR_COOLING];
  
  // --------------------------------------------
  // 0. Copy ALL variables so that density is 
  //    defined when calling Radiat.
  // --------------------------------------------

  for (nv = 0; nv < NVAR_COOLING; nv++) v1[nv] = v0[nv];

  // 2. Get K2 
  
  NVAR_COOLING_LOOP(nv) v1[nv] = v0[nv] + 0.5*dt*k1[nv];
  Radiat(v1, k2, run_conf);

  // 3. Obtain 1st and 2nd order solution -- */

  NVAR_COOLING_LOOP(nv) v2nd[nv] = v0[nv] + dt*k2[nv];
  NVAR_COOLING_LOOP(nv) v1st[nv] = v0[nv] + dt*k1[nv];

  // 4. Compute error

  vscal[RHOE] = fabs(v0[RHOE]) + dt*fabs(k1[RHOE]);
 
  NIONS_LOOP(nv) vscal[nv] = 1.0;  // fabs(v0[nv]) + dt*fabs(k1[nv]) + 1.e-6;     

  err = 0.0;

  
  return err;
}

////////////////////////////////////////////////////////////////////////
/*!
  Explicit 3rd order integrator with 2nd-order embedded solution.
  Attempt to integrate the cooling network with a single step
  if the error is within the given tolerance.
  \f[ \begin{array}{l}
  y^{\rm 3rd} = y_0 + \Delta t (k_1 + k_2 + 4k_3)/6  \\
  k_1 = RHS(y(0))                             \\
  k_2 = RHS(y(0) + \Delta t k_1)                     \\
  k_3 = RHS(y(0) + \Delta t(k_1 + k_2)/4)         \\  
  y^{\rm 2nd} = y(0) + \Delta t(k_1 + k_2)/2
  \end{array}
  \f]
  
  The error is computed from the difference between 2nd and 3rd order
  solutions.
  
  This method is used as low-order integrator when the rhs is not stiff. 
  
  \param [in]    v0    initial condition array
  \param [in]    k1    right hand side, R(v0)
  \param [out]   v3rd  final solution array
  \param [in]    dt    the time step we want to achieve
  \param [in]    tol   tolerance used to compute the error
  
  \return The (scaled) error between 2nd and 3rd-order solutions
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
uint32_t SolveODE_RKF23 (double *v0, double *k1, double *v3rd, 
                        double dt, double tol, RunConfig& run_conf)
{ 
  int  nv;
  double err= 0.0, scrh, dt_4, dt_6;
  double v1[NVAR_COOLING], v2nd[NVAR_COOLING], vscal[NVAR_COOLING];
  double k2[NVAR_COOLING], k3[NVAR_COOLING];

  dt_4 = 0.25*dt;
  dt_6 = dt/6.0;

  // --------------------------------------------
  // 0. Copy ALL variables so that density is 
  //    defined when calling Radiat.
  // --------------------------------------------

  for (nv = 0; nv < NVAR_COOLING; nv++) v1[nv] = v0[nv];

  // 1a. Obtain k2

  NVAR_COOLING_LOOP(nv) v1[nv] = v0[nv] + dt*k1[nv];
  Radiat(v1, k2, run_conf);

  // 1b. Obtiain k3

  NVAR_COOLING_LOOP(nv) v1[nv] = v0[nv] + dt_4*(k1[nv] + k2[nv]);
  Radiat(v1, k3, run_conf);

  // 1c. Compute 3rd- and 2nd-order solutions

  NVAR_COOLING_LOOP(nv) {
    v3rd[nv] = v0[nv] + dt_6*(k1[nv] + k2[nv] + 4.0*k3[nv]);
    v2nd[nv] = v0[nv] + 0.5*dt*(k1[nv] + k2[nv]);
  }

  // 2. Compute error

  vscal[RHOE] = fabs(v0[RHOE]) + dt*fabs(k1[RHOE]);
  NIONS_LOOP(nv) vscal[nv] = 1.0;  // fabs(v0[nv]) + dt*fabs(k1[nv]) + 1.e-6;

  err = 0.0;
  NVAR_COOLING_LOOP(nv){
    scrh = fabs(v3rd[nv] - v2nd[nv])/fabs(vscal[nv]);
    err  = MAX(err, scrh);
  }
  err /= tol;

  return err > 1;
}


#ifndef ROS34_MAXSTEPS
  #define ROS34_MAXSTEPS 40
#endif
#define DEBUG  NO

//////////////////////////////////////////////////////////////////////// 
/* Integrates the local cooling + chemistry ODE system over a single
   hydrodynamic timestep using a 4th-order Rosenbrock method with an
   embedded 3rd-order error estimator (ROS3/4).
   The solver uses adaptive sub-stepping and a numerical Jacobian to
   handle stiff regimes (rapid cooling and/or NEI effects).
   On success, returns 0 and stores the 4th-order solution in v4th.      */
////////////////////////////////////////////////////////////////////////
uint32_t SolveODE_Ros34 (double *v0, double *k1, double *v4th, 
                         double dt, double tol, RunConfig& run_conf)
{
  int      i, j, n, nv, ksub, kfail, count;
  int      done = 0;
  double   err, scrh, t, tend;
  double   v1[NVAR_COOLING], vscal[NVAR_COOLING];
  double   k2[NVAR_COOLING], k3[NVAR_COOLING];
  double   tsub[4096];
  double dt_grow, dt_shrink;
  double   a[NIONS+1][NIONS+1];
  double   J[NIONS+1][NIONS+1];
  double   vbeg[NVAR_COOLING];
  const double GAM = (1.0/2.0);
  const double A21 = 2.0;
  const double A31 = (48.0/25.0);
  const double A32 = (6.0/25.0);
  const double C21 = -8.0;
  const double C31 = (372.0/25.0);
  const double C32 = (12.0/5.0);
  const double C41 = (-112.0/125.0);
  const double C42 = (-54.0/125.0);
  const double C43 = (-2.0/5.0);
  const double B1  = (19.0/9.0);
  const double B2  = (1.0/2.0);
  const double B3  = (25.0/108.0);
  const double B4  = (125.0/108.0);
  const double E1  = (17.0/54.0);
  const double E2  = (7.0/36.0);
  const double E3  = 0.0;
  const double E4  = (125.0/108.0);

#if DEBUG == YES
dt = 2.e-5;
v0[RHO] = 2.084967279012e+00;
v0[VX1] = 7.068272651760e-05;
v0[VX2] = 0.000000000000e+00;
v0[VX3] = 5.816701934575e+01;
v0[PRS] = 1.975306290846e+03;
v0[BX1] = 0.000000000000e+00;
v0[BX2] = -1.615612985541e-01;
v0[BX3] = 0.000000000000e+00;
v0[RHOE] = v0[PRS]/(run_conf.gamma-1.0);
v0[X_HI] = 9.112619309036e-01;
Radiat (v0, k1, run_conf);
StiffnessDetection(v0, k1, dt, run_conf, 1);

double k10[NVAR_COOLING];
for (nv = 0; nv < NVAR_COOLING; nv++) k10[nv] = k1[nv];
double dt0 = dt;
#endif


  // 0. Copy ALL variables so that density is
  //    defined before calling Radiat.  

  for (nv = 0; nv < NVAR_COOLING; nv++) v1[nv] = vbeg[nv] = v0[nv];
  
  n = NIONS + 1;
 
  double g1[NIONS+1], g2[NIONS+1], g3[NIONS+1], g4[NIONS+1];
  int indx[NIONS+1];

  for (i = 0; i < n - 1; i++) vscal[i + IONS_BEG] = 1.0;

  // -- do not take error on Fe -- 

  #if COOLING == MINEq && Fe_IONS > 0
    vscal[FeI] = vscal[FeII] = vscal[FeIII] = 1.e18;
  #endif

  vscal[RHOE] = fabs(v0[RHOE]);

  Numerical_Jacobian (v0, J,run_conf); 

  t    = 0.0;
  tend = dt;
  ksub = kfail = 0;

  // ------------------------------------------------------
  // 1. Take one step (with adaptive step size control) 
  //    using Rosenbrock 4th order (with embedded 3rd order)
  // ------------------------------------------------------

  for (count=0; count < ROS34_MAXSTEPS; count++){

    // 1a. Compute (I - hJ)

    for (i = 0; i < n; i++) {   
      for (j = 0; j < n; j++) a[i][j] = -J[i][j];
      a[i][i] += 1.0/(GAM*dt);
    }

    LUDecompose (&a[0][0], n, indx, &scrh);  // LU decomposition of the matrix
   
    // 1b. Set right hand side and solve for g1 

    NVAR_COOLING_LOOP(i) g1[i-IONS_BEG] = k1[i];
    LUBackSubst (&a[0][0], n, indx, g1);    // Solve for g1
    NVAR_COOLING_LOOP(i) v1[i] = vbeg[i] + A21*g1[i - IONS_BEG];

    Radiat (v1, k2, run_conf);      

    // 1c. Set right hand side and solve for g2

    NVAR_COOLING_LOOP(i) g2[i-IONS_BEG] = k2[i] + C21*g1[i-IONS_BEG]/dt;
    LUBackSubst (&a[0][0], n, indx, g2);   // Solve for g2 
    NVAR_COOLING_LOOP(i){
      v1[i] = vbeg[i] + A31*g1[i - IONS_BEG] + A32*g2[i - IONS_BEG];
    }

    Radiat (v1, k3, run_conf);
    
    // 1d. Set right hand side and solve for g3 

    NVAR_COOLING_LOOP(i){
      g3[i-IONS_BEG] = k3[i] + (C31*g1[i-IONS_BEG] + C32*g2[i-IONS_BEG])/dt;
    }
    LUBackSubst (&a[0][0], n, indx, g3);    // Solve for g3

    // 1e. Set right hand side and solve for g4 

    NVAR_COOLING_LOOP(i){
      g4[i-IONS_BEG] = k3[i] + (  C41*g1[i-IONS_BEG] + C42*g2[i-IONS_BEG] 
                                + C43*g3[i-IONS_BEG])/dt;
    }
    LUBackSubst (&a[0][0], n, indx, g4);    // Solve for g4

    // 1f.  4th order solution & error estimation 

    err = 0.0;
    NVAR_COOLING_LOOP(i){   // Get 4th-order estimate of y and error estimate.
      int ib = i - IONS_BEG;
      v4th[i] = vbeg[i] + B1*g1[ib] + B2*g2[ib] + B3*g3[ib] + B4*g4[ib];
      scrh    = E1*g1[ib] + E2*g2[ib] + E3*g3[ib] + E4*g4[ib];
      err     = MAX(err, fabs(scrh)/vscal[i]);
    }

    err /= tol;   // Scale relative to required tolerance.

    // 1g. Check error

    if (err < 1.0) {  // Step succeeded. Compute size of next step and return.

      ksub++;      
      err = MAX(err, 1.e-18);

      t          += dt;
      tsub[ksub]  = t;

      // Provide an estimate for next dt

      dt_grow = 0.9*dt*pow(err, -0.25);
      dt_grow = MIN(dt_grow, 5.0*dt);
     // dt_grow = (err > 0.1296 ? 0.9*dt*pow(err,-0.25) : 1.5*dt); [NR Original]

      dt      = dt_grow;

      // Exit loop if final time has been reached

      if (fabs(t/tend - 1.0) < 1.e-9) done = 1;

      // Adjust dt not to exceed tend -- */

      if (dt > (tend - t)) dt = tend - t;

      // Initialize solution vector before continuing
   
      if (!done){
        NVAR_COOLING_LOOP(nv) vbeg[nv] = v4th[nv];
        Radiat (vbeg, k1, run_conf);
        Numerical_Jacobian (vbeg, J, run_conf); 
      }

    }else{   // Truncation error too large. Shrink dt and redo time step

      kfail++;
      dt_shrink = 0.9*dt*pow(err, -1.0/3.0);
      dt        = MAX(dt_shrink, 0.05*dt);
    }

    #if DEBUG == YES
    printf ("  - count = %d; err = %f; dt = %8.3e; t/tend = %f\n",
             count,err, dt, t/tend);
    #endif

    if (done) break;
  }

  #if DEBUG == YES
  if (count > MAX(ROS34_MAXSTEPS-5,1)){
    printf ("Ros34: count = %d, dt0 = %12.6e\n",count, dt0);

    printf ("k1 = %12.6e  %12.6e\n",k10[RHOE], k10[X_HI]);
    double rho = v0[RHO];
    double mu  = MeanMolecularWeight(v0); 
    double prs = v0[RHOE]*(run_conf.gamma-1.0);
    double T   = prs/rho*Kelvin(run_conf)*mu; 
    printf ("T = %8.3e; X_HI = %18.12e\n",T, v0[X_HI]);
    StiffnessDetection(v0, k10, dt0, run_conf, 1);
    QuitPLUTO ("Stop in Ros34 ",1);
  }

  printf ("Exit ok. count = %d\n",count);
  QuitPLUTO("",0); 
  #endif

  if (count == ROS34_MAXSTEPS) return ERR_TOO_MANY_ITERATIONS;
  else return 0;
}                                                                        
////////////////////////////////////////////////////////////////////////
/*!
  Compute Jacobian matrix numerically.
  
  Note: This function is time consuming.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void Numerical_Jacobian (double *v, double J[][(NIONS+1)], RunConfig& run_conf)
{
  int  n, nv, k, l;
  double eps;
  double vp[NVAR_COOLING], vm[NVAR_COOLING];
  double Rp[NVAR_COOLING], Rm[NVAR_COOLING];

  eps = 1.e-5;
  n   = NIONS + 1;

  // 1. Take partial derivatives with respect to pressure using central 
  //    differencincing.

  for (nv = 0; nv < NVAR_COOLING; nv++){
    vp[nv] = vm[nv] = v[nv];
  }
  vp[RHOE] *= 1.0 + eps;
  vm[RHOE] *= 1.0 - eps;
 
  Radiat (vp, Rp, run_conf);
  Radiat (vm, Rm, run_conf);

  for (k = 0; k < n; k++){
    int k1 = IONS_BEG + k;
    J[k][n - 1] = (Rp[k1] - Rm[k1])/(2.0*v[RHOE]*eps); //  drhs(*) / drhoe
  }

  // 2. Take partial derivatives with respect to ions

  for (l = 0; l < n - 1; l++){

    int l1 = IONS_BEG + l;

    for (nv = 0; nv < NVAR; nv++) vp[nv] = vm[nv] = v[nv];

    vp[l1] = v[l1] + eps;
    vm[l1] = v[l1] - eps;
                      
    vp[l1] = MIN(vp[l1], 1.0);
    vm[l1] = MAX(vm[l1], 0.0);
                         
    Radiat (vp , Rp,run_conf);
    Radiat (vm , Rm,run_conf);

    for (k = 0; k < n; k++){
      int k1 = IONS_BEG + k;
      J[k][l] = (Rp[k1] - Rm[k1])/(vp[l1] - vm[l1]); 
    }

  }
}

#undef DEBUG