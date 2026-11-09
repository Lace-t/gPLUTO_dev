////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Integrate cooling and reaction source terms.

  Solve the system of ordinary differential equations describing
  optically thin cooling and ionization network (if any).
  We use an adaptive step size integration method and follow the 
  strategy outlined in section 5.1 of Tesileanu et al. (2008) for handling
  stiff cells.

  On output, a time-step estimate for the next time level is computed using
  the relative or absolute variation obtained during the integration of the ODE 
  (from t(n) --> t(n+1))  
  \f[
     \Delta t_c = \min_{ijk}\left[\frac{\Delta t^n M_R}{\epsilon}\right]
     \,\quad\rm{where}\quad
     \epsilon = \max\left(\left|\frac{p^{n+1}}{p^n} - 1\right|,\,
                |X^{n+1}-X^n|\right)
  \f]
  where \f$ M_R \f$ is the maximum cooling rate (defined by the global variable  
  ::max_cooling_rate) and X are the chemical species.
  
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

////////////////////////////////////////////////////////////////////////
/*!
  Integrate cooling and reaction source terms.
  
  \param [in,out]  d     pointer to Data structure
  \param [in]     dt     the time step to be taken
  \param [out]    Dts    pointer to the Time_Step structure
  \param [in]     grid   pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void CoolingSource (Data *data, double dt, timeStep *Dts, Grid *grid)
{
  int nv, i , j ,k;
  uint32_t err_kernel = 0;
  RunConfig &run_config = *(data->run_config);

  double tol       = COOLING_ODE_TOL;
  double gamma_eos = run_config.gamma;
  double small_pressure   = run_config.small_prs; 
  double min_cooling_temp = run_config.minCoolingTemp; 
  double max_cooling_rate = run_config.maxCoolingRate;  
  double err, scrh;
  double mu0, T0, T1, mu1, prs;
  // double v0[NVAR_COOLING], v1[NVAR_COOLING], k1[NVAR_COOLING];
  double dt_cool_tmp;
  
  RBox box;

  // 0. Create box, list of time-dependent variables for this
  //    step. Initialize coefficients.

  RBoxDefine(grid->lbeg, grid->lend, &box);

  // for (nv = NVAR_COOLING; nv--;  ) {
  //   v0[nv] = v1[nv] = k1[nv] = 0.0;
  // }
  //  Dts->dt_cool = dt_cool_tmp = 1.e38;
  dt_cool_tmp = Dts->dt_cool;
  
  // ------------------------------------------------------ 
  // 1. Main loop on interior computational zones
  // ------------------------------------------------------ 
  mynvtxstart_("COOLING_SOURCE",RED);
  
  // pragma acc still does not work. 
  // #pragma acc parallel loop collapse(3) present(data, Dts, grid) copy (dt, dt_cool_tmp) reduction(min:dt_cool_tmp)
  //#pragma acc parallel loop collapse(3) present(data, Dts, grid) copyin (dt_cool_tmp) private(dt) reduction(min:dt_cool_tmp)
  
  #pragma omp target teams distribute parallel for collapse(3) reduction(min:dt_cool_tmp) map(to: dt)
  #pragma acc parallel loop collapse(3) reduction(min:dt_cool_tmp) copyin(dt) 
  BOX_LOOP(box,i,j,k){  
    
    // 1A. Copy solution data to local arrays, compute temperature, 
    //     internal energy, reset flag

    double v0[NVAR_COOLING], v1[NVAR_COOLING], k1[NVAR_COOLING];
    //#pragma acc loop seq
    NVAR_LOOP(nv) v0[nv] = v1[nv] = data->Vc(i,j,k,nv);
    
    prs = v0[PRS];
    mu0 = MeanMolecularWeight(v0);
    T0  = v0[PRS]/v0[RHO]*Kelvin(run_config)*mu0;

    data->err(i,j,k) = 0;

    #if EOS == IDEAL
    v0[RHOE] = v1[RHOE] = prs/(gamma_eos-1.0);
    #else
      //v1[RHOE] = v0[RHOE] = InternalEnergy(v0, T0); //insert this in future
    #endif

    // insert inside ERROR function 
    /*
    if (T0 <= 0.0){
      printLog ("! CoolingSource(): negative initial temperature\n");
      printLog (" %12.6e  %12.6e\n",v0[RHOE], v0[RHO]);
      printLog (" at: %f %f\n",gCooling_x1, gCooling_x2);
      QUIT_PLUTO(1);
    }
    */
   
    // 1B. Take one step on the ODEs from t1 = 0...dt

    Radiat(v0, k1, run_config);

    uint32_t local_err = 0;
    uint32_t errRKF23;
    #if COOLING_TIME_STEPPING == RKF23
    errRKF23  = SolveODE_RKF23(v0, k1, v1, dt, tol, run_config);
    #elif COOLING_TIME_STEPPING == ROS34
    local_err = SolveODE_Ros34(v0, k1, v1, dt, tol, run_config);
    #elif COOLING_TIME_STEPPING == SELECTIVE
    int stiff = 0;

    stiff     = StiffnessDetection(v0, k1, dt, run_config, 0);
    // v1[TRC]   = stiff;
    // NTRC_LOOP(nv) data->Vc(i,j,k,nv) = v1[nv];

    if (stiff == 0){
      errRKF23 = SolveODE_RKF23(v0, k1, v1, dt, tol, run_config);
      if (errRKF23 != 0){
         //stiff = StiffnessDetection(v0, k1, dt, run_config, 1);
         stiff = 1;
      }
    }

    if (stiff){
      local_err = SolveODE_Ros34(v0, k1, v1, dt, tol, run_config);
    }

    #else 
      #error Invalid COOLING_TIME_STEPPING
    #endif

    // 1C. Update error flag after conversion

    ATOMIC_ERR_UPDATE(err_kernel, local_err);
    KERNEL_ERR_UPDATE(data->err(i,j,k), local_err);
    if (local_err & ERR_CONS2PRIM_FAIL) {
      ATOMIC_ERR_UPDATE(err_kernel,  ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_QUIT);
    }

    // 1D. As asafety measure, constrain the ions to
    //     lie in [0,1] and pressure to be positive.
      
    //#pragma acc loop seq
    NIONS_LOOP(nv){
      v1[nv] = MAX(v1[nv], 0.0);
      v1[nv] = MIN(v1[nv], 1.0);
    }

    #if COOLING == H2_COOL // not implemented 
    v1[X_H2] = MIN(v1[X_H2], 0.5);
    #endif
    
    mu1 = MeanMolecularWeight(v1);

    #if EOS == IDEAL
    prs = v1[RHOE]*(gamma_eos - 1.0);
    T1  = prs/v1[RHO]*Kelvin(run_config)*mu1;
    #elif EOS == PVTE_LAW
    status = GetEV_Temperature(v1[RHOE], v1, &T1);
    prs    = v1[RHO]*T1/(Kelvin(run_config)*mu1);
    #endif

    if (prs < 0.0) prs = small_pressure;  // Is this necessary ? 

    // 1E. If gas cooled below lower threshold, redefine
    //     pressure so that T = min_cooling_temp.

    if (T1 < min_cooling_temp && T0 > min_cooling_temp){
      prs = min_cooling_temp*v1[RHO]/(Kelvin(run_config)*mu1);
    }

    // 1F. Suggest next time step based on fractional variaton
    //     of pressure and ions

    err = fabs(prs/v0[PRS] - 1.0);
    
    #if COOLING == MINEq //not implemented
    for (nv = NFLX; nv < NFLX + NIONS - Fe_IONS; nv++) {
      err = MAX(err, fabs(v1[nv] - v0[nv]));
    }
    #else
    //#pragma acc loop seq
    NIONS_LOOP(nv) err = MAX(err, fabs(v1[nv] - v0[nv]));
    #endif 

    scrh = dt*max_cooling_rate/(err);  // Prevent division by 0  

    dt_cool_tmp = MIN(dt_cool_tmp, scrh); //this should be added for reduction for GPU
    //Dts->dt_cool  = MIN(Dts->dt_cool, scrh);

    // 1G. Update main solution array

    data->Vc(i,j,k,PRS) = prs;
    NIONS_LOOP(nv) data->Vc(i,j,k,nv) = v1[nv];
  } // End BOX_LOOP() 

  mynvtxstop_();

  // ------------------------------------------------------
  // 2. Update cooling time step, convert to conservative
  // ------------------------------------------------------

  Dts->dt_cool = dt_cool_tmp; 
  // #pragma acc update device(Dts->dt_cool)

  PrimToCons3D(data, 0, grid, &box);

  // ------------------------------------------------------
  // 3. Error Check condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("CoolingSource()", data, err_kernel, &box, grid);
  #endif

  // 3. Quit check: by default, quit if ERR_QUIT is raised
  //    unless FALLBACK is activated and we're currently 
  //    running using the DEFAULT method.

  #if WARNING_MESSAGES >= 0 
  int quit_condition = (err_kernel & ERR_QUIT);

  #if FALLBACK_MODE != NO
  quit_condition *= (data->integration_mode != DEFAULT);
  Fallback_Trigger(data->err, data->flag, box);
  #endif

  if (quit_condition) QuitPLUTO("from CoolingSource()",err_kernel);
  #endif
}
 

////////////////////////////////////////////////////////////////////////
/* StiffnessDetection:
   Detects whether the local cooling + chemistry ODE system is stiff
   in the current cell. The decision is based on:
     (i)  the fractional variation of energy and ion fractions over
          the hydrodynamic timestep, and
     (ii) the deviation from ionization equilibrium (NEI).
   Returns 0 if an explicit integrator (e.g. RK3) is safe,
   returns 1 if a stiff (semi-implicit) solver is required.           */
////////////////////////////////////////////////////////////////////////
int StiffnessDetection(double *v, double *rhs, double dt_hydro, 
                       RunConfig &run_conf, int debug)
{
  #if COOLING == SNEq
  double veq[NVAR_COOLING];
  double rho = v[RHO];
  double prs, N, T, st, N_H_rho;
  double unit_time = run_conf.unit_length/run_conf.unit_velocity;
  double mu  = MeanMolecularWeight(v); 
  #endif
  
  // 1. Compute Temperature, number density and equilibrium values

  #if EOS == IDEAL
  if (v[RHOE] < 0.0) v[RHOE] = run_conf.small_prs/(run_conf.gamma-1.0);
  #endif

  #if COOLING == SNEq
  prs = v[RHOE]*(run_conf.gamma-1.0);
  T   = prs/rho*Kelvin(run_conf)*mu; 
  st  = sqrt(T);
  N_H_rho    = (run_conf.unit_density/CONST_amu)*(H_MASS_FRAC/CONST_AH);

  double fn  = v[X_HI];
  N   = N_H_rho*rho;  // Number density of hydrogen N_H = N(HI) + N(HII)

  CompEquil (N, T, veq);

  // Measures the distance from equilibrium values. 
  // This is intended only for cooling with chemical species
  // NEI = Non-Equilibrium Ionization.
  double delta_NEI = 0.0; 
  delta_NEI = fabs(v[X_HI] - veq[X_HI]);
  #endif

  // 2. Compute eigenvalues and stiffness parameter (intensive)

  #if 0
  double J[2][2];  
  Numerical_Jacobian (v, J, run_conf);
  double a = 1.0;
  double b = -(J[0][0] + J[1][1]);
  double c =   J[0][0]*J[1][1] - J[0][1]*J[1][0];
  double delta2 = b*b - 4.0*a*c;
  double lambda_max = 0.0, lambda_min = 1.0;
  if (delta2 > 0.0){
    double xp = (-b + sqrt(delta2))/(2*a);
    double xm = (-b - sqrt(delta2))/(2*a);
    xp = MIN(xp,0.0);
    xm = MIN(xm,0.0);
    lambda_max = MAX(fabs(xp), fabs(xm));
    lambda_min = MIN(fabs(xp), fabs(xm));
  }
  double rstiff = lambda_max/lambda_min;
  double hlambda =  dt_hydro*lambda_max;
  #endif

  // 4. Use Fractional-variation criterion to determine "stiffness"

  double dt_cool   = v[RHOE]/(fabs(rhs[RHOE]) + 1.e-15*v[RHOE]);
  double eps_chem = 0.0, eps_cool, eps;
  #if COOLING == SNEq
  double n_el = N*(1.0 -  v[X_HI] + FRAC_Z);  // electron number density, in cm^{-3} 
  double cr = 2.6e-11/st;
  double ci = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6);
  double dt_ion  = 1.0/(unit_time*n_el*ci);  // ionization    time scale
  double dt_rec  = 1.0/(unit_time*n_el*cr);  // recombinatoin time scale
  double dt_chem = MIN(dt_ion, dt_rec);

  eps_chem = dt_hydro/dt_chem;
  #endif

  eps_cool = dt_hydro/dt_cool;
  eps      = MAX(eps_cool, eps_chem);

  // 5. Print some diagnostic (on demand)
  #if 0 
  if (debug){
    ShowState(v,1);
    printf ("  T = %8.3e; fn = %f; eps_cool = %8.3e; eps_chem = %8.3e, dNEI = %f;\n", 
              T, fn, eps_cool, eps_chem, delta_NEI);
    //printf ("  r = %f, h*lambda = %f\n", rstiff, hlambda);
  
    // if (rstiff > 99000){

    //   int n = 2;
    //   printLog ("----------------------------------------------------------------\n");
    //   for (int k1 = 0; k1 < n; k1++){
    //     for (int k2 = 0; k2 < n; k2++){
    //       printLog ("%12.3e   ", fabs(J[k1][k2]) > 1.e-12 ? J[k1][k2]:0.0);
    //     }
    //     printLog ("\n");
    //   }
    //   printLog ("----------------------------------------------------------------\n");
    // }
  }
  #endif
//  if (eps < 0.1 && delta_NEI < 0.05) return 0;
  if (eps < 0.1) return 0;


  return 1; 
}
