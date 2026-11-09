// /////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Radiation implicit step
  
  Main routines used for the implicit integration of the radiation-matter interaction
  terms. Four different implicit methods are here implemented:
  	
  - RADIATION_FIXEDPOINT_RAD: Solves the implicit step following Takahashi & Oshuga
		(2013), by iterating the radiation fields.
		
  - RADIATION_FIXEDPOINT_GAS: fixed-point algorithm based on iterations of the MHD
		fields.

  - RADIATION_NEWTON_GAS: Solves the implicit step by means of the Newton's method,
    performing iterations of the matter fields (p_g,u^i), i.e., the matter
    pressure and spatial components of its 4-velocity. The components of the Jacobian
    are computed numerically by carrying small variations of these fields.

  - RADIATION_NEWTON_RAD: Same as RADIATION_NEWTON_RAD, but performing iterations of
		the radiation fields.

	In all methods the convergence criterion considers both radiation and MHD fields
	unless RADIATION_FULL_CONVERGENCE == NO, in which case only the variation of the
	iterated fields is controlled. 
  
  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 09, 2026
  
	\b References
	 -  Melon Fuksman, J. D., and Mignone, A. 2019, "A Radiative
			Transfer Module for Relativistic Magnetohydrodynamics in the
			PLUTO Code" ApJS, 242, 20.
   -  Melon Fuksman, J. D., Klahr, H., Flock, M. and Mignone, A. 2021,
      "A Two-moment Radiation Hydrodynamics Scheme Applicable to Simulations
      of Planet Formation in Circumstellar Disks" ApJ, 906, 78.
   -  Nonrelativistic implicit step: Melon Fuksman, J. D., & Klahr, H. 2022, ApJ, 936, 16. Appendix B.
*/
// /////////////////////////////////////////////////////////////////////

#include "gpluto.hpp"
#if RADIATION
#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Local implicit step for radiation-matter interaction.
  Update both primitive and conserved quantities. Primitive and conserved
  fields must be consistent before every call to this function.

  \param[in,out] u           array of conservative variables
  \param[in,out] v           array of primitive variables
  \param[in,out] source      array of source terms
  \param[in,out] flag        array of flags
  \param[in]     dt          time step
  \param[in]     x1          x1 cell-center coordinate
  \param[in]     x2          x2 cell-center coordinate
  \param[in]     x3          x3 cell-center coordinate
  \param[in]     run_config  run-time configuration
*/
////////////////////////////////////////////////////////////////////////
uint32_t RadStep (double *u, double *v, double *source, uint32_t *flag,
                  double dt, double x1, double x2, double x3,
                  RunConfig &run_config)
{
  int i, j, m ;
  double err, gamma ;
  uint32_t radstep_err = 0 ;

  const int comps = RADIATION_NEQS - 1 ;
  Rad_data rad_data;

  double x[RADIATION_NEQS], dx[RADIATION_NEQS], mf[RADIATION_NEQS];
  double J[RADIATION_NEQS][RADIATION_NEQS];
  State state;
  STATE_SET(state, v, u, NULL, x1, x2, x3);
  state.flag = flag;

  // Set flag and time step
  rad_data.flag = flag ;
  rad_data.dt = dt ;
  rad_data.run_config = &run_config ;

  // Store primitive and conserved variables
  rad_data.pv = v ;
  rad_data.cv = u ;

  // Dispatch to the NR-specific Newton solver when requested
  #if RADIATION_IMPLICIT_NR
  return RadImplicitNR(u, v, source, flag, dt, x1, x2, run_config);
  #endif

  // Set initial energies
  #if RADIATION_NR
    // NR conservation law: Ttot[0] = E_gas + (c/c_r)*E_rad
    const double cratio = run_config.rad_reduced_c / run_config.rad_c ;
  #endif

  #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
    // Store (Egas,mgas) in Rini
    rad_data.Rini[0] = rad_data.Ttot[0] = u[ENG] ;
    #if RADIATION_NR
    rad_data.Ttot[0] += u[ENR] / cratio ;
    #else
    rad_data.Ttot[0] += u[ENR] ;
    #endif
  #else
    // Store (Erad,Frad) in Rini and Rprev
    rad_data.Rini[0] = u[ENR] ;
    #if RADIATION_NR
    rad_data.Ttot[0] = u[ENG] + u[ENR] / cratio ;
    #else
    rad_data.Ttot[0] = u[ENR] + u[ENG] ;
    #endif
    #if RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD
    rad_data.Rprev[0] = u[ENR] ;
    #endif
  #endif

  // Set initial extra variables if full convergence is imposed
  #if (RADIATION_IMPL == RADIATION_NEWTON_RAD || RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD) \
    && RADIATION_FULL_CONVERGENCE == YES
  rad_data.exv_prev = v[PRS] ;
  #elif (RADIATION_IMPL == RADIATION_NEWTON_GAS || RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS ) \
      && RADIATION_FULL_CONVERGENCE == YES
  rad_data.exv_prev = v[ENR] ;
  #endif

  // Set fluxes and total momentum
  for (j=0; j<comps; j++ ){
    #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
      // Store (Egas,mgas) in Rini
      rad_data.Rini[j+1] = rad_data.Ttot[j+1] = u[MX1+j] ;
      #if RADIATION_NR
      rad_data.Ttot[j+1] += u[FR1+j] / run_config.rad_reduced_c ;
      #else
      rad_data.Ttot[j+1] += u[FR1+j] ;
      #endif
    #else
      // Store (Erad,Frad) in Rini and Rprev
      rad_data.Rini[j+1] = u[FR1+j] ;
      #if RADIATION_NR
      rad_data.Ttot[j+1] = u[MX1+j] + u[FR1+j] / run_config.rad_reduced_c ;
      #else
      rad_data.Ttot[j+1] = u[FR1+j] + u[MX1+j] ;
      #endif
      #if RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD
      rad_data.Rprev[j+1] = u[FR1+j] ;
      #endif
    #endif
  }

  // Initial guess for the iterated fields 
  #if RADIATION_IMPL == RADIATION_NEWTON_GAS \
    || RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
    gamma = u[RHO]/v[RHO] ;
    x[0] = v[PRS] ;
    #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
      rad_data.Rprev[0] = x[0] ; // Store gas pressure in Rprev
    #endif
    for (j = 0 ; j < comps ; j++ ) {
      #if RADIATION_NR
      x[j+1] = v[VX1+j] ;
      #else
      x[j+1] = gamma*v[VX1+j] ;
      #endif
      #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
        rad_data.Rprev[j+1] = x[j+1] ; // Store proper velocity in Rprev
      #endif
    }
  #elif RADIATION_IMPL == RADIATION_NEWTON_RAD
    x[0] = u[ENR] ;
    for (j = 0 ; j < comps ; j++ ) x[j+1] = u[FR1+j] ;
  #endif
  
// -----------------------------------------------------------------
//    Implicit step until convergence or maximum number of iterations
//    -----------------------------------------------------------------
  m = 0 ; err = 1.0 ;
  while ( err > run_config.rad_tolerance && m++ < run_config.rad_max_iter ){

    /*********************************
          Update iterated fields
      *********************************/
    #if RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD

      // Set coefficients of the system C.(E,F^i)^{(m+1)} == b   	
      RadFPMatrices<RADIATION_NEQS>(&rad_data, dx, J);
      
      // Update iterated fields and store them in x 
      if ( GaussianSolve<RADIATION_NEQS>(J, dx, x) ) {
        (*flag) |= FLAG_GAUSSIAN_SOLVE;
        radstep_err |= ERR_GAUSSIAN_SOLVE;
      }
      
    #elif RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
    
      // Compute source terms and update iterated fields 
      RadNewtonMinusF(&rad_data, x, mf);
      
    #else

      // Compute -F and the Jacobian 
      RadNewtonJacobian<RADIATION_NEQS> (x, mf, J, &rad_data) ; 

      // Solve the system J*dx == -F 
      if ( GaussianSolve<RADIATION_NEQS>(J, mf, dx) ){
        (*flag) |= FLAG_GAUSSIAN_SOLVE;
        radstep_err |= ERR_GAUSSIAN_SOLVE;
        for (j=0; j < RADIATION_NEQS ; j++ ) dx[j] = mf[j];
      }

      // Update iterated fields 
      for (j=0; j < RADIATION_NEQS ; j++ ) x[j] += dx[j];

    #endif
    
    /*********************************
             Error calculation
      *********************************/
    #if RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD
    
      // Update conserved variables 
      u[ENR] = x[0] ;
      #if RADIATION_NR
      u[ENG] = rad_data.Ttot[0] - x[0] / cratio ;
      #else
      u[ENG] = rad_data.Ttot[0] - x[0] ;
      #endif
      for (j=0; j<comps; j++){
        u[FR1+j] = x[j+1] ;	
        #if RADIATION_NR
        u[MX1+j] = rad_data.Ttot[j+1] - x[j+1] / run_config.rad_reduced_c ;
        #else
        u[MX1+j] = rad_data.Ttot[j+1] - x[j+1] ;
        #endif
      }
      
      uint32_t local_err = ConsToPrim (state, run_config);
      radstep_err |= local_err;

      // Compute relative differences
      err = RadErr(v, NULL, &rad_data) ;
      
    #elif RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS

      // Compute relative differences 
      err = RadErr(x, NULL, &rad_data) ;
    
    #else

      // Compute relative differences 
      err = RadErr(x, dx, &rad_data) ;

    #endif

  }//End of iterations   

  // Final update if needed 
  #if RADIATION_IMPL == RADIATION_NEWTON_GAS
    RadIterToPrim (x, &rad_data);
    // if ( RadIterToPrim (x, &rad_data) ) QUIT_PLUTO(1);
  #elif RADIATION_IMPL == RADIATION_NEWTON_RAD
    u[ENR] = x[0] ;
    #if RADIATION_NR
    u[ENG] = rad_data.Ttot[0] - x[0] / cratio ;
    #else
    u[ENG] = rad_data.Ttot[0] - x[0] ;
    #endif
    for (j=0; j<comps; j++){
      u[FR1+j] = x[j+1] ;	
      #if RADIATION_NR
      u[MX1+j] = rad_data.Ttot[j+1] - x[j+1] / run_config.rad_reduced_c ;
      #else
      u[MX1+j] = rad_data.Ttot[j+1] - x[j+1] ;
      #endif
    }
    uint32_t local_err = ConsToPrim (state, run_config);
    radstep_err |= local_err;
  #endif

  // Check number of iterations
  if ( m > run_config.rad_max_iter ) {    radstep_err |= ERR_TOO_MANY_ITERATIONS;
  }

  // Compute and store source terms
  #if RADIATION_IMEX_SSP2
  RadSourceFunction(v, source, run_config);
  #endif

  // Update conserved fields if needed
  #if RADIATION_IMPL == RADIATION_NEWTON_GAS
  PrimToCons(v, u, run_config);
  #endif

  return radstep_err ;
}
#pragma omp end declare target
#endif

#if RADIATION_NR
#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Perform the implicit step for the nonrelativistic radiation module.
  Update both primitive and conserved quantities. Primitive and conserved
  fields must be consistent before every call to this function.

  \param[in,out] u           array of conservative variables
  \param[in,out] v           array of primitive variables
  \param[in,out] source      array of source terms
  \param[in,out] flag        array of flags
  \param[in]     dt          time step
  \param[in]     run_config  runtime configuration

  \b References
   - Melon Fuksman, J. D., & Klahr, H. 2022, ApJ, 936, 16. (Appendix B)
*/
////////////////////////////////////////////////////////////////////////
uint32_t RadImplicitNR (double *u, double *v, double *source, uint32_t *flag,
                        double dt, double x1, double x2, RunConfig &run_config)
{
  int j, m ;
  double abs_op, scat_op, tot_op, s, sf, err, cratio ;
  double mtot[3], frad[3], vel[3], vel0[3], vel2, ekin, etot, etek, erad,
         rho0, A, B, C, gmm1, x, x0, xsq, xaux, xaux2, y, y0, y2, f, fp ;
  uint32_t radstep_err = 0 ;

  #if EOS == IDEAL
  gmm1 = run_config.gamma - 1.0;
  #endif
  cratio = run_config.rad_reduced_c / run_config.rad_c;

  #if RADIATION_VAR_OPACITIES == NO
  abs_op  = run_config.rad_absorption_coeff;
  scat_op = run_config.rad_scattering_coeff;
  tot_op  = abs_op + scat_op;
  #endif

  // Get modified total energy and momentum
  ekin = 0.5*(u[MX1]*u[MX1]+u[MX2]*u[MX2]+u[MX3]*u[MX3])/u[RHO];
  #if PHYSICS == MHD
  ekin += 0.5*(u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3]);
  #endif
  etot = u[ENG] + u[ENR]/cratio ;
  #if IRRADIATION && !RADIATION_IMEX_SSP2
  etot -= dt*u[FIR];
  #endif
  etek = etot - ekin ;
  erad = u[ENR];
  for (j=0; j<3; j++){
    mtot[j] = u[MX1+j] + u[FR1+j]/run_config.rad_reduced_c ;
    vel0[j]  = v[VX1+j] ;
  }

  rho0 = u[RHO] ;

  #if RADIATION_IMPL == RADIATION_NEWTON_NR_GAS
  // x = -(E_gas - ekin)/rho0, first guess uses Er = Er(n)
  x = -(u[ENG]-ekin)/rho0;
  #if RADIATION_FULL_CONVERGENCE
  y = erad/rho0;
  #endif
  #else
  // x = Er/rho, first guess uses Er = Er(n)
  x = erad/rho0;
  #if RADIATION_FULL_CONVERGENCE
  y = -(u[ENG]-ekin)/rho0;
  #endif
  #endif

  // Solve A*x^4 + B*x + C == 0 with Newton's method
  m = 0; err = 1.;
  while(err > run_config.rad_tolerance && m++ < run_config.rad_max_iter ){

    // Get opacities
    #if RADIATION_VAR_OPACITIES
    UserDefOpacities (v, &abs_op, &scat_op, run_config);
    tot_op = abs_op + scat_op ;
    #endif

    // Get polynomial coefficients independent of ekin
    // (computed only once when opacities are constant)
    #if RADIATION_VAR_OPACITIES == NO
    if(m==1){
    #endif
      double dt_rho_c = dt*rho0*run_config.rad_reduced_c;
      s = dt_rho_c*abs_op ;
      sf = dt_rho_c*tot_op + 1. ;
      A = s*run_config.rad_constant
           *pow(run_config.rad_ideal_gas_constant*gmm1, 4);
      s += 1. ;
      B = -s*cratio*rho0 ;
    #if RADIATION_VAR_OPACITIES == NO
    }
    #endif

    // Update momentum and kinetic energy
    #if RADIATION_VAR_OPACITIES
    for (j=0; j<3; j++){
      frad[j]  = u[FR1+j]/sf ;
      u[MX1+j] = mtot[j] - frad[j]/run_config.rad_reduced_c ;
      vel[j]   = u[MX1+j]/u[RHO] ;
    }
    ekin = 0.5*(u[MX1]*u[MX1]+u[MX2]*u[MX2]+u[MX3]*u[MX3])/u[RHO] ;
    vel2 = ekin/u[RHO] ;
    #if PHYSICS == MHD
    ekin += 0.5*(u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3]);
    #endif
    etek = etot - ekin ;
    #endif

    // Compute remaining polynomial coefficient
    C = erad + B/rho0*etek ;

    // Update x guess with Newton's method
    x0 = x ;
    xsq = x*x ;
    #if RADIATION_IMPL == RADIATION_NEWTON_NR_GAS
    f  = A*xsq*xsq + B*x + C ;
    fp = 4.*A*xsq*x + B ;
    #else
    xaux  = x/cratio - etek/rho0 ;
    xaux2 = xaux*xaux ;
    f  = erad - s*rho0*x + A*xaux2*xaux2 ;
    fp = - s*rho0 + 4.*A*xaux2*xaux/cratio ;
    #endif
    x = x - f/fp ;

    // Compute pressure error
    double dx2 = x-x0 ;
    dx2 *= dx2 ;
    err = (xsq > 1e-40) ? dx2/xsq : dx2/1e-40 ;

    #if RADIATION_FULL_CONVERGENCE
    y0 = y ;
    #if RADIATION_IMPL == RADIATION_NEWTON_NR_GAS
    // Er/rho
    y = (etek/rho0 + x)*cratio ;
    #else
    // Specific internal energy
    y = etot/rho0 - x/cratio ;
    #endif
    y2 = y*y ;
    dx2 = y-y0 ;
    dx2 *= dx2 ;
    err += (y2 > 1e-40) ? dx2/y2 : dx2/1e-40 ;
    #endif

    // Add velocity error (except for constant opacities)
    #if RADIATION_VAR_OPACITIES
    for (j=0; j<3; j++){
      dx2 = vel[j]-vel0[j] ;
      dx2 *= dx2 ;
      err += (vel2 > 1e-40) ? dx2/vel2 : dx2/1e-40 ;
      vel0[j] = vel[j] ;
    }
    #endif

    // Update pressure for opacities calculation
    #if RADIATION_VAR_OPACITIES
    #if RADIATION_IMPL == RADIATION_NEWTON_NR_GAS
    v[PRS] = -gmm1*x*rho0;
    #else
    v[PRS] = gmm1*(etek - x*rho0/cratio);
    #endif
    if (v[PRS] < 0.0){
      radstep_err |= ERR_NEGATIVE_PRS;
      v[PRS] = run_config.small_prs;
    }
    #endif

  }

  // Check number of iterations
  if ( m > run_config.rad_max_iter ) {
    (*flag) |= FLAG_RAD_MAX_ITER ;
    radstep_err |= ERR_TOO_MANY_ITERATIONS;
  }

  // Update radiation flux and momentum
  s = 1. + dt*u[RHO]*run_config.rad_reduced_c*tot_op ;
  for (j=0; j<3; j++){
    mtot[j] = u[MX1+j] + u[FR1+j]/run_config.rad_reduced_c ;
    u[FR1+j] /= s ;
    u[MX1+j] = mtot[j] - u[FR1+j]/run_config.rad_reduced_c ;
  }

  // Update radiation and gas energy
  #if RADIATION_IMPL == RADIATION_NEWTON_NR_GAS
  ekin = 0.5*(u[MX1]*u[MX1]+u[MX2]*u[MX2]+u[MX3]*u[MX3])/u[RHO];
  #if PHYSICS == MHD
  ekin += 0.5*(u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3]);
  #endif
  u[ENG] = ekin - x*rho0 ;
  u[ENR] = (etot - u[ENG])*cratio ;
  #else
  u[ENR] = x*rho0 ;
  u[ENG] = (etot - u[ENR]/cratio) ;
  #endif

  // Update primitive fields
  {
    State state;
    STATE_SET(state, v, u, NULL, x1, x2, 0.0);
    state.flag = flag;
    uint32_t local_err = ConsToPrim(state, run_config);
    radstep_err |= local_err;
  }

  // Compute and store source terms if needed
  #if RADIATION_IMEX_SSP2
  RadSourceFunction(v, source, run_config);
  #endif

  return radstep_err ;
}
#pragma omp end declare target
#endif // RADIATION_NR

////////////////////////////////////////////////////////////////////////
/*!
  Apply the local radiation implicit step over the active 3D domain.

  \param[in,out] data  Pointer to the simulation data structure.
  \param[out]    S     Source-term array used by IMEX-SSP2.
  \param[in]     grid  Pointer to the grid structure.
  \param[in]     dt    Time step.
*/
////////////////////////////////////////////////////////////////////////
void RadStep3D(Data *data, Ary4D &S, Grid *grid, double dt)
{
  int i,j,k,nv;
  static RBox lbox;
  static bool first_call = true;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag

  Ary4D &U = data->Uc;
  Ary4D &V = data->Vc;
  Ary3Duint32_t &flag = data->flag;
  RunConfig &run_config = *(data->run_config);

  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];

  if (first_call){
    RBoxDefine(grid->lbeg[IDIR], grid->lend[IDIR],
               grid->lbeg[JDIR], grid->lend[JDIR],
               grid->lbeg[KDIR], grid->lend[KDIR], &lbox);
    #pragma acc enter data copyin(lbox)
    #pragma omp target enter data map(to: lbox)
    first_call = false;
  }

  mynvtxstart_("RadStep3D",RED);

  double v[NVAR], u[NVAR], s[NVAR];

  #if RADIATION_IMEX_SSP2
  #pragma omp target teams distribute parallel for collapse(3) private(v, u, s)
  #pragma acc parallel loop collapse(3) \
          present(data, grid, run_config, S, lbox) \
          private(v[:NVAR], u[:NVAR], s[:NVAR]) copy(err_kernel)
  #else
  #pragma omp target teams distribute parallel for collapse(3) private(v, u, s)
  #pragma acc parallel loop collapse(3) \
          present(data, grid, run_config, lbox) \
          private(v[:NVAR], u[:NVAR], s[:NVAR]) copy(err_kernel)
  #endif
  BOX_LOOP(lbox, i,j,k){

    data->err(i,j,k) = 0;

    #pragma acc loop seq
    NVAR_LOOP(nv) {
      u[nv] = U(i,j,k,nv);
      v[nv] = V(i,j,k,nv);
      if (isnan(u[nv]) || isnan(v[nv])) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_NAN, ERR_QUIT);
      }
    }

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    ConsToOrtho(u, metric);
    #endif

    uint32_t local_err = RadStep(u, v, s, &flag(i,j,k), dt,
                                 x1[i], x2[j], x3[k], run_config);

    #if GEOMETRY != CARTESIAN
    OrthoToCons(u, metric);
    #endif
    
    ATOMIC_ERR_UPDATE(err_kernel, local_err);
    KERNEL_ERR_UPDATE(data->err(i,j,k), local_err);
    if (local_err & ERR_CONS2PRIM_FAIL) {
      ATOMIC_ERR_UPDATE(err_kernel, ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_QUIT);
    }

    #pragma acc loop seq
    NVAR_LOOP(nv) {
      U(i,j,k,nv) = u[nv];
      V(i,j,k,nv) = v[nv];
      #if RADIATION_IMEX_SSP2
      S(i,j,k,nv) = s[nv];
      #endif
    }

  }

  #if WARNING_MESSAGES > 0
  #pragma acc data present(data,grid,lbox)
  {
    Where ("RadStep3D()", data, err_kernel, &lbox, grid);
  }
  #endif

  #if WARNING_MESSAGES >= 0
  int quit_condition = (err_kernel & ERR_QUIT);

  #if FALLBACK_MODE != NO
  quit_condition *= (data->integration_mode != DEFAULT);
  Fallback_Trigger(data->err, data->flag, &lbox);
  #endif

  if (quit_condition) QuitPLUTO("from RadStep3D()", err_kernel);
  #endif
  
  mynvtxstop_();
}

////////////////////////////////////////////////////////////////////////
/*!
  Apply the local radiation implicit step over the active 3D domain.

  \param[in,out] data  Pointer to the simulation data structure.
  \param[in]     grid  Pointer to the grid structure.
  \param[in]     dt    Time step.
*/
////////////////////////////////////////////////////////////////////////
void RadStep3D(Data *data, Grid *grid, double dt)
{
  Ary4D S;
  RadStep3D(data, S, grid, dt);
}
