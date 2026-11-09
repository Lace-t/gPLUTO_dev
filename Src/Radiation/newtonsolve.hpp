// /////////////////////////////////////////////////////////////////////
/*!
  \file 
  \brief  Newton's method, main functions

  Functions used to compute the Jacobian matrix during Newton's method, and
  to update the iterated fields in case RADIATION_FIXEDPOINT_GAS == YES.
  
  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    Dec 12, 2019
*/
// /////////////////////////////////////////////////////////////////////

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Conver iterated to primitive fields if iterations are carried on MHD fields

  \param[in]     x         vector of iterated fields
  \param[in,out] rad_data  pointer to rad_data structure that contains
                           the primitive fields
  \return Returns false if success, and true if prs == NaN is found.
*/
////////////////////////////////////////////////////////////////////////
inline bool RadIterToPrim (double x[], Rad_data *rad_data)
{
  int i ;
  const int comps = RADIATION_NEQS - 1 ;
  State state;
  double * v ;
  double u[3] ; 
  uint32_t * flag = rad_data-> flag ;
  double u2, gamma, gamma2, gm1 ;
  double rhogamma, rhoh, vel2, vB, B2 ;
  double th, th2, wt ;
  #if EOS == IDEAL
   const double gm = rad_data->run_config->gamma ;
   const double gmmr = gm/(gm - 1.0) ;
  #endif

  // Get primitive variables
  v = rad_data-> pv ;

  // Set pressure 
  v[PRS] = x[0];
  STATE_SET(state, v, rad_data->cv, NULL, 0.0, 0.0, 0.0);
  state.flag = flag;

  // Check p != NaN 
  if (v[PRS] != v[PRS]) {
    (*flag) |= ERR_NAN;
    #if WARNING_MESSAGES
    printf("! RadIterToPrim: NaN found while setting pressure.\n");
    #endif
    return true ;
  } 

  // Check p > 0 
  if (v[PRS] < 0.0) {
    (*flag) |= ERR_NEGATIVE_PRS;
    #if WARNING_MESSAGES
    printf("! RadIterToPrim: negative pressure (%8.2e) during implicit step.\n", v[PRS]);
    #endif
    v[PRS]   = 1.e-20;
  } 
  
  #if RADIATION_NR

  u2 = 0.; vel2 = 0.;
  for (i=0; i<comps; i++){
    v[VX1+i] = x[i+1];
    u[i] = x[i+1] / rad_data->run_config->rad_c;
    u2   += u[i]*u[i];
    vel2 += x[i+1]*x[i+1];
  }
  rad_data->u2    = u2;
  rad_data->gamma  = 1.0;
  rad_data->gamma2 = 1.0;

  #if EOS == IDEAL
  {
    double kinb2;
    #if PHYSICS == HD
      kinb2 = 0.5*v[RHO]*vel2;
    #elif PHYSICS == MHD
      kinb2 = 0.5*(v[RHO]*vel2 + v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
    #endif
    v[ENR] = (rad_data->Ttot[0] - kinb2 - v[PRS]/(gm - 1.0))
           * rad_data->run_config->rad_reduced_c / rad_data->run_config->rad_c;
  }
  #endif

  for (i=0; i<comps; i++)
    v[FR1+i] = rad_data->run_config->rad_reduced_c
             * (rad_data->Ttot[i+1] - v[RHO]*v[VX1+i]);

  #else

  // Set 4-velocity and Lorentz factor
  u2 = 0. ;
  for (i=0; i<comps; i++){
    u[i] = x[i+1] ;
    u2 += u[i]*u[i];
  }
  gamma2 = u2 + 1.0 ;
  gamma = sqrt(gamma2) ;
  gm1 = 1.0/gamma ;

  // Store gamma, gamma2 and u2
  rad_data-> u2 = u2 ;
  rad_data-> gamma = gamma ;
  rad_data-> gamma2 = gamma2 ;

  // Set coordinate velocity
  for (i=0; i<comps; i++) v[VX1+i] = gm1 * u[i] ;

  // Set density
  rhogamma = rad_data->cv[RHO] ;
  v[RHO] = gm1 * rhogamma ;

  // Set enthalpy density
  #if EOS == IDEAL
    rhoh = v[RHO] + gmmr*v[PRS] ;
  #elif EOS == TAUB
    rhoh = 2.5*v[PRS] + sqrt(2.25*v[PRS]*v[PRS] + v[RHO]*v[RHO]);
  #endif

  // Set radiation fields
  #if PHYSICS == RHD

  v[ENR] = rad_data-> Ttot[0] - rhoh*gamma2 + v[PRS] ;
  for (i=0; i<comps; i++)
    v[FR1+i] = rad_data-> Ttot[i+1] - rhoh*gamma*u[i] ;

  #elif PHYSICS == RMHD

  vel2  = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  vB    = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  B2    = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];

    #if RMHD_REDUCED_ENERGY == YES
      #if EOS == IDEAL
      v[ENR] = rad_data-> Ttot[0] - v[PRS]*(gamma2*gmmr - 1.0)
             - rhogamma*u2/(gamma+1.0) - 0.5*(B2*(1.0 + vel2) - vB*vB) ;
      #elif EOS == TAUB
      th = v[PRS]/v[RHO];
      th2 = th*th;
      v[ENR] = rad_data-> Ttot[0] - v[PRS]*(gamma2*2.5 - 1.0)
             - rhogamma*(2.25*th2*gamma2 + u2)/(gamma*sqrt(1.0 + 2.25*th2) + 1.0)
             - 0.5*(B2*(1.0 + vel2) - vB*vB) ;
      #endif
    #else
      v[ENR] = rad_data-> Ttot[0] - rhoh*gamma2 + v[PRS]
             - 0.5*(B2*(1.0 + vel2) - vB*vB) ;
    #endif

  wt = rhoh*gamma2 + B2 ;
  for (i=0; i<comps; i++)
    v[FR1+i] = rad_data-> Ttot[i+1] - wt*v[VX1+i] + vB*v[BX1+i] ;

  #endif

  #endif // RADIATION_NR

  // Check Er > 0 
  if (v[ENR] < 0.0) {
    (*flag) |= ERR_NEGATIVE_ENR;
    #if WARNING_MESSAGES
    printf("! RadIterToPrim: negative radiation energy (%8.2e) during implicit step.\n", v[ENR]);
    #endif
    v[ENR]   = rad_data->run_config->small_enr;
  } 

  return false ;
}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  If Newton's method is used, compute and store in mf the components
  of -F, being F the function whose roots are being looked for. Otherwise,
  if RADIATION_FIXEDPOINT_GAS, perform a fixed-point iteration.

  \param[in,out] rad_data  pointer to rad_data structure
  \param[in]     x         vector of iterated fields
  \param[in,out] mf        vector that stores the value of -F(x)
*/
////////////////////////////////////////////////////////////////////////
inline void RadNewtonMinusF(Rad_data * rad_data, double x[], double mf[])
{
  int i, j ;
  const int comps = RADIATION_NEQS - 1 ;
  State state;
  double u[3], u2, gamma, gamma2 ;
  double D, uD[3], uuD, uF ;
  double B, rhogamma, q, dt, T ;
  RunConfig &run_config = *(rad_data->run_config);

  #if RADIATION_VAR_OPACITIES
    double abs_op, scat_op, tot_op ;
  #else
    const double abs_op = rad_data->run_config->rad_absorption_coeff;
    const double scat_op = rad_data->run_config->rad_scattering_coeff;
    const double tot_op = abs_op + scat_op ;
  #endif

  double *v, *consv ;

  v = rad_data-> pv;
  consv = rad_data-> cv;  

  #if RADIATION_NR
  dt = rad_data->run_config->rad_reduced_c * rad_data->dt ;
  #else
  dt = rad_data-> dt ;
  #endif

  #if RADIATION_IMPL == RADIATION_NEWTON_GAS || RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS

    // Convert iterated fields to primitive fields 
    RadIterToPrim (x, rad_data);

    // Compute 4-velocity and derived quantities 
    gamma = rad_data-> gamma;
    gamma2 = rad_data-> gamma2;
    u2 = rad_data-> u2;

  #else

    // Compute conserved fields
    consv[ENR] = x[0];
    #if RADIATION_NR
    consv[ENG] = rad_data->Ttot[0] - x[0] / run_config.rad_reduced_c * run_config.rad_c ;
    #else
    consv[ENG] = rad_data->Ttot[0] - x[0];
    #endif
    for (i=0; i<comps; i++){
      consv[FR1+i] = x[i+1] ;
      #if RADIATION_NR
      consv[MX1+i] = rad_data->Ttot[i+1] - x[i+1] / run_config.rad_reduced_c ;
      #else
      consv[MX1+i] = rad_data->Ttot[i+1] - x[i+1] ;
      #endif
    }

    // Convert conserved to primitive fields 
    STATE_SET(state, rad_data->pv, rad_data->cv, NULL, 0.0, 0.0, 0.0);
    state.flag = rad_data->flag;
    ConsToPrim(state, run_config);

    // Compute 4-velocity and derived quantities
    #if RADIATION_NR
    u2 = 0.;
    gamma = 1.0;
    gamma2 = 1.0;
    #else
    gamma = consv[RHO]/v[RHO] ;
    gamma2 = gamma*gamma ;
    u2 = gamma2 - 1.0 ;
    #endif

  #endif

  // Compute u[i] = beta (NR) or proper velocity (relativistic)
  #if RADIATION_NR
  u2 = 0. ;
  for (i=0; i<comps; i++) {
    u[i] = v[VX1+i] / run_config.rad_c ;
    u2 += u[i]*u[i] ;
  }
  #else
  for (i=0; i<comps; i++) u[i] = gamma*v[VX1+i] ;
  #endif

  // Set opacities 
  #if RADIATION_VAR_OPACITIES == YES
    UserDefOpacities (v, &abs_op, &scat_op, run_config);
    tot_op = abs_op + scat_op ; 
  #endif

  // Compute products involving the proper velocity 
  uuD = 0.;
  uF = 0.;
  for (i=0; i<comps; i++ ){
    uD[i] = 0.;
    uF += u[i]*v[FR1+i] ;
    for (j=0; j<comps; j++ ){
      D = EddTensor(v, FR1+j, FR1+i);
      uD[i] += u[j]*D;
      uuD += u[i]*u[j]*D ;
    }
  }
  
  // Compute some useful quantities 
  rhogamma = consv[RHO] ;
  T = GetTemperature(v[RHO], v[PRS], *(rad_data->run_config));
  B = Blackbody(T, *(rad_data->run_config)) ;
  q = - abs_op*v[RHO]*B ;

  // Compute -F[x] or update iterated fields if
  // RADIATION_FIXEDPOINT_GAS is used
  #if RADIATION_NR
  mf[0] = q
       + v[RHO] * ( abs_op - tot_op * (u2 + uuD) ) * v[ENR]
       - v[RHO] * ( abs_op - scat_op ) * uF ;
  #else
  mf[0] = q * gamma
       + rhogamma * ( abs_op - scat_op*( u2 + uuD ) ) * v[ENR]
       - v[RHO] * ( abs_op  - scat_op * (u2 + gamma2) ) * uF ;
  #endif

  #if RADIATION_IMPL != RADIATION_FIXEDPOINT_GAS
    mf[0] = rad_data->Rini[0] - v[ENR] - dt*mf[0] ;
  #else
    consv[ENG] = rad_data-> Rini[0] + dt*mf[0] ;
    #if RADIATION_NR
    consv[ENR] = (rad_data->Ttot[0] - consv[ENG])
               * run_config.rad_reduced_c / run_config.rad_c ;
    #else
    consv[ENR] = rad_data-> Ttot[0] - consv[ENG] ;
    #endif
  #endif

	for (i=0; i<comps; i++){
    #if RADIATION_NR
    mf[i+1] = q * u[i]
            - v[RHO] * ( tot_op * uD[i]
                       + scat_op * u[i] ) * v[ENR]
            + v[RHO] * ( tot_op * v[FR1+i]
                       - abs_op * 2.0 * uF * u[i] ) ;
    #else
    mf[i+1] = q * u[i]
            - v[RHO] * ( tot_op * uD[i]
                       + scat_op * u[i] * (gamma2 + uuD) ) * v[ENR]
            + rhogamma * ( tot_op * v[FR1+i]
                         + scat_op * 2.0 * uF * u[i] ) ;
    #endif

    #if RADIATION_IMPL != RADIATION_FIXEDPOINT_GAS
      mf[i+1] = rad_data->Rini[i+1] - v[FR1+i] - dt*mf[i+1] ;
    #else
      consv[MX1+i] = rad_data-> Rini[i+1] + dt*mf[i+1] ;
      #if RADIATION_NR
      consv[FR1+i] = run_config.rad_reduced_c
                   * (rad_data->Ttot[i+1] - consv[MX1+i]) ;
      #else
      consv[FR1+i] = rad_data-> Ttot[i+1] - consv[MX1+i] ;
      #endif
    #endif
  }
  
  #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
    // Update Erad for error computation 
    #if RADIATION_FULL_CONVERGENCE == YES
    rad_data-> exv = consv[ENR] ;
    #endif
    
    // Update iterated fields 
    STATE_SET(state, rad_data->pv, rad_data->cv, NULL, 0.0, 0.0, 0.0);
    state.flag = rad_data->flag;
    ConsToPrim(state, *(rad_data->run_config)) ;
    x[0] = v[PRS] ;
    #if RADIATION_NR
    for (i=0; i<comps; i++) x[i+1] = v[VX1+i] ;
    #else
    gamma = consv[RHO]/v[RHO] ;
    for (i=0; i<comps; i++) x[i+1] = gamma * v[VX1+i] ;
    #endif
  #endif

}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
template<int N>
////////////////////////////////////////////////////////////////////////
/*!
  Compute and store in J the Jacobian matrix (dF^i/dx^j) of the system F(x) == 0,
  calculated using difference quotients with small variations of the variables.
  Store in mf the current value of -F(x).

  \tparam        N        Size of the nonlinear system.
  \param[in]     x         vector of iterated fields
  \param[in,out] mf        vector that stores the value of -F(x)
  \param[in,out] J         Jacobian matrix
  \param[in,out] rad_data  pointer to rad_data structure, used in this case
                           to store the initial values of the radiation fields
                           as well as primitive fields
*/
////////////////////////////////////////////////////////////////////////
inline void RadNewtonJacobian (double x[], double mf[], double J[][N], Rad_data * rad_data)
{
  int i, j ;
  double y[N], mf2[N], dvars[N] ; 

  // - Define increments to compute derivatives -
  for (i=0 ; i < N ; i++ )
    dvars[i] = (fabs(x[i]) > 1e-20) ? -1e-4*x[i] : 1e-20 ;

  // - Store the current value of -F(x) -
  RadNewtonMinusF(rad_data,x,mf);

  // - Update auxiliary variables for error computation -
  #if RADIATION_FULL_CONVERGENCE == YES
    #if RADIATION_IMPL == RADIATION_NEWTON_GAS
    rad_data->exv = rad_data->pv[ENR] ; 
    #elif RADIATION_IMPL == RADIATION_NEWTON_RAD
    rad_data->exv = rad_data->pv[PRS] ;
    #endif
  #endif

  // - Compute and store the elements of the Jacobian -
  for (i=0 ; i < N ; i++ ) y[i] = x[i] ;

  for (j=0; j< N ; j++ ){
    y[j] += dvars[j] ;

    RadNewtonMinusF(rad_data,y,mf2);

    for (i=0 ; i < N ; i++ ) J[i][j] = ( mf[i] - mf2[i] ) / dvars[j] ;
    y[j] = x[j] ;
  }

}
#pragma omp end declare target
