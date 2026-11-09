
// /////////////////////////////////////////////////////////////////////
/*!
  \file 
  \brief  Source functions for the radiation transport equations.
  
  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 8, 2026

*/
// /////////////////////////////////////////////////////////////////////

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Compute the Eulerian frame components of the source function G^\mu,
  that accounts for the interaction between radiation and matter.

  \param[in]  primvar     vector of primitive fields
  \param[out] G           vector that stores the value of G^\mu
  \param[in]  run_config  runtime configuration for opacities
*/
////////////////////////////////////////////////////////////////////////
inline void RadSourceFunction(double *primvar, double *G, RunConfig &run_config)
{
  int i, j, nv ;
  const int comps = RADIATION_NEQS - 1 ;
  double B, rho, T ;
  double u[3], uD[3], uF, q ;

  // Set opacities
  #if RADIATION_VAR_OPACITIES == YES
    double abs_op, scat_op, tot_op ;
    UserDefOpacities (primvar, &abs_op, &scat_op, run_config);
    tot_op = abs_op + scat_op ;
  #else
    const double abs_op = run_config.rad_absorption_coeff;
    const double scat_op = run_config.rad_scattering_coeff;
    const double tot_op = abs_op + scat_op ;
  #endif

  NVAR_LOOP(nv) G[nv] = 0.;

  #if RADIATION_NR

  rho = primvar[RHO] ;
  T   = GetTemperature(rho, primvar[PRS], run_config);
  B   = Blackbody(T, run_config) ;

  #if RADIATION_IMPLICIT_NR
  G[0] = abs_op * rho * (primvar[ENR] - B) ;
  for (i=0; i<comps; i++) G[i+1] = tot_op * rho * primvar[FR1+i] ;
  return ;
  #else

  double s, comovEr, comovFr[3] ;

  for (i=0; i<comps; i++) u[i] = primvar[VX1+i] / run_config.rad_c ;

  uF = 0.;
  for (i=0; i<comps; i++ ){
    uD[i] = 0.;
    uF += u[i]*primvar[FR1+i] ;
    for (j=0; j<comps; j++ )
      uD[i] += u[j] * EddTensor(primvar, FR1+j, FR1+i);
  }

  comovEr = primvar[ENR] - 2.*uF ;
  for (i=0; i<comps; i++)
    comovFr[i] = primvar[FR1+i] - primvar[ENR] * (u[i] + uD[i]) ;

  q = abs_op * rho * (comovEr - B) ;
  s = tot_op * rho ;
  uF = 0.; for (i=0; i<comps; i++) uF += u[i] * comovFr[i] ;

  G[0] = q + s * uF ;
  for (i=0; i<comps; i++) G[i+1] = q * u[i] + s * comovFr[i] ;
  #endif

  #else

  double u2, gamma, gamma2, D, uuD ;
  double rhogamma ;

  // Compute gamma, gamma^2, u and u^2
  gamma = primvar[VX1]*primvar[VX1]
        + primvar[VX2]*primvar[VX2]
        + primvar[VX3]*primvar[VX3] ;
  u2 = gamma ;
  gamma2 = 1.0/(1.0-gamma) ;
  u2 *= gamma2 ;
  gamma = sqrt(gamma2) ;
  for (i=0; i<comps; i++) u[i] = gamma * primvar[VX1+i] ;

  //Compute products involving the proper velocity
  uuD = 0.;
  uF = 0.;
  for (i=0; i<comps; i++ ){
    uD[i] = 0.;
    uF += u[i]*primvar[FR1+i] ;
    for (j=0; j<comps; j++ ){
      D = EddTensor(primvar, FR1+j, FR1+i);
      uD[i] += u[j]*D;
      uuD += u[i]*u[j]*D ;
    }
  }

  // Compute some useful quantities
  rho = primvar[RHO] ;
  rhogamma = rho*gamma ;
  T = GetTemperature(rho, primvar[PRS], run_config);
  B = Blackbody(T, run_config) ;
  q = - abs_op*rho*B ;

  // Compute source function in the Eulerian frame
  G[0] = q * gamma
       + rhogamma * ( abs_op - scat_op*( u2 + uuD ) ) * primvar[ENR]
       - rho * ( abs_op  - scat_op * (u2 + gamma2) ) * uF ;

	for (i=0; i<comps; i++){
    G[i+1] = q * u[i]
           - rho * ( tot_op * uD[i]
                   + scat_op * u[i] * (gamma2 + uuD) ) * primvar[ENR]
           + rhogamma * ( tot_op * primvar[FR1+i]
                        + scat_op * 2.0 * uF * u[i] ) ;
  }

  #endif

}
#pragma omp end declare target

#if RADIATION_NR
#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Explicitly integrate order v/c source terms.

  \param[in]     primvar     vector of primitive variables
  \param[in,out] rhs         vector of right-hand sides
  \param[in]     dt          time step
  \param[in]     run_config  runtime configuration for reduced_c
*/
////////////////////////////////////////////////////////////////////////
inline void AddRelTerms (double *primvar, double *rhs, double dt,
                         RunConfig &run_config)
{
  int i, j ;
  double u[3], uD[3], uF, uFc, G[4] ;
  double B, rho, q, s, comovEr, comovFr[3] ;
  const int comps = RADIATION_NEQS - 1 ;
  const double radC        = run_config.rad_c ;
  const double reducedC    = run_config.rad_reduced_c ;

  // Set opacities
  double abs_op, scat_op, tot_op ;
  #if RADIATION_VAR_OPACITIES == YES
  UserDefOpacities (primvar, &abs_op, &scat_op, run_config);
  tot_op = abs_op + scat_op ;
  #else
  abs_op  = run_config.rad_absorption_coeff;
  scat_op = run_config.rad_scattering_coeff;
  tot_op  = abs_op + scat_op ;
  #endif

  // Compute beta = v/c
  for (i=0; i<comps; i++) u[i] = primvar[VX1+i] / radC ;

  // Compute products involving beta
  uF = 0.;
  for (i=0; i<comps; i++ ){
    uD[i] = 0.;
    uF += u[i]*primvar[FR1+i] ;
    for (j=0; j<comps; j++ )
      uD[i] += u[j] * EddTensor(primvar, FR1+j, FR1+i);
  }

  // Compute comoving radiation fields
  comovEr = primvar[ENR] - 2.*uF ;
  for (i=0; i<comps; i++)
    comovFr[i] = primvar[FR1+i] - primvar[ENR] * (u[i] + uD[i]) ;

  // Compute some auxiliary quantities
  rho = primvar[RHO] ;
  B = Blackbody( GetTemperature(rho, primvar[PRS], run_config), run_config ) ;
  q = abs_op * rho * (comovEr - B) ;
  s = tot_op * rho ;
  uFc = 0.; for (i=0; i<comps; i++) uFc += u[i] * comovFr[i] ;

  // Compute relativistic component of source function in the lab frame
  G[0] = -2.*abs_op*rho*uF + s*uFc ;
  for (i=0; i<comps; i++) G[i+1] = q*u[i] - s*primvar[ENR]*(u[i] + uD[i]) ;

  // Update rhs
  rhs[ENG] += dt*radC*G[0] ;
  rhs[ENR] -= dt*reducedC*G[0] ;
  for(i=0; i<comps; i++){
    rhs[MX1+i] += dt*G[i+1] ;
    rhs[FR1+i] -= dt*G[i+1]*reducedC ;
  }

}
#pragma omp end declare target
#endif
