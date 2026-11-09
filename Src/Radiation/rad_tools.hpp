// /////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Some functions used by the radiation module.
  
  Contains some functions used to compute radiation-matter interaction
  terms, and to solve systems of linear equations during the implicit step.
  
  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 9, 2026
*/
// /////////////////////////////////////////////////////////////////////

#pragma omp declare target
#pragma acc routine seq
template<int N>
////////////////////////////////////////////////////////////////////////
/*!
  Solve the system of linear equations given by C.x=b, using
  gaussian elimination with scaled partial pivoting (see e.g.
  Burden R. L. and Faires J. D., Numerical Analysis).

  \tparam     N  Dimension of the system.
  \param[in]  C  matrix of coefficients.
  \param[in]  b  vector of constant terms.
  \param[out] x  vector of solutions.

  \return Returns 1 if the system is solved, and 0 if no unique
          solution exists.
*/
////////////////////////////////////////////////////////////////////////
inline int GaussianSolve (double C[][N], double b[N], double x[N])
{
	int i, j, k, p;
	int nrow[N] ; // Row vector, used for permutations of rows.
	double m[N] ; // Vector of coefficients used during the elimination.
	double s[N] ; // Vector of scaling factors.
	double A[N][N+1] ; // Augmented matrix of the system.
	double max_fabsA, sum, scale, q;

  // -- Store C and b values in A --
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) A[i][j] = C[i][j];
		A[i][N] = b[i];
	}

  // -- Initialize row pointer --
	for (i=0; i<N; i++) nrow[i] = i;

  // -- Compute scaling factors --
	for (i=0; i<N; i++) {
		max_fabsA = 0. ;
		for (j=0; j<N; j++) {
			scale = fabs(A[i][j]) ;
			if ( scale > max_fabsA ) max_fabsA = scale;
		}
		if ( max_fabsA == 0 ) {
			return 1;
		} else {
			s[i] = max_fabsA ;	
		}
	}

  // -------------------------------------
  //    	Elimination process
  //  -------------------------------------
	for (i=0; i<N-1; i++) {

		p = i;
		max_fabsA = fabs(A[nrow[i]][i]) / s[nrow[i]];	

    // -- Get maximum quotient |A[nrow[i]][j]|/s[nrow[j]] for scaled pivoting --	
    for (j = i + 1; j < N; j++) {
      q = fabs(A[nrow[j]][i]) / s[nrow[j]];
      if (q > max_fabsA) {
          max_fabsA = q;
          p = j;
      }
    }

    // -- Check if solution is unique --
		if ( A[nrow[p]][i] == 0. ) return 1;

    // -- Switch rows --
		if ( nrow[i]!=nrow[p] ) {
			j = nrow[i] ;
			nrow[i] = nrow[p] ;
			nrow[p] = j ;
		}

    // -- Elimination --
		for (j=i+1; j<N; j++) {
			m[j] = A[nrow[j]][i]/A[nrow[i]][i] ;
			for (k=i; k<N+1; k++) {
				A[nrow[j]][k] -= m[j]*A[nrow[i]][k] ;
			}
		}

	}
  // -------------------------------------
  //    	End of elimination
  //  -------------------------------------

  // -- Check if solution is unique for the last time --
	if ( A[nrow[N-1]][N-1] == 0. ) return 1;
		
  // -- Backward substitution --
	x[N-1] = A[nrow[N-1]][N]/A[nrow[N-1]][N-1] ;

	for (i=N-2; i>=0; i-- ) {
		sum=0;
		for(j=i+1; j<N; j++) {
			sum += A[nrow[i]][j]*x[j] ;
		}
		x[i] = ( A[nrow[i]][N] - sum ) / A[nrow[i]][i] ;
	}

  // -- Output when no errors occur  --
	return 0;

}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
inline double EddTensor (double *v, int i, int j){
// --------------------------------------------------------
//		Compute the (i,j) component of the Eddington tensor
//		using the M1 closure, taking as input the state v.
//		(i,j) values must always be within FR1 and FR3.
//   --------------------------------------------------------

  double ni, nj, f2, chi_2, edt ;

  double flux_module = v[FR1]*v[FR1] + v[FR2]*v[FR2] + v[FR3]*v[FR3] ;

  if (flux_module < 1e-150 || v[ENR] < 1e-150 ) {
    if ( i==j ) return 0.33333333333333 ; 
    return 0.0 ;
  }

  flux_module = sqrt(flux_module) ;

  ni = v[i]/flux_module;
  nj = v[j]/flux_module ;		

	if( flux_module > v[ENR] ) return ni*nj ;

  f2 = flux_module/v[ENR] ;
  f2 = f2*f2 ;
  chi_2 = ( 1.5 + 2.0*f2 )/( 5.0 + 2.0*sqrt(4.0-3.0*f2) ) ;
  edt = (3.0*chi_2-0.5)*ni*nj ;

  if ( i==j ) edt += 0.5 - chi_2 ;

  return edt ;
}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Return the blackbody intensity corresponding to the input temperature.

  \param[in] temperature  Input temperature.
  \param[in] run_config   Run-time configuration for constants.

  \return Return (4*PI)B(T)=radiationConst*T^4.
*/
////////////////////////////////////////////////////////////////////////
inline double Blackbody (double temperature, const RunConfig &run_config)
{
  // const double radiationConst = 4.*CONST_sigma/CONST_c/(UNIT_VELOCITY*UNIT_VELOCITY)/UNIT_DENSITY ;
  const double radiationConst = run_config.rad_constant;
	double T4 = temperature*temperature;
	T4 = T4*T4;
	return radiationConst * T4;
}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Return the (ideal) gas temperature corresponding to given gas pressure
  and density.

  \param[in] rho         Input mass density.
  \param[in] prs         Input (gas) pressure.
  \param[in] run_config  Runtime configuration.

  \return Return T =  idealGasConst*prs/rho, where idealGasConst is the
          equilibrium ratio prs/(T*rho) for the considered ideal gas.
*/
////////////////////////////////////////////////////////////////////////
inline double GetTemperature (double rho, double prs, RunConfig &run_config)
{
  const double idealGasConst = run_config.rad_ideal_gas_constant;
	return idealGasConst * prs / rho;
}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
template<int N>
////////////////////////////////////////////////////////////////////////
/*!
  Compute the components of the matrix C and the column vector in eq.
  (52) of Takahashi et al. 2013, ApJ 772(2):127. Used if RADIATION_IMPL
  == RADIATION_FIXEDPOINT_RAD.

  \tparam     N         Size of the linear system.
  \param[in]  rad_data  pointer to rad_data structure
  \param[out] b         column vector
  \param[out] C         matrix C
*/
////////////////////////////////////////////////////////////////////////
inline void RadFPMatrices (Rad_data * rad_data, double b[], double C[][N])
{
  int i, j ;
  const int comps = N-1 ;
  double * primvar ;
  double u[3], u2, gamma, gamma2;
  double D, uD[3], uuD ;
  double rho, B, rhogamma, dt, T, mdt_rho, coeff1, coeff2 ;

  primvar = rad_data-> pv ;

  // -- Set opacities --
  #if RADIATION_VAR_OPACITIES == YES
    double abs_op, scat_op, tot_op ;
    UserDefOpacities (primvar, &abs_op, &scat_op, *(rad_data->run_config));
    tot_op = abs_op + scat_op ; 
  #else
    const double abs_op = rad_data->run_config->rad_absorption_coeff;
    const double scat_op = rad_data->run_config->rad_scattering_coeff;
    const double tot_op = abs_op + scat_op ;
  #endif

  // -- Get dt --
  #if RADIATION_NR
  dt = rad_data->run_config->rad_reduced_c * rad_data->dt ;
  #else
  dt = rad_data->dt ;
  #endif

  #if RADIATION_NR

  // -- Compute beta = v/c and beta^2 --
  u2 = 0.;
  for (i=0; i<comps; i++){
    u[i] = primvar[VX1+i] / rad_data->run_config->rad_c ;
    u2  += u[i]*u[i] ;
  }

  // -- Set products of beta and the Eddington tensor --
  uuD = 0.;
  for (i=0; i<comps; i++ ){
    uD[i] = 0.;
    for (j=0; j<comps; j++ ){
      D = EddTensor(primvar, FR1+j, FR1+i);
      uD[i] += u[j]*D;
      uuD += u[i]*u[j]*D ;
    }
  }

  // -- Set some constants used to calculate C and b  --
  rho = primvar[RHO] ;
  mdt_rho = -dt*rho ;
  T = GetTemperature(rho, primvar[PRS], *(rad_data->run_config));
  B = Blackbody(T, *(rad_data->run_config)) ;
  coeff1 = dt * abs_op * rho * B ;

  // -- Set coefficients of the system C.(E,F^i)=b  --
  b[0] = rad_data->Rini[0] + coeff1 ;
  C[0][0] = 1.0 + mdt_rho * ( scat_op*( u2 + uuD ) - abs_op );

  for (i=0; i<comps; i++){
    b[i+1] = rad_data->Rini[i+1] + coeff1 * u[i] ;

    C[0][i+1] = mdt_rho * u[i] * ( abs_op - scat_op ) ;
    C[i+1][0] = mdt_rho * ( tot_op * uD[i] + scat_op * u[i] ) ;

    for (j=i; j<comps; j++)
      C[i+1][j+1] = C[j+1][i+1] = mdt_rho * abs_op * 2.0 * u[i] * u[j] ;

    C[i+1][i+1] += 1.0 - mdt_rho * tot_op ;
  }

  #else

  // -- Set gamma, gamma^2, u and u^2 --
	gamma = primvar[VX1]*primvar[VX1]
        + primvar[VX2]*primvar[VX2]
				+ primvar[VX3]*primvar[VX3] ;
	u2 = gamma ;
	gamma2 = 1.0/(1.0-gamma) ;
	u2 *= gamma2 ;
  gamma = sqrt(gamma2) ;
	for (i=0; i<comps; i++) u[i] = gamma * primvar[VX1+i] ;

  // -- Set products of u and the Eddington tensor --
  uuD = 0.;
  for (i=0; i<comps; i++ ){
    uD[i] = 0.;
    for (j=0; j<comps; j++ ){
      D = EddTensor(primvar, FR1+j, FR1+i);
      uD[i] += u[j]*D;
      uuD += u[i]*u[j]*D ;
    }
  }

  // -- Set some constants used to calculate C and b  --
  rho = primvar[RHO] ;
  rhogamma = rad_data->cv[RHO] ;
  mdt_rho = -dt*rho ;
  T = GetTemperature(rho, primvar[PRS], *(rad_data->run_config));
  B = Blackbody(T, *(rad_data->run_config)) ;
  coeff1 = dt * abs_op * rho * B ;
  coeff2 = gamma * mdt_rho ;

  // -- Set coefficients of the system C.(E,F^i)=b  --
  b[0] = rad_data-> Rini[0] + coeff1 * gamma ;
  C[0][0] = 1.0 - dt * rhogamma * ( scat_op*( u2 + uuD ) - abs_op );

	for (i=0; i<comps; i++){
    b[i+1] = rad_data-> Rini[i+1] + coeff1 * u[i] ;

    C[0][i+1] = mdt_rho * u[i] * ( abs_op - scat_op * (u2 + gamma2) ) ;
    C[i+1][0] = mdt_rho * ( tot_op * uD[i] + scat_op * u[i] * (gamma2 + uuD) ) ;

    for (j=i; j<comps; j++)
      C[i+1][j+1] = C[j+1][i+1] = - coeff2 * scat_op * 2.0 * u[i] * u[j] ;

    C[i+1][i+1] += 1.0 - coeff2 * tot_op ;
  }

  #endif
}
#pragma omp end declare target
