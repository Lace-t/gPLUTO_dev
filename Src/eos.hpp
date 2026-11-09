////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief gPLUTO header file for Equation of State definition.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   March 13, 2025
*/
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline double Enthalpy (double *v, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double h;

  // ----------------------------------------------
  // 1. Classical flows 
  // ----------------------------------------------

#if (PHYSICS == HD) || (PHYSICS == MHD)
  h = 1.0;
  return h;
#endif

  // ----------------------------------------------
  // 2. Relativistic Flows
  // ----------------------------------------------

#if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == GRMHD) || (PHYSICS == ResRMHD)
  
  double theta;

 // --> MODIFY AND COMPUTE rho*h instead of h

  #if EOS == IDEAL
  double gamma = run_config.gamma;
  double gmmr  = gamma/(gamma - 1.0);
  theta = v[PRS]/v[RHO];
  h     = 1.0 + gmmr*theta;
  #elif EOS == TAUB
  theta = v[PRS]/v[RHO];
  h     = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
  #endif

  return h;
#endif

}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline double Entropy (double *v, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double s;
  double rho = v[RHO];

  #if EOS == IDEAL  
  {
    double gamma = run_config.gamma;
    s = v[PRS]/pow(rho,gamma);  
  }
  #elif EOS == TAUB
  {
    double th = v[PRS]/rho;
    s  = v[PRS]/pow(rho,5./3.)*(1.5*th + sqrt(2.25*th*th + 1.0));  
  }
  #else
  s = 0.0; // not defined
  #endif

  return s;
}

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline void SoundSpeed2 (State &state, RunConfig &run_config)
/*!
  Return the square of the sound speed
*/
////////////////////////////////////////////////////////////////////////
{

  // ----------------------------------------------
  // 1. Classical flows 
  // ----------------------------------------------

#if (PHYSICS == HD) || (PHYSICS == MHD)

  #if EOS == IDEAL  
  state.a2 = run_config.gamma*state.v[PRS]/state.v[RHO];
  #elif EOS == ISOTHERMAL
  double a = run_config.isoSoundSpeed;
  state.a2 = a*a; 
  #else
  #error ! SoundSpeed2: Unknown EOS
  #endif

#endif

  // ----------------------------------------------
  // 2. Relativistic Flows
  // ----------------------------------------------

#if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == GRMHD) || (PHYSICS == ResRMHD)

  double h = Enthalpy(state.v, run_config);
  double theta = state.v[PRS]/state.v[RHO];
  state.h = h;

  #if EOS == IDEAL
  state.a2 = run_config.gamma*theta/h;
  #elif EOS == TAUB
  state.a2 = theta/(3.0*h)*(5.0*h - 8.0*theta)/(h - theta);
  #else
  #error ! SoundSpeed2: Unknown EOS
  #endif

#endif

}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline double Kelvin (RunConfig &run_config)
/*!
  Define the conversion constant between dimensionless temperature 
  prs/rho and physical temperature in Kelvin,  
  T = (prs/rho)*KELVIN*mu
*/
////////////////////////////////////////////////////////////////////////
{ 
  double uvel = run_config.unit_velocity;
  double Tk   = (uvel*uvel*CONST_amu/CONST_kB);
  return Tk;
}
