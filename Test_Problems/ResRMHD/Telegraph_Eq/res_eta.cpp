////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Define the scalar resistivity for the RMHD equations. 

  \authors G. Mattia \n
           A. Mignone (mignone@ph.unito.it)
  \date    June 09, 2017
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double Resistive_eta(const State& state, const RunConfig& run_config)
/*!
  \param [in]  state      state structure
  \param [in]  run_config RunConfig Structure
*/
////////////////////////////////////////////////////////////////////////
{
  const double eta = run_config.inputParam[ETA];
  return eta;
}


