////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Define the scalar resistivity for the RMHD equations.

  \authors G. Mattia \n
           A. Mignone (andrea.mignone@unito.it)
  \date    Dec 01, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

 ////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
double Resistive_eta(const State& state, const RunConfig& run_config)
/*!
  \param [in]  state      state structure
  \param [in]  run_config RunConfig Structure
*/
////////////////////////////////////////////////////////////////////////
{
  return 1.0;
}
#pragma omp end declare target
