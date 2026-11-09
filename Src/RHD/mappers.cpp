////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to
  an array of conservative variables for the RHD equations.

  The ConsToPrim() converts an array of conservative quantities to
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:

      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Dec 06, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define  ENERGY_SOLVE   1
#define  ENTROPY_SOLVE  2
#define  PRESSURE_FIX   3

////////////////////////////////////////////////////////////////////////
/*!
  Convert from a local state from conservative to primitive variables.
  
  \param [in,out]  state       reference to State structure
  \param [in]      run_config  reference to RunConfig structure
  
  \return Return 0 if conversion was successful.
  Otherwise an err code is returned if conversion did not
  go through or fixes were applied.
*/
////////////////////////////////////////////////////////////////////////
uint32_t ConsToPrim (State &state, RunConfig &run_config)
{
  int      method = ENERGY_SOLVE;
  uint32_t err_mapper = 0, err  = 0;
  double  *u = state.u;
  double  *v = state.v;

  #if ENTROPY_SWITCH
  uint32_t *flag = state.flag;
  if (*flag & FLAG_ENTROPY) method = ENTROPY_SOLVE;
  #endif

  // ------------------------------------------------------
  // 1. Check density and energy positivity
  // ------------------------------------------------------

  if (u[RHO] < 0.0) {
    u[RHO]      = run_config.small_rho;
    err_mapper |= ERR_NEGATIVE_RHO;
  }

  if (u[ENG] < 0.0) {
    method      = PRESSURE_FIX;
    err_mapper |= ERR_NEGATIVE_ENG;
  }

  #if RADIATION
  if (u[ENR] < 0.0) {
    u[ENR]      = run_config.small_enr;
    err_mapper |= ERR_NEGATIVE_ENR;
  }
  #endif

  // ------------------------------------------------------
  // 2a. Try to invert using energy eq.
  // ------------------------------------------------------

  if (method == ENERGY_SOLVE){
    err = RHD_EnergySolve(u, v, run_config);
    if (err){   // -- Didn't work out ? --> Try pressure fix --
      method  = PRESSURE_FIX;
    }
    err_mapper |= err;
  }

  // ------------------------------------------------------
  // 2b. Try to invert using entropy equation.
  // ------------------------------------------------------

  #if ENTROPY_SWITCH
  if (method == ENTROPY_SOLVE) {
    err = RHD_EntropySolve(u, v, run_config);
    if (err) {   // -- Something's wrong ? --> Try pressure fix --
      method = PRESSURE_FIX;
    }
    err_mapper |= err;
  }
  #endif

  // ------------------------------------------------------
  // 2c. Apply pressure fix if necessary
  // ------------------------------------------------------

  if (method == PRESSURE_FIX){
    err = RHD_PressureFix(u, v, run_config);
    if (err != 0) {
      err_mapper |= err;
      err_mapper |= ERR_CONS2PRIM_FAIL;
    }
  }

  // ------------------------------------------------------
  // 3. Complete conversion
  // ------------------------------------------------------

  #if RADIATION
  v[ENR] = u[ENR];
  v[FR1] = u[FR1]; 
  v[FR2] = u[FR2]; 
  v[FR3] = u[FR3];
  #endif

  #if NSCL > 0
  int nv;
  NSCL_LOOP(nv)  v[nv] = u[nv]/u[RHO];
  #endif

  return err_mapper;
}

#undef  ENERGY_SOLVE
#undef  ENTROPY_SOLVE
#undef  PRESSURE_FIX
