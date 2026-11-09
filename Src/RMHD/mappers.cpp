////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Conservative to primitive conversion for RMHD

  Convert conservative variables to primitive variables for the RMHD
  module.

  \authors A. Mignone  (andrea.mignone@unito.it)
           M. Rossazza (marco.rossazza@unito.it)
           S. Truzzi   ()
  \date    Aug 21, 2025
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
  double   scrh;
  double  *u = state.u;
  double  *v = state.v;

  #if ENTROPY_SWITCH
  uint32_t *flag = state.flag;
  if (*flag & FLAG_ENTROPY) method = ENTROPY_SOLVE;
  #endif

  // ------------------------------------------------------
  // 1. Convert mag. field (trivial)
  // ------------------------------------------------------

  v[BX1] = u[BX1];
  v[BX2] = u[BX2];
  v[BX3] = u[BX3];

  // ------------------------------------------------------
  // 2. Check density and energy positivity
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
  // 3. Pressure recovery
  // ------------------------------------------------------

  // 3a. Try to invert using energy equation.

  if (method == ENERGY_SOLVE){
    err = RMHD_EnergySolve(u, v, run_config);
    #if PCP_NR == YES
    // Match PLUTO: Cai recovery is attempted after native recovery fails.
    if(err && !RMHD_PCP_NR(u,v,run_config)) err=0;
    #endif
    if (err){   // -- Something's wrong ? --> Try pressure fix -- 
      method  = PRESSURE_FIX;
    }
    err_mapper |= err;
  }

  // 3b. Try to invert using entropy equation.

  #if ENTROPY_SWITCH
  if (method == ENTROPY_SOLVE) {
    err = RMHD_EntropySolve(u, v, run_config);
    if (err) {   // -- Something's wrong ? --> Try pressure fix -- 
      method = PRESSURE_FIX;
    }
    err_mapper |= err;
  }
  #endif

  // 3c. Apply pressure fix

  if (method == PRESSURE_FIX){
    err  = RMHD_PressureFix(u, v, run_config);
    if (err != 0) {
      err_mapper |= err;
      err_mapper |= ERR_CONS2PRIM_FAIL;
    }
  }

  // 3d. Convert radiation fields
  #if RADIATION
  v[ENR] = u[ENR];
  v[FR1] = u[FR1]; 
  v[FR2] = u[FR2]; 
  v[FR3] = u[FR3];
  #endif

  // ------------------------------------------------------
  // 4. Check velocity is subluminal (this check
  //    may not be needed actually...)
  // ------------------------------------------------------

  scrh = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  if (scrh >= 1.0) {
    err_mapper |= ERR_SUPERLUMINAL_SPEED;
    err_mapper |= ERR_CONS2PRIM_FAIL;
  }

  #if RADIATION
  scrh = v[FR1]*v[FR1] + v[FR2]*v[FR2] + v[FR3]*v[FR3];
  scrh /= (v[ENR]*v[ENR]);
  if (scrh >= 1.0) {
    err_mapper |= ERR_SUPERLUMINAL_RADIATION_FLUX;
    err_mapper |= ERR_CONS2PRIM_FAIL;
  }
  #endif
  
  // ------------------------------------------------------
  // 5. Complete conversion, update flag
  // ------------------------------------------------------

  #if NSCL > 0
  int nv;
  NSCL_LOOP(nv) v[nv] = u[nv]/u[RHO];
  #endif

  #ifdef GLM_MHD
  v[PSI_GLM] = u[PSI_GLM];
  #endif
   
  return err_mapper;
}

#undef  ENERGY_SOLVE
#undef  ENTROPY_SOLVE
#undef  PRESSURE_FIX
