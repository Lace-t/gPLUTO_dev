////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Conservative to primitive conversion for RMHD

  Convert conservative variables to primitive variables for the RMHD
  module.

  \authors A. Mignone  (andrea.mignone@unito.it)
           M. Rossazza (marco.rossazza@unito.it)
  \date    Dec 07, 2025
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
  int method = ENERGY_SOLVE;
  uint32_t err_mapper = 0;
  double tau;
  double  *u = state.u;
  double  *v = state.v;

  #if EOS == IDEAL
  double gmm1 = run_config.gamma - 1.0;
  #endif

  #if ENTROPY_SWITCH
  uint32_t *flag = state.flag;
  if (*flag & FLAG_ENTROPY) method = ENTROPY_SOLVE;
  #endif

  // ------------------------------------------------------
  // 1. Check density positivity
  // ------------------------------------------------------

  if (u[RHO] < 0.0) {
    u[RHO]      = run_config.small_rho;
    err_mapper |= ERR_NEGATIVE_RHO;
    err_mapper |= ERR_CONS2PRIM_FAIL;
  }

  // ------------------------------------------------------
  // 2. Compute density, velocity and mag. field
  // ------------------------------------------------------

  v[RHO] = u[RHO];
  tau    = 1.0/u[RHO];
  v[VX1] = u[MX1]*tau;
  v[VX2] = u[MX2]*tau;
  v[VX3] = u[MX3]*tau;

  v[BX1] = u[BX1];
  v[BX2] = u[BX2];
  v[BX3] = u[BX3];

  // ------------------------------------------------------
  // 3. Compute pressure from total energy or entropy
  // ------------------------------------------------------

  // 3a. Recover pressure using total energy

#if EOS == IDEAL
  double m2, b2, kinb2;
  m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
  b2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];
  kinb2 = 0.5*(m2*tau + b2);
  if (method == ENERGY_SOLVE){
    if (u[ENG] > 0.0) {
      v[PRS] = gmm1*(u[ENG] - kinb2);
      if (v[PRS] < 0.0){
        v[PRS]      = run_config.small_prs;
        u[ENG]      = v[PRS]/gmm1 + kinb2; // Recompute energy
        err_mapper |= ERR_NEGATIVE_PRS;
      }
    }else {
      u[ENG]      = run_config.small_prs/gmm1 + kinb2;
      v[PRS]      = run_config.small_prs;
      err_mapper |= ERR_NEGATIVE_ENG;
    }
    #if ENTROPY_SWITCH
    double rho = v[RHO];
    u[ENTR] = v[PRS]/pow(rho,gmm1);  // Recompute entropy
    #endif
  }

  // 3b. Recover pressure using entropy

  #if ENTROPY_SWITCH
  if (method == ENTROPY_SOLVE){
    double rho    = v[RHO];
    double rhog1  = pow(rho, gmm1);
    v[PRS] = u[ENTR]*rhog1; 
    if (v[PRS] < 0.0 || u[ENTR] < 0.0){
      v[PRS]      = run_config.small_prs;
      u[ENTR]     = v[PRS]/rhog1;  
      err_mapper |= ERR_NEGATIVE_ENTR;
      err_mapper |= ERR_CONS2PRIM_FAIL;
    }
    u[ENG] = v[PRS]/gmm1 + kinb2; // recompute energy 
  }
  #endif // ENTROPY_SWITCH 

#endif //  EOS == IDEAL 

  // ------------------------------------------------------
  // 4. Complete conversion
  // ------------------------------------------------------

  #ifdef GLM_MHD
  v[PSI_GLM] = u[PSI_GLM];
  #endif

  #if (NSCL > 0) || RADIATION
  int nv;
  #endif

  #if NSCL > 0
  NSCL_LOOP(nv) v[nv] = u[nv]/u[RHO];
  #endif

  #if RADIATION
  if (u[ENR] < 0.0) {
    u[ENR]      = run_config.small_enr;
    err_mapper |= ERR_NEGATIVE_ENR;
  }
  NRAD_LOOP(nv) v[nv] = u[nv];
  #if IRRADIATION
  v[FIR] = u[FIR];
  #endif
  #endif

  return err_mapper;
}

#undef  ENERGY_SOLVE
#undef  ENTROPY_SOLVE
#undef  PRESSURE_FIX
