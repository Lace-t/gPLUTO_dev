////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Conservative to primitive conversion for ResRMHD

  Convert conservative variables to primitive variables for the ResRMHD
  module.

  \authors A. Mignone  (andrea.mignone@unito.it)
           M. Rossazza (marco.rossazza@unito.it)
           S. Truzzi   ()
  \date    Jan 3, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define  ENERGY_SOLVE   1
#define  ENTROPY_SOLVE  2
#define  PRESSURE_FIX   3

////////////////////////////////////////////////////////////////////////
uint32_t ConsToPrim (State &state, RunConfig &run_config)
/*!
  Convert from a local state from conservative to primitive variables.
  
  \param [in,out]  state       reference to State structure
  \param [in]      run_config  reference to RunConfig structure
  
  \return Return 0 if conversion was successful.
  Otherwise an err code is returned if conversion did not
  go through or fixes were applied.
                                                                      */
////////////////////////////////////////////////////////////////////////
{
  int      method = ENERGY_SOLVE;
  uint32_t err_mapper = 0;

  double Emag2, Bmag2;
  double uh[NVAR];

  double  *u = state.u;
  double  *v = state.v;

  #if ENTROPY_SWITCH
  uint32_t *flag = state.flag;
  if (*flag & FLAG_ENTROPY) method = ENTROPY_SOLVE;
  #endif

  // ------------------------------------------------------
  // 1. Trivial stuffs go first...
  // ------------------------------------------------------

  v[BX1] = u[BX1];
  v[BX2] = u[BX2];
  v[BX3] = u[BX3];

  v[EX1] = u[EX1];
  v[EX2] = u[EX2];
  v[EX3] = u[EX3];

  #ifdef CRG
  v[CRG] = u[CRG];
  #endif

  Emag2 = u[EX1]*u[EX1] + u[EX2]*u[EX2] + u[EX3]*u[EX3];
  Bmag2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];

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

  // ------------------------------------------------------
  // 3. Build "hydro" conserved variables by subtracting
  //    EM contribution.
  // ------------------------------------------------------

  uh[RHO] = u[RHO];
  uh[MX1] = u[MX1] - (u[EX2]*u[BX3] - u[EX3]*u[BX2]);
  uh[MX2] = u[MX2] - (u[EX3]*u[BX1] - u[EX1]*u[BX3]);
  uh[MX3] = u[MX3] - (u[EX1]*u[BX2] - u[EX2]*u[BX1]);
  uh[ENG] = u[ENG] - 0.5*(Emag2 + Bmag2);
  #if ENTROPY_SWITCH
  uh[ENTR] = u[ENTR];
  #endif

  // 4a. Attempt to recover primitive variables
  //     by inverting the energy equation (use RHD solvers).

  if (method == ENERGY_SOLVE){
    int err = RHD_EnergySolve(uh, v, run_config);
    if (err){   // Didn't go through ? --> Try entropy solve or pressure fix
      method      = PRESSURE_FIX;
      err_mapper |= err;
    }
    #if ENTROPY_SWITCH   // Entropy has been re-defined in RHD_EnergySolve()
    u[ENTR] = uh[ENTR];  // and should be copied back to u[]
    #endif
  }

  // 4a. Attempt to recover primitive variables
  //     by inverting the entropy equation (use RHD solvers).

  #if ENTROPY_SWITCH
  if (method == ENTROPY_SOLVE) {
    int err = RHD_EntropySolve(uh, v, run_config);
    if (err){   // Didn't go through ? --> Try pressure fix 
      method      = PRESSURE_FIX;
      err_mapper |= err;
    }
    u[ENG] = uh[ENG] + 0.5*(Emag2 + Bmag2); // uh[ENG] updated in RHD_EntropySolve()
  }
  #endif

  // 4b. Apply pressure fix if necessary

  if (method == PRESSURE_FIX){
    int err = RHD_PressureFix(uh, v, run_config);
    if (err) {
      err_mapper |= err;
      err_mapper |= ERR_CONS2PRIM_FAIL;
    }
    u[ENG] = uh[ENG] + 0.5*(Emag2 + Bmag2);
    #if ENTROPY_SWITCH 
    u[ENTR] = uh[ENTR];
    #endif
  } // end if pressure_fix
  
  // ------------------------------------------------------
  // 5. Complete conversion
  // ------------------------------------------------------

  #if NSCL > 0
  int nv;
  NSCL_LOOP(nv) v[nv] = u[nv]/u[RHO];
  #endif

  #ifdef PSI_GLM
  v[PSI_GLM] = u[PSI_GLM];
  #endif

  #ifdef PHI_GLM
  v[PHI_GLM] = u[PHI_GLM];
  #endif

  return err_mapper;
}

#undef  ENERGY_SOLVE
#undef  ENTROPY_SOLVE
#undef  PRESSURE_FIX
