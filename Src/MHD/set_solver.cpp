////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Return a pointer to a Riemann solver function.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Mar 05, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
RiemannSolver *SetSolver (const char *solver)
/*!
  Depending on the choice of the Riemann solver specified in
  pluto.ini, return a pointer to the corresponding Riemann solver
  function
*/
////////////////////////////////////////////////////////////////////////
{

  // 0. Set Pointers to different Riemann solvers

#if HALL_MHD
  if (strcmp(solver,"hll")) {
    printLog ("! SetSolver(): %s not compatible with HALL_MHD, use hll instead\n",
            solver);
    QuitPLUTO (" ", 1);
  }
#endif

  // Set pointers to Riemann solvers

  #if EOS == IDEAL
    if (!strcmp(solver, "tvdlf"))        return (&LF_Solver_T);
    else if (!strcmp(solver, "gforce"))  return (&GFORCE_Solver_T);
    else if (!strcmp(solver, "roe"))     return (&Roe_Solver_T);
    else if (!strcmp(solver, "hllc"))    return (&HLLC_Solver_T);
    else if (!strcmp(solver, "hlld"))    return (&HLLD_Solver_T);
    else if (!strcmp(solver, "hll"))     return (&HLL_Solver_T);
  #elif EOS == ISOTHERMAL
    if (!strcmp(solver, "tvdlf"))        return (&LF_Solver_T);
    else if (!strcmp(solver, "gforce"))  return (&GFORCE_Solver_T);
    else if (!strcmp(solver, "roe"))     return (&Roe_Solver_T);
    else if (!strcmp(solver, "hllc"))    return (&HLLC_Solver_T);
    else if (!strcmp(solver, "hlld"))    return (&HLLD_Solver_T);
    else if (!strcmp(solver, "hll"))     return (&HLL_Solver_T);
  #endif

  string msg = string("! SetSolver(): solver '") + solver + "' not available";
  QuitPLUTO(msg, 1);

  return NULL;
}
