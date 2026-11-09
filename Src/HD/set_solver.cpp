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

  #if EOS == IDEAL

    if      (!strcmp(solver, "tvdlf")) return (&LF_Solver_T);
    else if (!strcmp(solver, "roe"))   return (&Roe_Solver_T);
    else if (!strcmp(solver, "hll"))   return (&HLL_Solver_T);
    else if (!strcmp(solver, "hllc"))  return (&HLLC_Solver_T);
    else if (!strcmp(solver, "two_shock")) return (&TwoShock_Solver_T);
  #elif EOS == ISOTHERMAL
    if      (!strcmp(solver, "tvdlf")) return (&LF_Solver_T);
    else if (!strcmp(solver, "roe"))   return (&Roe_Solver_T);
    else if (!strcmp(solver, "hll"))   return (&HLL_Solver_T);
    else if (!strcmp(solver, "hllc"))  return (&HLLC_Solver_T);
  #endif

  string msg = string("! SetSolver(): solver '") + solver + "' not available";
  QuitPLUTO(msg, 1);
  return NULL;
}
