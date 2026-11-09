#include"gpluto.hpp"

////////////////////////////////////////////////////////////////////////
RiemannSolver *SetSolver (const char *solver)
/*
 *   return a pointer to a riemann solver function
*/
////////////////////////////////////////////////////////////////////////
{

  // 0. Set Pointers to different Riemann solvers

   if (!strcmp(solver, "tvdlf"))     return (&LF_Solver_T);
   else if (!strcmp(solver, "hll"))  return (&HLL_Solver_T);
   else if (!strcmp(solver, "hllc")) return (&HLLC_Solver_T);
   else if (!strcmp(solver, "two_shock")) return (&TwoShock_Solver_T);

  string msg = string("! SetSolver(): solver '") + solver + "' not available";
  QuitPLUTO(msg, 1);
  return NULL;
}
