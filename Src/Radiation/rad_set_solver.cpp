#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Depending on the choice of the radiation Riemann solver specified in
  pluto.ini, return a pointer to the corresponding rad. Riemann solver
  function

  \param[in] solver  Radiation Riemann solver name.
  \return Pointer to the selected radiation Riemann solver.
*/
////////////////////////////////////////////////////////////////////////
RiemannSolver *RadSetSolver (const char *solver)
{
  
  if (!strcmp(solver, "tvdlf"))         return (&Rad_LF_Solver_T);
  else if (!strcmp(solver, "hll"))      return (&Rad_HLL_Solver_T);
  else if (!strcmp(solver, "hllc"))     return (&Rad_HLLC_Solver_T);

  string msg = string("! RadSetSolver(): solver '") + solver + "' not available";
  QuitPLUTO(msg, 1);
  return NULL;

}
