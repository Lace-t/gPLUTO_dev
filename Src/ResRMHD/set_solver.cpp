////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Return a pointer to a Riemann solver function.

  \author A. Mignone (mignone@to.infn.it)
  \date   Oct 16, 2020
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

   if      (!strcmp(solver, "tvdlf"))  return (&LF_Solver_T);
   else if (!strcmp(solver, "mhll"))   return (&MHLL_Solver_T);
   else if (!strcmp(solver, "mhllc"))  return (&MHLLC_Solver_T);

  string msg = string("! SetSolver(): solver '") + solver + "' not available";
  QuitPLUTO(msg, 1);
  return NULL;
}
