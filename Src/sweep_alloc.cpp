#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Allocate memory areas for arrays inside the sweep
  structure.
*/
////////////////////////////////////////////////////////////////////////
void AllocateSweepArrays (Sweep *sweep, Grid *grid)
{
  long int size_tot = grid->ntot[KDIR]*grid->ntot[JDIR]*grid->ntot[IDIR];
  long int size_tot_aug = size_tot;

#ifdef FARGO

  // 0. Arrays allocated using size_tot_aug are recycled in
  //   FARGO_ShiftSolution to save memory. In parallel mode, 
  //   one-dimensional arrays must be large enough to contain 
  //   data values coming from neighbour processors.

  #if GEOMETRY == SPHERICAL && DIMENSIONS == 3
  size_tot_aug = (grid->ntot[KDIR]+2*FARGO_BUF)*grid->ntot[JDIR]*grid->ntot[IDIR];
  #else
  size_tot_aug = grid->ntot[KDIR]*(grid->ntot[JDIR]+2*FARGO_BUF)*grid->ntot[IDIR];
  #endif

#endif

  // 1. Allocate memory for structure members

  sweep->uC      = ARRAY_1D(size_tot*NVAR, double);
  sweep->vL      = ARRAY_1D(size_tot_aug*NVAR, double);
  sweep->vR      = ARRAY_1D(size_tot_aug*NVAR, double);
  sweep->vC      = ARRAY_1D(size_tot_aug*NVAR, double);
  sweep->fL      = ARRAY_1D(size_tot_aug*NVAR, double);
  sweep->fR      = ARRAY_1D(size_tot_aug*NVAR, double);
  sweep->flux    = ARRAY_1D(size_tot_aug*NVAR, double);
  sweep->rhs     = ARRAY_1D(size_tot*NVAR, double);

  sweep->cs2     = ARRAY_1D(size_tot, double);
  sweep->press   = ARRAY_1D(size_tot, double);
  sweep->pL      = ARRAY_1D(size_tot, double);
  sweep->pR      = ARRAY_1D(size_tot, double);
  sweep->Bn      = ARRAY_1D(size_tot, double);
  sweep->En      = ARRAY_1D(size_tot, double);
  sweep->PSIn    = ARRAY_1D(size_tot, double);
  sweep->SL      = ARRAY_1D(size_tot, double);
  sweep->SR      = ARRAY_1D(size_tot, double);

  #if RADIATION
  sweep->SrL     = ARRAY_1D(grid->nmax, double);
  sweep->SrR     = ARRAY_1D(grid->nmax, double);
  #endif

  #if (SELF_GRAV)
  sweep->g_c     = ARRAY_1D(size_tot, double);;  
  #endif

  sweep->SaL      = ARRAY_1D(size_tot, double);
  sweep->SaR      = ARRAY_1D(size_tot, double);
  sweep->Sc       = ARRAY_1D(size_tot, double);
  sweep->dL       = ARRAY_1D(size_tot, double);
  sweep->dR       = ARRAY_1D(size_tot, double);
  sweep->aL       = ARRAY_1D(size_tot, double);
  sweep->aR       = ARRAY_1D(size_tot, double);

  sweep->flag    = ARRAY_1D(size_tot, uint32_t);
}

