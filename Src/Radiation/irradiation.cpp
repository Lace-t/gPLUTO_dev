////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Template for the irradiation flux divergence.

  Defines the divergence of the irradiation flux field stored in
  v[FIR] (= u[FIR]). This value is treated as a passive scalar
  advected with the fluid; its spatial distribution is set here
  and optionally updated every time step (IRRADIATION_UPDATE == YES).

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 8 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if IRRADIATION
////////////////////////////////////////////////////////////////////////
/*!
  Compute the divergence of the irradiation flux field and store it
  in \c data->Vc(i,j,k,FIR) (primitive) and \c data->Uc(i,j,k,FIR)
  (conservative, identical).

  \param[in,out] data  Pointer to Data structure.
  \param[in]     grid  Pointer to Grid structure.
*/
////////////////////////////////////////////////////////////////////////
void RadIrradiationFlux (Data *data, Grid *grid)
{
  int i, j, k;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    // Store the divergence of the irradiation flux (energy density / time).
    data->Vc(i,j,k,FIR) = 0.0;
    data->Uc(i,j,k,FIR) = data->Vc(i,j,k,FIR);
  }}}
}
#endif
