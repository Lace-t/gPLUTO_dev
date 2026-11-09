// /////////////////////////////////////////////////////////////////////
/*! 
 \file
 \brief Compute the flux for the M1 radiation transport equations.

  Compute the flux of the conservative radiation transport equations
  in the direction given by ::g_dir.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    Dec 12, 2019

*/
// /////////////////////////////////////////////////////////////////////

template<int dir> void RadFlux (State &, double *);

#pragma omp declare target
#pragma acc routine seq
template<int dir>
////////////////////////////////////////////////////////////////////////
/*!
  Compute the 1D flux for the radiation block.

  \tparam     dir   Sweep direction (IDIR, JDIR, or KDIR).
  \param[in]  state  Input radiation state.
  \param[out] flux   Flux vector.
*/
////////////////////////////////////////////////////////////////////////
inline void RadFlux (State &state, double *flux)
{
  constexpr int FRn = FR1 + dir;
  double *u = state.u;

  flux[ENR] = u[FRn] ;
  flux[FR1] = EddTensor (u, FR1, FRn) * u[ENR] ;
  flux[FR2] = EddTensor (u, FR2, FRn) * u[ENR] ;
  flux[FR3] = EddTensor (u, FR3, FRn) * u[ENR] ;
}
#pragma omp end declare target