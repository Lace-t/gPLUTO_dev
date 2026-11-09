////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Return the (default) number of ghost zones

  Return the number of ghost zones employed by the selected numerical
  algorithm.
  The minimum number for a 2nd-order algorithm is 2.
  Higher-order interpolation scheme may require more zones.

  User may override the number of ghost zones by defining the 
  macro \c NGHOST_USR

  \authors A. Mignone (andrea.mignone@unito.it)

  \date    Apr 17, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Return the number of ghost zones, depending on the selected
  scheme.
*/
////////////////////////////////////////////////////////////////////////
int GetNghost (void)
{
  #if NGHOST_USR > 0
  return NGHOST_USR;
  #else
  
  int nghost = 2;

  // --------------------------------------------------------
  // Choose stencil based on reconstruction.
  // --------------------------------------------------------

  #if    (RECONSTRUCTION == MP5)    || (RECONSTRUCTION == PARABOLIC)  \
      || (RECONSTRUCTION == WENOZ) \
      || (RECONSTRUCTION == MP5_FD) || (RECONSTRUCTION == WENOZ_FD)  \
      || (RECONSTRUCTION == LINEAR && LIMITER == FOURTH_ORDER_LIM) \
      || (RING_AVERAGE_REC > 2)
  nghost = 3;
  #endif  

  #if SHOCK_FLATTENING == MULTID

  // --------------------------------------------------------
  // The MULTID shock flattening only need 2 ghost zones.
  // However for axisymmetric simulations with CTU
  // 3 zones will ensure that flag[][][] will remain
  // symmetric around the axis.
  // --------------------------------------------------------

   nghost = MAX(3, nghost);
  #endif

  // --------------------------------------------------------
  // FARGO PPM needs at least 3 ghost zones
  // --------------------------------------------------------

  #ifdef FARGO
   #if FARGO_ORDER == 3
    nghost = MAX(3,nghost);
   #endif
  #endif

  // --------------------------------------------------------
  // Particles are required to propagate for no more than 2
  // zones. However, with option #2 of the CR predictor,
  // a single particle is advanced for dt*v (instead of
  // dt*v/2) thus potentially ending up in the very last
  // ghost zones. In this case, the FORCE cannot be properly
  // computed, unless an extra zone is added.
  // --------------------------------------------------------

  #if (PARTICLES == PARTICLES_CR)
  #if PARTICLES_CR_PREDICTOR == 2
  nghost = MAX(3, nghost);
  #endif
  #endif

    // --------------------------------------------------------
    // The high-order module requires minimum 4 ghost zones 
    // (3 reconstruction + 1 laplacians) when unlimited and 
    // with pointwise reconstructions.
    // All the other cases require a 5 point stencil.
    // --------------------------------------------------------

  #ifdef HIGH_ORDER
  nghost = 3;
  #endif

  return nghost;
  
#endif // NGHOST_USR > 0 
}
