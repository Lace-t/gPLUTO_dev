////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Include source terms using operator splitting.

  The SplitSource() function handles source terms in a separate
  step using operator splitting.
  It is called from Integrate() between standard hydro advances.
  At present these source terms are one or more of the following:

  - optically thin radiative losses (cooling)
  - Diffusion operators: 
    - resistivity 
    - Thermal conduction
    - Viscosity
  - additional user-defined terms may also be included here.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Apr 17, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void SplitSource (Data *data, double dt, timeStep *Dts, Grid *grid)
/*! 
  Take one step on operator-split source terms.
  
  \param [in,out] d   pointer to PLUTO Data structure containing 
  the solution array updated from the most 
  recent call
  \param[in]      dt  the time step used to integrate the source 
  terms
  \param[in]     Dts  pointer to the time step structre
  \param[in]    grid  pointer to an array of grid structures
*/
////////////////////////////////////////////////////////////////////////
{
  #if (!defined(GLM_MHD)) && (COOLING == NO) && !RADIATION_NR
  #endif
  
  // ------------------------------------------------------
  // 0. Radiation substepping (NR module only)
  // ------------------------------------------------------

  #if RADIATION_NR
  RadSubstepping(data, dt, Dts, grid);
  #endif

  // ------------------------------------------------------
  // 1. Include source terms due GLM
  // ------------------------------------------------------

  #ifdef GLM_MHD 
  GLM_Source (data, dt, grid);
  #endif

  // ------------------------------------------------------
  // 2. Include source terms due to Cooling/Heating losses
  // ------------------------------------------------------

  #if COOLING != NO
    #if COOLING == POWER_LAW  
    /** PowerLawCooling function for POWER LAW COOLING */
    PowerLawCooling (*data, dt, *Dts, *grid);
    #elif COOLING == SYNCH
    SynchrotronCooling (*data, dt, *Dts, *grid);

    #else 

    CoolingSource (data, dt, Dts, grid);

    #endif
  #else
  #if !RADIATION_NR
  (void)Dts;
  #endif
  #endif

  // ------------------------------------------------------
  // 3. Include parabolic terms via STS/RKL
  //  - resistivity 
  //  - thermal conduction
  //  - viscosity 
  // ------------------------------------------------------
/*
#if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
  STS (data, dt, Dts, grid);
#endif

#if (PARABOLIC_FLUX & RK_LEGENDRE)
  RKL (data, dt, Dts, grid);
#endif
*/
}