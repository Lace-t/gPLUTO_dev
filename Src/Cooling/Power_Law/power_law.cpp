////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief  Take a source step using power-law cooling.

  Integrate the ODE 
  \f[
       dp_{\rm cgs}/dt_{\rm cgs}
      = -(\Gamma - 1) \Lambda(\rho_{\rm cgs}, T_{\rm cgs})
       \qquad{\rm where}\qquad
       \Lambda = \frac{a_{br}}{(\mu m_H)^2} \rho_{\rm cgs}^2 \sqrt{T_{\rm cgs}}
      \,,\qquad
         a_{br} = 2.e-27   c.g.s
  \f]
  which accounts for bremmstrahlung cooling.
  Here the subscript 'cgs' means that the corresponding quantity is given
  in c.g.s units. 
  We denote with \c mu the molecular weight, \c mH the hydrogen mass (in c.g.s).
  
  The previous integration is carried out analytically since the density
  does not change during this step.
  In non-dimensional form:
  \f[
      \frac{dp}{dt} = -{\rm cost} \;\rho  (\rho p)^{\HALF}
  \f]
 
   [notice that since  \c p/rho=T/KELVIN this is equivalent to:
    \c dp/dt=-cost rho^2 (T/KELVIN)^(1/2) ]
 
  The quantity \c cost is determined by transforming the dimensional
  equation into the non-dimensional one.
  If p, rho and t are in code (non-dimensional) units and if \c L_0,
  \c rho_0, and \c V_0 are the unit length, density and velocity,
  then \c cost is found to be:
  \verbatim
                   a_br * (gamma - 1) * L_0 * rho_0
        cost = -------------------------------------------
                sqrt(kB * mu * mH) * mu * mH * V_0^2
  \endverbatim
   where a_{br} = 2.e-27 (in c.g.s), kB is the Boltmann constant
  
 
  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Nov 18, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#define KELVIN 5.4e12
////////////////////////////////////////////////////////////////////////
void PowerLawCooling (Data &data, double dt, timeStep &Dts, Grid &grid)
/*!
  \param [in,out]  Vc_temp    a pointer to the PLUTO 3D data array containing
  pimitive variables.
  \param [in]      dt    the current integration time step
  \param [in]      Dts   a pointer to the timeStep structure
  \param [in]      grid  pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
  
  int     i, j, k;
  double  cost, dE;
  double  rho, p, T, p_f, T_f;

  double g_minCoolingTemp = 50.0; //this should put inside run_config
  double g_maxCoolingRate = 0.1;  //this should put inside run_config

  const Ary4D& Vc_temp      = data.Vc; 
  double g_gamma = data.run_config->gamma;

  cost  = UNIT_LENGTH*UNIT_DENSITY/(UNIT_VELOCITY*UNIT_VELOCITY);
  cost *= 2.e-27*(g_gamma-1.0)/ (0.5*CONST_mH*sqrt(0.5*CONST_mH*CONST_kB));
  
  // -------------------------------------------------------------
  // Integrate analytically
  // -------------------------------------------------------------

  dE = 1.e-18;

  RBox box;
  box.ibeg = grid.gbeg[0];
  box.jbeg = grid.gbeg[1];
  box.kbeg = grid.gbeg[2];

  box.iend = grid.gend[0];
  box.jend = grid.gend[1];
  box.kend = grid.gend[2];

  BOX_LOOP(box,k,j,i){

//  ----  Find initial temperature in Kelvin  ----  

    rho = Vc_temp(i,j,k,RHO);
    p   = Vc_temp(i,j,k,PRS);

    T   = (p/rho*Kelvin(*data.run_config));

    if (T < g_minCoolingTemp) continue;

//  ----  Find final energy  ----  

    p_f = sqrt(p) - 0.5*cost*rho*sqrt(rho)*dt;
    p_f = MAX(p_f, 0.0);
    p_f = p_f*p_f;
    T_f = p_f/rho*Kelvin(*data.run_config);
    T_f = MAX (T_f, g_minCoolingTemp);

//  ----  Update Energy  ----  

    p_f = T_f*rho/Kelvin(*data.run_config);

    Vc_temp(i,j,k,PRS) = p_f;
    dE = fabs(1.0 - p_f/p) + 1.e-18;

    Dts.dt_cool = MIN(Dts.dt_cool, dt*g_maxCoolingRate/dE);
  }

  
}

