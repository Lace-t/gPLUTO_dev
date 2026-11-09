////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Compute right hand side for Tabulated cooling

  \authors A. Mignone (andrea.mignone@unito.it)\n
           S. Truzzi

 \b References

  \date   Nov 22, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"


////////////////////////////////////////////////////////////////////////
void Radiat (double *v, double *rhs, RunConfig& run_conf)
/*!
  Provide r.h.s. for tabulated cooling.
*/
////////////////////////////////////////////////////////////////////////
{ 
  double ntab = run_conf.ntab;
  int    klo, khi, kmid;
  double  mu, T, Tmid, scrh, dT, prs;
  double E_cost;
  double unit_density  = run_conf.unit_density;
  double unit_length   = run_conf.unit_length;
  double unit_velocity = run_conf.unit_velocity;

  double nH = unit_density/CONST_amu*H_MASS_FRAC/CONST_AH*v[RHO];
  double ne = nH*(1.0 + 0.5*CONST_AZ*FRAC_Z);

  double small_prs= run_conf.small_prs;
  double  g_minCoolingTemp = run_conf.minCoolingTemp;
  double g_gamma = run_conf.gamma;

  // -------------------------------------------
  // Read tabulated cooling function
  // -------------------------------------------

  
  // ---------------------------------------------
  // Get pressure and temperature 
  // ---------------------------------------------

  E_cost = unit_length/unit_density/POW3(unit_velocity);
  prs = v[RHOE]*(g_gamma-1.0);

  if (prs < 0.0) {
    v[RHOE] = small_prs/(g_gamma - 1.0);
  }

  //mu  = MeanMolecularWeight(v);
  mu = (CONST_AH + FRAC_He*CONST_AHe + FRAC_Z*CONST_AZ) /
  (2.0 + 3.0*FRAC_He + FRAC_Z*(1.0 + CONST_AZ*0.5));
  
  T   = prs/v[RHO]*Kelvin(run_conf)*mu;

  if (isnan(T)){
    printf ("! Radiat(): Nan found: rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QuitPLUTO(" ",1);
  }

/*
  if (T < g_minCoolingTemp) { 
    rhs[RHOE] = 0.0;
    return;
  }
*/
  // ----------------------------------------------
  // Table lookup by binary search  
  // ----------------------------------------------

  klo = 0;
  khi = ntab - 1;

  if (T > run_conf.T_tab[khi] || T < run_conf.T_tab[klo]){
    rhs[RHOE] = 0.0;
    return;
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = run_conf.T_tab[kmid];

    T <= Tmid ? (khi =kmid) : (klo=kmid);
    
    /*
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
    */
  }

  // -----------------------------------------------
  // Compute r.h.s
  // -----------------------------------------------

  dT        = run_conf.T_tab[khi] - run_conf.T_tab[klo];
  scrh      = run_conf.L_tab[klo]*(run_conf.T_tab[khi] - T)/dT + run_conf.L_tab[khi]*(T - run_conf.T_tab[klo])/dT;
  rhs[RHOE] = -ne*nH*scrh*E_cost;
  
  // ----------------------------------------------
  // Temperature cutoff
  // ----------------------------------------------

  rhs[RHOE] *= 1.0 - 1.0/cosh( pow( T/g_minCoolingTemp, 12));

}
