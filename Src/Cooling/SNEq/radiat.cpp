///////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Right hand side for SNEq cooling.

  \authors A. Mignone (andrea.mignone@unito.it)\n
           O. Tesileanu

  \b References
     - "Simulating radiative astrophysical flows with the PLUTO code:
        a non-equilibrium, multi-species cooling function" \n
       Tesileanu, Mignone and Massaglia, A&A (2008) 488, 429

  \date   Nov 23, 2025
*/
// ///////////////////////////////////////////////////////////////////// //
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void Radiat (double *v, double *rhs, RunConfig& run_conf)
/*!
  Cooling for neutral or singly ionized gas: good up to about 35,000 K
  in equilibrium or shocks in neutral gas up to about 80 km/s.
  Assumed abundances in ab
  Uses t : Kelvin
  dene :  electron density cm*-3
  fneut : hydrogen neutral fraction    (adimensionale)
  ci,cr : H ionization and recombination rate coefficients 
  
  
  em(1) = TOTAL EMISSIVITY : (ergs cm**3 s**-1)
  em(2) = Ly alpha  + two photon continuum: Aggarwal MNRAS 202,
  10**4.3 K
  em(3) = H alpha: Aggarwal, Case B
  em(4) = He I 584 + two photon + 623 (all n=2 excitations): Berrington
  &Kingston,JPB 20:
  
  em(5) = C I 9850 + 9823: Mendoza, IAU 103, 5000 K
  em(6) = C II, 156 micron: Mendoza, 10,000 K
  em(7) = C II] 2325 A: Mendoza, 15,000 K
  em(8) = N I 5200 A: Mendoza, 7500 K
  em(9) = N II 6584 + 6548 A: Mendoza
  em(10) = O I 63 micron: Mendoza,2500 K
  em(11) = O I 6300 A + 6363 A: Mendoza, 7500 K
  em(12) = O II 3727: Mendoza
  em(13) = Mg II 2800: Mendoza
  em(14) = Si II 35 micron: Dufton&Kingston, MNRAS 248
  em(15) = S II 6717+6727: Mendoza
  em(16) = Fe II 25 micron: Nussbaumer&Storey
  em(17) = Fe II 1.6 micron
  em(18) = thermal energy lost by ionization
  em(19) = thermal energy lost by recombination (2/3 kT per
  recombination.
  The ionization energy lost is not included here.
*/
////////////////////////////////////////////////////////////////////////
{
  int     kk;
  double  T, mu, st, rho, prs, fn;
  double  N_H, n_el, cr, ci, rlosst, em[20];

   double t00[18] = {0.0   , 0.0   , 1.18e5, 1.40e5, 2.46e5, 1.46e4, 
                           92.1  , 6.18e4, 2.76e4, 2.19e4, 228.0 , 2.28e4,
                           3.86e4, 5.13e4, 410.0 , 2.13e4, 575.0 , 8980.0};

   double  ep[18] = {0.0   , 0.0 , 10.2  , 1.89, 21.2  , 1.26, 
                           0.0079, 5.33, 2.38  , 1.89, 0.0197, 1.96,
                           3.33  , 4.43, 0.0354, 1.85, 0.0495, 0.775};

   double critn[18] = {0.0  , 0.0   , 1.e10, 1.e10, 1.e10,  312.0, 
                             0.849, 1.93e7, 124.0, 865.0, 1090.0, 3950.0, 
                             177.0, 1.e10 ,  16.8,  96.0,  580.0, 1130.0};

   double om[18] = {0.0 , 0.0 , 0.90, 0.35, 0.15  , 0.067,
                          0.63, 0.52, 0.90, 0.30, 0.0055, 0.19,
                          0.33, 8.0 , 2.85, 1.75, 0.3   ,0.39};

   double ab[18] = {0.0   , 0.0    , 1.0    , 1.0    , 0.1    , 0.0003,
                          0.0003, 0.0003 , 0.0001 , 0.0001 , 0.0006 , 0.0006,
                          0.0006, 0.00002, 0.00004, 0.00004, 0.00004, 0.00004};

   double fn1[20] = {0.0, 0.0, 0.0, 0.0, 0.0,
                           0.1, 1.0, 1.0, 0.0, 1.0, 
                           0.0, 0.0, 1.0, 1.0, 1.0, 
                           1.0, 1.0, 1.0, 0.0, 1.0};

   double fn2[20] = {0.0, 0.0, 1.0, 1.0, 1.0,
                           0.0, 0.0, 0.0, 1.0,-1.0,
                           1.0, 1.0,-1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0,-1.0};
   //int first_call = 1;
   double E_cost, unit_time, N_H_rho;
  const double unit_length   = run_conf.unit_length;
  const double unit_velocity = run_conf.unit_velocity;
  const double unit_density  = run_conf.unit_density;

  // 0. Compute constants & conversion factor from
  //    total density to hydrogen number density, i.e.
  //    nH = N_H_rho * rho        

  //if (first_call) {
    E_cost     = unit_length/unit_density/POW3(unit_velocity);
    unit_time  = unit_length/unit_velocity;
    N_H_rho    = (unit_density/CONST_amu)*(H_MASS_FRAC/CONST_AH);
   // first_call = 0;
  //}

  // 1. Force fneut to stay between [0,1] and
  //    compute temperature from rho, rhoe and fractions.
  
  v[X_HI] = MAX(v[X_HI], 0.0);
  v[X_HI] = MIN(v[X_HI], 1.0);

  rho = v[RHO];
  mu  = MeanMolecularWeight(v); 
  if (mu < 0.0){
    //QuitPLUTO("from Radiat(): Negative mu in radiat",1);
  }
  #if EOS == IDEAL
   if (v[RHOE] < 0.0) v[RHOE] = run_conf.small_prs/(run_conf.gamma-1.0);
   prs = v[RHOE]*(run_conf.gamma-1.0);
   T   = prs/rho*Kelvin(run_conf)*mu; 
  #else
   int status = GetEV_Temperature(v[RHOE], v, &T);
   if (status != 0) {
     T = T_CUT_RHOE;
     v[RHOE] = InternalEnergy(v, T);
   }
  #endif

  fn = v[X_HI];

  // ------------------------------------------------------
  // 2. Compute right hand sides for internal energy &
  //    ionization fraction.
  // ------------------------------------------------------

  st    = sqrt(T);
  N_H   = N_H_rho*rho;  // Number density of hydrogen N_H = N(HI) + N(HII)

  // Ionization & Recomb. coefficients for Hydrogen 

  cr = 2.6e-11/st;
  ci = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6);

  n_el = N_H*(1.0 -  v[X_HI] + FRAC_Z);  // electron number density, in cm^{-3} 
  rhs[X_HI] = unit_time*n_el*(-(ci + cr)*v[X_HI] + cr); 
  em[1] = 0.0;
  for (kk = 2; kk <= 17; kk++){
    em[kk] = 1.6e-12*8.63e-6*om[kk]*ep[kk]*exp(-t00[kk]/T)/st;
    em[kk] = em[kk]*critn[kk]*st/(n_el + critn[kk]*st);
    em[kk] = em[kk]*ab[kk]*(fn1[kk] +  v[X_HI]*fn2[kk]);
    em[1] = em[1] + em[kk];
  }

  em[18] = ci*13.6*1.6e-12* v[X_HI];
  em[19] = cr*0.67*1.6e-12*(1.0 -  v[X_HI])*T/11590.0;
  em[1]  = em[1] + em[18] + em[19];

  // 2a. rlosst is the energy loss in units of erg/cm^3/s;
  //    it must be multiplied by cost_E in order to match 
  //    non-dimensional units.
  //    Source term for the neutral fraction scales with unit_time


  // 2b. Suppress rhs when T falls below run_conf.minCoolingTemp

  //if (fn < 1.0e-10) fn = 0.0;
  fn = (fn < 1.0e-10) ? 0.0 : fn;
  fn = 1.0 - 1.0/cosh( pow( T/run_conf.minCoolingTemp, 12));
  rlosst    =  em[1]*n_el*N_H;
  rhs[RHOE] = -E_cost*rlosst*fn;
}

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
double CompEquil (double N, double T, double *v)
/*!
  Compute the equilibrium ionization balance for (rho,T)   
*/
////////////////////////////////////////////////////////////////////////
{
  double cr, ci, st, n_e;

  st = sqrt(T);
  cr = 2.6e-11/st;             // recombination / ionization coefficients  
  ci = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6);
  v[X_HI] = cr / (ci + cr);   // compute fraction of neutrals at equilibrium 
  n_e = (1.0 - v[X_HI] + FRAC_Z)*N;

  return n_e;
}
#pragma omp end declare target