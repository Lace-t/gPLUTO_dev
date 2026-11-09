////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Cooling main header file.

  Contains function prototypes and basic definitions for the cooling
  modules.
  
  \author A. Mignone (andrea.mignone@unito.it)
  \date   Nov 30, 2025
*/
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
/* The number of variables used in the cooling modules
    must also include internal energy (RHOE) for which 
    cooling functions are best defined.                               */
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/* Define IONS_BEG to be the array index for the first ion.           */
////////////////////////////////////////////////////////////////////////

#define IONS_BEG   (NFLX + NTRACER + (ENTROPY_SWITCH != 0)) 

// Cooling time Stepping labels

#define RKF12        1 
#define RKF23        2
#define ROS34        3

// Tunable macros

#ifndef COOLING_ODE_TOL           // Accruacy used to estimate the error 
                                  // in the ode time stepping 
  #define COOLING_ODE_TOL 2.e-5    // Same as PLUTO
#endif

#ifndef COOLING_TIME_STEPPING            // Ode time-stepping
  #define COOLING_TIME_STEPPING  RKF23
#endif

////////////////////////////////////////////////////////////////////////
/* Cooling module specific definitions                                */
////////////////////////////////////////////////////////////////////////
#if  COOLING == TABULATED || COOLING == POWER_LAW || COOLING == SYNCH
  
  #define NIONS   0
  
#elif COOLING == SNEq

  #define NIONS  1
  #define X_HI   IONS_BEG

#endif 

#if 0 // not implemented yet

#if COOLING == H2_COOL

  void   H2RateTables(double, double *);

  #define NIONS    3
  #define X_HI     (IONS_BEG)
  #define X_H2     (IONS_BEG + 1)
  #define X_HII    (IONS_BEG + 2)

#elif COOLING == MINEq

// -------------------------------------------------------------
// Notice: the order is absolutely important and it
// MUST NOT be changed !!!
// -------------------------------------------------------------

  #define C_IONS  3   // in [1,5] 
  #define N_IONS  3   // in [1,5] 
  #define O_IONS  3   // in [1,5] 
  #define Ne_IONS 3   // in [1,5] 
  #define S_IONS  3   // in [1,5] 
  #define Fe_IONS 0   // in [0,3] 
  
  #if C_IONS == 0
   #define C_EXPAND(a,b,c,d,e)  
  #elif C_IONS == 1       
   #define C_EXPAND(a,b,c,d,e)  ,a
  #elif C_IONS == 2     
   #define C_EXPAND(a,b,c,d,e)  ,a,b
  #elif C_IONS == 3       
   #define C_EXPAND(a,b,c,d,e)  ,a,b,c 
  #elif C_IONS == 4
   #define C_EXPAND(a,b,c,d,e)  ,a,b,c,d
  #elif C_IONS == 5      
   #define C_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
  #endif
  
  #if N_IONS == 0
   #define N_EXPAND(a,b,c,d,e)  
  #elif N_IONS == 1      
   #define N_EXPAND(a,b,c,d,e)  ,a
  #elif N_IONS == 2     
   #define N_EXPAND(a,b,c,d,e)  ,a,b
  #elif N_IONS == 3       
   #define N_EXPAND(a,b,c,d,e)  ,a,b,c 
  #elif N_IONS == 4
   #define N_EXPAND(a,b,c,d,e)  ,a,b,c,d
  #elif N_IONS == 5      
   #define N_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
  #endif
  
  #if O_IONS == 0
   #define O_EXPAND(a,b,c,d,e)  
  #elif O_IONS == 1      
   #define O_EXPAND(a,b,c,d,e)  ,a
  #elif O_IONS == 2     
   #define O_EXPAND(a,b,c,d,e)  ,a,b
  #elif O_IONS == 3       
   #define O_EXPAND(a,b,c,d,e)  ,a,b,c 
  #elif O_IONS == 4
   #define O_EXPAND(a,b,c,d,e)  ,a,b,c,d
  #elif O_IONS == 5      
   #define O_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
  #endif
  
  #if Ne_IONS == 0
   #define Ne_EXPAND(a,b,c,d,e)  
  #elif Ne_IONS == 1      
   #define Ne_EXPAND(a,b,c,d,e)  ,a
  #elif Ne_IONS == 2     
   #define Ne_EXPAND(a,b,c,d,e)  ,a,b
  #elif Ne_IONS == 3       
   #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c 
  #elif Ne_IONS == 4
   #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c,d
  #elif Ne_IONS == 5      
   #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
  #endif
  
  #if S_IONS == 0
   #define S_EXPAND(a,b,c,d,e)  
  #elif S_IONS == 1      
   #define S_EXPAND(a,b,c,d,e)  ,a
  #elif S_IONS == 2     
   #define S_EXPAND(a,b,c,d,e)  ,a,b
  #elif S_IONS == 3       
   #define S_EXPAND(a,b,c,d,e)  ,a,b,c 
  #elif S_IONS == 4
   #define S_EXPAND(a,b,c,d,e)  ,a,b,c,d
  #elif S_IONS == 5      
   #define S_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
  #endif
  
  #if Fe_IONS == 0
   #define Fe_EXPAND(a,b,c)  
  #elif Fe_IONS == 1      
   #define Fe_EXPAND(a,b,c)  ,a
  #elif Fe_IONS == 2     
   #define Fe_EXPAND(a,b,c)  ,a,b
  #elif Fe_IONS == 3       
   #define Fe_EXPAND(a,b,c)  ,a,b,c 
  #endif
  
  ////////////////////////////////////////////////////////////////////////
/* Ions are labeled progressively, depending on how many ionization 
       stages are effectively included in the network through the previous 
       X_EXPAND macros. 
       Elements are ordered as {H, He, C, N, O, Ne, S, Fe} and must be 
       carefully respected everywhere in the code. 
       Hydrogen and Helium are always included.                       */
////////////////////////////////////////////////////////////////////////
  
  enum {
    X_HI = IONS_BEG, X_HeI, X_HeII
    C_EXPAND(X_CI, X_CII, X_CIII, X_CIV, X_CV)
    N_EXPAND(X_NI, X_NII, X_NIII, X_NIV, X_NV)
    O_EXPAND(X_OI, X_OII, X_OIII, X_OIV, X_OV)
    Ne_EXPAND(X_NeI, X_NeII, X_NeIII, X_NeIV, X_NeV)
    S_EXPAND(X_SI, X_SII, X_SIII, X_SIV, X_SV)
    Fe_EXPAND(X_FeI, X_FeII, X_FeIII)
  };
  
  #define NIONS  (3+C_IONS+N_IONS+O_IONS+Ne_IONS+S_IONS+Fe_IONS)
  
  double H_MassFrac (void);
  double find_N_rho ();
#endif
#endif // end if 0 not implemented yet


typedef struct COOLING_COEFFICIENTS
{
  double Lrate[NIONS];  //  Abundance increase rates. This is to be multiplied with N(x-1) - it's the ionization from x-1 to x 
  double Crate[NIONS];  //  Abundance decrease rates. This is to be multiplied with N(x) - it's the ionization/recombination from x  
  double Rrate[NIONS];  //  Abundance increase rates. This is to be multiplied with N(x+1) - it's the recombination from x+1 to x  
  double N;                    //  Total number density                 
  double Ne;                   //  Electron number density ( cm^{-3} )  
  double T;                    //  Temperature (K)                      
  double La[NIONS];     //  The coefficient of N_el in Lrate     
  double Lb[NIONS];     //  The coefficient of X_H in Lrate      
  double Lc[NIONS];     //  The coefficient of X_He in Lrate     
  double Ca[NIONS];     //  The coefficient of N_el in Crate     
  double Cb[NIONS];     //  The coefficient of X_H in Crate      
  double Cc[NIONS];     //  The coefficient of X_He in Crate     
  double Ra[NIONS];     //  The coefficient of N_el in Rrate     
  double Rb[NIONS];     //  The coefficient of X_H in Rrate      
  double Rc[NIONS];     //  The coefficient of X_He in Rrate     
  double muN, muD;        //  The numerator and denominator of the mean molecular weight - used for d\mu computations  
  double de[NIONS];     //  The radiative losses read from cooling tables, interpolated for T and N_el               
  double de_dne[NIONS]; //  The slope coefficient in the radiative losses interpolation function of N_el             

// -------- NEW ENTRIES -------- 
  
  double dnel_dX[NIONS];
  double fCH, fRH;
  double dmuN_dX[NIONS];
  double dmuD_dX[NIONS];
  double dLIR_dX[NIONS];

} COOL_COEFF;

////////////////////////////////////////////////////////////////////////
/* Function prototypes                                                */
////////////////////////////////////////////////////////////////////////

void ReadTabFile(RunConfig &);

#pragma omp declare target
#pragma acc routine seq
double CompEquil  (double, double, double *); 
#pragma omp end declare target

void CoolingSource (Data *, double, timeStep *, Grid *);

#if COOLING == POWER_LAW
void PowerLawCooling (Data &, double, timeStep &, Grid &);
#elif COOLING == SYNCH
void  SynchrotronCooling(Data &, double, timeStep &, Grid &);
#endif

// double GetMaxRate (double *, double *, double);
// void   Jacobian (double *, double *, double **);
// void   NormalizeIons (double *);
#pragma acc routine seq
void   Numerical_Jacobian (double*, double J[][(NIONS+1)], RunConfig &);

//void   Numerical_Jacobian (double*, double [NIONS+1][NIONS+1], RunConfig);

#pragma omp declare target
  #pragma acc routine seq
  void  Radiat (double *, double *, RunConfig &);
  #pragma acc routine seq
  uint32_t SolveODE_RKF12 (double *, double *, double *, double, double, RunConfig &);
  #pragma acc routine seq
  uint32_t SolveODE_RKF23 (double *, double *, double *, double, double, RunConfig &);
  #pragma acc routine seq
  uint32_t SolveODE_Ros34 (double *, double *, double *, double, double, RunConfig &);
  #pragma acc routine seq
  int StiffnessDetection(double *, double *, double, RunConfig &, int);
#pragma omp end declare target
