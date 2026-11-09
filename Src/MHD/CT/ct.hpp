////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Header file for Constrained-Transport (CT) module.

  Provides macros, function prototypes and structure definitions for
  the constrained transport (CT) MHD module to control the divergence-free
  condition.

  \author A. Mignone (andrea.mignone@unito.it)\n
          S.Truzzi (stefano.truzzi@unito.it)
  \date   Nob 6, 2025
*/
////////////////////////////////////////////////////////////////////////
#define STAGGERED_MHD

// ---- set labels for CT_EMF_AVERAGE ---- 

#define ARITHMETIC        1
#define UCT0              2
#define CT_CONTACT        3
#define CT_FLUX           4
#define UCT_HLL           5
#define UCT_HLLD          6
#define UCT_GFORCE        7
#define CT_MAXWELL        8

// ---- staggered component labels ---- 

#define BX1s  0
#if INCLUDE_JDIR
  #define BX2s  1
#endif
#if INCLUDE_KDIR
  #define BX3s  2
#endif

#if (PHYSICS == ResRMHD)
  #define   EX1s  3 
  #if INCLUDE_JDIR
    #define EX2s  (EX1s + 1)
  #endif
  #if INCLUDE_KDIR
    #define EX3s  (EX1s + 2)
  #endif
#endif

// ---- backward compatibility ---- 

#define BXs BX1s
#define BYs BX2s
#define BZs BX3s

////////////////////////////////////////////////////////////////////////
/* Macros                                                              */
////////////////////////////////////////////////////////////////////////

/*! \name CT_ASSIGN macro.
    The following macros is used to avoid too many pre-processor 
    directives when the same instruction has to replicated for 
    the electric field in ResRMHD. 
    When both B and E use constrained transport, the macro evaluates
    both argument. Otherwise only the 1st one is used.
  */

/**@{ */
#if PHYSICS == ResRMHD && DIVE_CONTROL == CONSTRAINED_TRANSPORT
  #define CT_ASSIGN(a, b)    a   b
#else  
  #define CT_ASSIGN(a, b)    a
#endif

#if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  #define CT_EXPAND(k, a) DIM_EXPAND(k = BX1s; a; , k = BX2s; a; , k = BX3s; a;)  \
                          DIM_EXPAND(k = EX1s; a; , k = EX2s; a; , k = EX3s; a;)

#else
  #define CT_EXPAND(k, a)   DIM_EXPAND(k = BX1s; a; , \
                                       k = BX2s; a; , \
                                       k = BX3s; a;)  
#endif
/**@} */

////////////////////////////////////////////////////////////////////////
/* Default values                                                      */
////////////////////////////////////////////////////////////////////////

#if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
#ifndef CT_EMF_AVERAGE
  #if PHYSICS == ResRMHD
    #define  CT_EMF_AVERAGE         CT_MAXWELL
  #else
    #define  CT_EMF_AVERAGE         UCT_HLL
  #endif
#endif
#endif

#ifndef  CT_ENERGY_CORRECTION
 #define  CT_ENERGY_CORRECTION       0
#endif

#ifndef CT_EMF_BOUNDARY 
 #define CT_EMF_BOUNDARY            NO
#endif

#ifndef  CT_MAGNETIC_ENERGY
 #define  CT_MAGNETIC_ENERGY         0 // see CT_MagneticEnergy()
#endif

#define eps_CT_CONTACT   1.e-6

////////////////////////////////////////////////////////////////////////
/* \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping                                               */
////////////////////////////////////////////////////////////////////////

void CT_AllocateEMF (EMF *, Grid *);
void CT_AverageStaggeredFields (Data *, int, RBox *, Grid *);
void CT_CheckDivB (const Ary3D (&)[6], Grid *, RBox *);
void CT_ComputeEMF (const Data *, Grid *, RBox *);
void CT_ComputeCenterEMF(const Data *, const Grid *);
void CT_EMF_Boundary(const Data *, Grid *);
void CT_EMF_ArithmeticAverage (const EMF *, const double, RBox *);
void CT_EMF_IntegrateToCorner (const Data *, const EMF *, RBox *);
void CT_EMF_Riemann2D (const Data *, const EMF *, Grid *, RBox *);
void CT_EMF_Flux(const Data *, const EMF *, Grid *, RBox *);

void   CT_EnergyCorrection (Data *, Ary4D &, double, Grid *, RBox *);
double HD_EnergyFlux(double *, double *, int, RunConfig &);
double CT_MagneticEnergy(const Data *, int, int, int);

void CT_FillStaggeredFields (const Ary3D (&)[6], const Ary3D (&), int, Grid *);
void CT_GetEMF (const Data *, Grid *);
#if PHYSICS == ResRMHD || PHYSICS == GRMHD
void CT_ComputeCharge(Data *, RBox *, Grid *);
void CT_InterfaceCurrent(Data *, Grid *);
void CT_IMEXImplicitUpdate(Data *, Ary4D &, double, Grid *);
void CT_MaxwellSolver(const Data *, const EMF *, Grid *, RBox *);
void CT_ResistiveUpdate(const Data *, double, Grid *);
#endif
void CT_ResistiveEMF (const Data *, int, Grid *, RBox *);
void CT_Update(const Data *, Ary3D (&)[6], double, Grid *, RBox *);

// \endcond 
