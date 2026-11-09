// /////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Header file for the radiation module.

  Set labels, indexes and basic prototyping for the radiation module.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    Nov 23, 2023
*/
// /////////////////////////////////////////////////////////////////////

// *************************************************
//     Define convenient and user-friendly
//     pointer labels for geometry setting
//   *************************************************

#if GEOMETRY == CYLINDRICAL 

  #define iFRR   FR1
  #define iFRPHI FR2
  #define iFRZ   FR3
  
#elif GEOMETRY == SPHERICAL 

  #define iFRR   FR1
  #define iFRTH  FR2
  #define iFRPHI FR3

#endif

// *************************************************
//     Radiation implicit method labels
//   *************************************************

 #define RADIATION_FIXEDPOINT_RAD    1
 #define RADIATION_NEWTON_GAS        2
 #define RADIATION_NEWTON_RAD        3
 #define RADIATION_FIXEDPOINT_GAS    4
 #define RADIATION_NEWTON_NR_GAS     5  // Newton on internal energy (NR only)
 #define RADIATION_NEWTON_NR_RAD     6  // Newton on radiation energy (NR only)


 #ifndef RADIATION_IMPL // Implicit method
  #if RADIATION_NR
   #define RADIATION_IMPL RADIATION_NEWTON_NR_GAS
  #else
   #define RADIATION_IMPL RADIATION_FIXEDPOINT_RAD
  #endif
 #endif

 // -- Convenience flag: true when using either NR Newton method --
 #if (RADIATION_IMPL == RADIATION_NEWTON_NR_GAS) \
  || (RADIATION_IMPL == RADIATION_NEWTON_NR_RAD)
  #define RADIATION_IMPLICIT_NR YES
 #else
  #define RADIATION_IMPLICIT_NR NO
 #endif

 #ifndef RADIATION_NEQS // Size of the system solved in the implicit step
  #define RADIATION_NEQS  4 
 #endif


// *************************************************
//     Additional options
//   *************************************************

 #ifndef RADIATION_IMEX_SSP2 
  #define RADIATION_IMEX_SSP2 NO
 #endif  // Apply IMEX-SSP2(2,2,2) scheme [Pareschi, L., & Russo, G. 2005,
         //            JSCom, 25, 129]. By default apply IMEX1 in [Melon Fuksman, J. D.,
         //            & Mignone, A. 2019, ApJS, 242, 20]

 #ifndef RADIATION_DIFF_LIMITING 
  #define RADIATION_DIFF_LIMITING YES
 #endif  // Limit the radiation speeds in the diffusion limit following
         //            Sadowski et al. MNRAS 429, 3533–3550 (2013)

 #ifndef RADIATION_VAR_OPACITIES
  #define RADIATION_VAR_OPACITIES NO
 #endif  // Define opacity coefficients as functions of the primitive fields.
         //            A function
         //            UserDefOpacities (double *v, double *abs_op, double *scat_op, RunConfig *run_config)
         //            must be defined in init.c  

 #ifndef RADIATION_FULL_CONVERGENCE
  #define RADIATION_FULL_CONVERGENCE NO
 #endif  // Require convergence of both radiation and RMHD

 #if IRRADIATION
  #ifndef IRRADIATION_UPDATE
   #define IRRADIATION_UPDATE NO
  #endif
 #endif  // Store space-dependent irradiation flux divergence.
         //            Compute at every step if IRRADIATION_UPDATE == YES and only at
         //            t = 0 otherwise. Must be defined in irradiation.cpp.

// *********************************************************************
/*! The Rad_data structure contains some information about the
    radiation fields, as well as some auxiliary variables used
    at the implicit radiation step.
   ********************************************************************* */
#pragma omp declare target
typedef struct RAD_DATA{
 double dt;           /**< Time step of the implicit step. */
 uint32_t * flag;     /**< Array of flags used to tag zones where convertion
                      failed or the HLLC solver needs to be replaced by HLL. */

 double * pv;         /**< Array of primitive fields. */
 double * cv;         /**< Array of conserved fields. */
 struct RunConfig_ *run_config; /**< Runtime configuration used during the implicit step. */

 double Ttot[4];       /**< Total (gas+radiation) (0,\mu) components of the
                      stress-energy tensor. */

 double Rini[4];       /**< (0,\mu) components of the radiation part of the
                      stress-energy tensor, before the first iteration of the
                      implicit step. */

 double Rprev[4];      /**< Auxiliary vector used to compute relative errors.
                      Stores the (0,\mu) components of the radiation part of the
                      stress-energy tensor, computed at the previous iteration. */

 double u[4];         /**< Spatial components of the 4-velocity. */
 double u2;           /**< Inner product u_iu^i. */
 double gamma;        /**< Lorentz factor. */
 double gamma2;       /**< Squared Lorentz factor. */

 double exv;          /**< Auxiliary variable used to compute relative errors.
                      Stores the gas pressure of the current iteration if
                      radiation fields are iterated, and the radiation energy
                      density if matter fields are iterated instead. */

 double exv_prev;     /**< Auxiliary variable used to compute relative errors.
                      Stores the gas pressure of the previous iteration if
                      radiation fields are iterated, and the radiation energy
                      density if matter fields are iterated instead. */

 char fill[20];       /**< Fill structure to power of 2 */
} Rad_data;
#pragma omp end declare target

#if RADIATION_VAR_OPACITIES == YES
#pragma omp declare target
#pragma acc routine seq
void UserDefOpacities (double *, double *, double *, RunConfig &);
#pragma omp end declare target
#endif

#include "Radiation/rad_tools.hpp"
#include "Radiation/flux_lim.hpp"
#include "Radiation/rad_errors.hpp"
#include "Radiation/rad_flux.hpp"
#include "Radiation/rad_sourcef.hpp"
#include "Radiation/rad_speed.hpp"
#include "Radiation/newtonsolve.hpp"

#pragma omp declare target
#pragma acc routine seq
uint32_t RadStep (double *, double *, double *, uint32_t *, double,
                  double, double, double, RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
uint32_t RadImplicitNR (double *, double *, double *, uint32_t *, double,
                        double, double, RunConfig &);
#pragma omp end declare target

void RadStep3D(Data *, Ary4D &, Grid *, double);
void RadStep3D(Data *, Grid *, double);

#if RADIATION_IMEX_SSP2
void AddRadSource1 (Data *, Grid *, Ary4D &, RBox *, double);
void AddRadSource2 (Data *, Grid *, Ary4D &, Ary4D &, RBox *, double, double);
#endif

#if RADIATION_NR
int AdvanceRadStep (Data *, double, timeStep *, Grid *);
void RadRightHandSide (Data *, timeStep *, const Grid *, double, RBox *, int);
void RadSubstepping   (Data *, double, timeStep *, Grid *);
void UpdateRadStage (Data *, Ary4D &, timeStep *, Grid *, RBox *, double);
#endif

#if IRRADIATION
void RadIrradiationFlux (Data *, Grid *);
#endif
void Rad_LF_Solver_T  (Data *, timeStep *, Grid *, RBox *, int);
void Rad_HLL_Solver_T  (Data *, timeStep *, Grid *, RBox *, int);
void Rad_HLLC_Solver_T (Data *, timeStep *, Grid *, RBox *, int);
