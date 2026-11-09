////////////////////////////////////////////////////////////////////////
/* gPLUTO main header file

   Contains prototypes, macros and global variables declaration for gPLUTO

   This code is distributed under the BSD 3-Clause License (see below).
  \date   May 28, 2026

  ----------------------------------------------------------------------
   BSD 3-Clause License
   Copyright (c) [2020], [Andrea Mignone]
   Copyright (c) [2020], NVIDIA CORPORATION
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  ----------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////
#define gPLUTO_VERSION  "0.91  (June 2026)"

#include <stdio.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
// #include <cstdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include <filesystem>

#include <iostream>
#include <chrono>
#include <sstream>
using namespace std;

#ifdef _OPENACC
#include <openacc.h>
#include <cuda.h>

#ifdef PARALLEL
#include <nccl.h>
#endif
#endif

#ifdef _OPENMP
  #ifdef PARALLEL
    #if (defined _AMDGPU)
      #pragma omp declare target
      #include <rccl.h>
      #pragma omp end declare target 
    #else
      #include <cuda.h>
      #include <nccl.h>
    #endif
  #endif
#endif

#define YES        1
#define NO         0
#define DEFAULT   -1
#define TRUE       YES
#define FALSE      NO

/*
#ifdef PARALLEL
  typedef MPI_Datatype  FileDataType;
  typedef MPI_Offset    FileOffset;
#else
  typedef int           FileDataType;
  typedef long int      FileOffset;
#endif
*/

// Geometry Labels ( > 0)

#define CARTESIAN    1
#define CYLINDRICAL  2
#define SPHERICAL    3

// Equation of state (EoS) labels 

#define IDEAL         1
#define PVTE_LAW      2
#define TAUB          3
#define BAROTROPIC    4
#define ISOTHERMAL    5

// Entropy switch options 

#define  SELECTIVE      1
#define  ALWAYS         2

// Time stepping labels 

#define EULER                      1
#define RK2                        2
#define RK3                        3
#define RK4                        4
#define EXP_MIDPOINT               9  // -- Used for dust time stepping -- 
#define SEMI_IMPLICIT             10  // -- Used for dust time stepping -- 
#define ARK4                      11
#define IERK45                    12

#define EXPLICIT             1 // -- just a number different from 0 !!!  -- 
#define SUPER_TIME_STEPPING  2 // -- just a number different from EXPLICIT -- 
#define RK_CHEBYSHEV         4
#define RK_LEGENDRE          8

#define IMEX        2   // Any number different from EXPLICIT (=1) 


// Output labels 

#define DBL_OUTPUT      1
#define FLT_OUTPUT      2
#define VTK_OUTPUT      3
#define VTK_ERR_OUTPUT  4
#define DBL_H5_OUTPUT   5
#define FLT_H5_OUTPUT   6
#define TAB_OUTPUT      7
#define PPM_OUTPUT      8
#define PNG_OUTPUT      9

#define PARTICLES_DBL_OUTPUT   10
#define PARTICLES_FLT_OUTPUT   11
#define PARTICLES_VTK_OUTPUT   12
#define PARTICLES_TAB_OUTPUT   13
#define PARTICLES_HDF5_OUTPUT  14

#define VTK_VECTOR  5  // -- any number but NOT 1  -- 

// Cooling labels

#define POWER_LAW    3
#define MINEq        4
#define SNEq         5
#define TABULATED    6
#define H2_COOL      7
#define KROME        8
#define SYNCH	      9

// Particle Labels

#define PARTICLES_LP    1
#define PARTICLES_CR    2
#define PARTICLES_DUST  3

// Physics modules labels

#define HD            2
#define RHD           3
#define MHD           4
#define RMHD          5
#define ResRMHD       6
#define GRMHD         7


#define MAX_STRING_SIZE    128   // Max number of characters in a string
#define MAX_OUTPUT_TYPES    16   // Max number of allowed data formats
#define MAX_OUTPUT_VARS     64   // Maximum number of variables dumped to disk
#define MAX_RESTART_PARAMS  12   // Maximum number of restart parameters

// Labels for div(B) Control: 
// If you move them to the MHD header,  definitions.hpp 
// (which is included before) cannot correctly use them 

#define EIGHT_WAVES            1
#define DIV_CLEANING           2
#define CONSTRAINED_TRANSPORT  3

// Body Force Labels: 
// Please do not change them since they are
// used in bitwise operations 

#define VECTOR     4   // corresponds to  100 in binary  
#define POTENTIAL  8   // corresponds to 1000 in binary  

// Boundary conditions labels

#define OUTFLOW          1  // any number except 0 !! 
#define REFLECTIVE       2
#define AXISYMMETRIC     3
#define EQTSYMMETRIC     4
#define PERIODIC         5
#define SHEARING         6
#define ROTATED          7
#define USERDEF          8
#define SINGULAR_AXIS    9
#define CONSTANT         10

/*! \name Labels identifying different boundary and domain regions.
   These are useful in Boundary() and when setting RBox structures
   in different points of the code.
*/
/**@{ */
#define X1_BEG        101   /**< Boundary region at X1 beg */
#define X1_END        102   /**< Boundary region at X1 end */
#define X2_BEG        103   /**< Boundary region at X2 beg */
#define X2_END        104   /**< Boundary region at X2 end */
#define X3_BEG        105   /**< Boundary region at X3 beg */
#define X3_END        106   /**< Boundary region at X3 end */
/**@} */

/*! \name Errors and flags bit labels.
    The following macros define the bits that can be turned on or off
    in an uint32_t (32 bit) variable.
    They are typically assign to the \c ***err or \c ***flag array in order to 
    bookeep, in each zone (i,j,k), which bits are actually switched on or off.

    In the case of \c ***err, the bits are reserved for error diagnostic while 
    for \c ***flag the bits are used to enable/disable certain actions in 
    specific computational cells.
    A simple bitwise operation is used to enable a flag, e.g.,
    <tt> flag(i,j,k) |= FLAG_XXX </tt>.
    For instance, by turning the ::FLAG_HLL bit on, we have
    <tt> flag = 00000100 </tt>, while by also enabling the ::FLAG_SPLIT_CELL
    one has <tt> flag = 00010100 </tt> and so on.
    Individual bits can also be turned off with the complement operator, 
    e.g., <tt> flag &= ~FLAG_XXX; </tt>
*/
/**@{ */
#define ERR_NEGATIVE_RHO                 (1<<0)  // rho < 0  (mappers) 
#define ERR_NEGATIVE_ENG                 (1<<1)  // eng < 0  (mappers) 
#define ERR_NEGATIVE_PRS                 (1<<2)  // prs < 0  (mappers) 
#define ERR_NEGATIVE_ENTR                (1<<3)  // sigma (entropy) < 0  (mappers) 
#define ERR_CONS2PRIM_FAIL               (1<<4)  // mappers 
#define ERR_PRESSURE_FIX                 (1<<5)  // Pressure fix fail (mappers) 
#define ERR_ENG_SOLVE_INACCURATE         (1<<6)  // XXX_EnergySolve()  
#define ERR_ENG_SOLVE_NO_SOLUTION        (1<<7)  // XXX_EnergySolve() 
#define ERR_SUPERLUMINAL_SPEED           (1<<8)  // v > 1 (mappers) 
#define ERR_TOO_MANY_ITERATIONS          (1<<9)  // Energy solve, IMEX 
#define ERR_RIEMANN_SOLVE_FAIL           (1<<10) // Riemann solver has failed  
#define ERR_NAN                          (1<<11) // NaN found  
#define ERR_POINT_VALUE_FAIL             (1<<12) // Point value cannot be obtained  
#define ERR_IMEX_GUESS_FAIL              (1<<13) // Guess in Imex source was difficult 
#define ERR_DIVB_FAIL                    (1<<14) // Divergence of B too large (CT) 
#define ERR_QUIT                         (1<<15)  
#define ERR_SUPERLUMINAL_RADIATION_FLUX  (1<<16)
#define ERR_GAUSSIAN_SOLVE               (1<<17) /* GaussianSolve() failed */
#define ERR_NEGATIVE_ENR                 (1<<18) /* enr < 0  (radiation) */

#define FLAG_FLAT                       (1<<0) 
#define FLAG_MINMOD                     (1<<1)
#define FLAG_HLL                        (1<<2)
#define FLAG_ENTROPY                    (1<<3)
#define FLAG_ENERGY                     (1<<4)  // Used to later unflag FLAG_ENTROPY
#define FLAG_INTERNAL_BOUNDARY          (1<<5)
#define FLAG_HO_LAP_LIMITER             (1<<6) 
#define FLAG_FALLBACK                   (1<<7) 
#define FLAG_RAD_MAX_ITER               (1<<8)  // Reached max number of iterations in rad. solver
#define FLAG_GAUSSIAN_SOLVE             (1<<9)  // GaussianSolve() returned error

#define FLAG_SHOCK                      (1<<11)

#define FLAG_COOLING_ROS34               (1<<12)  

#define ERR_NEGATIVE_SOUND_SPEED                 2  // MODULE MHD ROE 
#define ERR_MATRIX_NO_COND                       4  // MODULE MHD ROE 

/**@} */      

// Fallback options

#define FALLBACK_BLENDING         (1<<0)
#define FALLBACK_RECOMPUTE        (1<<1)

// Operations when exchanging flag() across processors 

#define FLAG_COPY   0   
#define FLAG_OR     1
#define FLAG_AND    2

#define IDIR     0     //   This sequence (0,1,2) should 
#define JDIR     1     //   never be changed             
#define KDIR     2     //                                
#define ALL_DIR -1
#define NO_DIR  -2

// Location of a variable inside the cell

#define CELL_CENTER 0  // Means (i, j, k)         
#define X1_FACE     1  // Means (i-1/2, j, k)    
#define X2_FACE     2  // Means (i, j-1/2, k)    
#define X3_FACE     3  // Means (i, j, k-1/2)    
#define X1_EDGE     4  // Means (i, j-1/2, k-1/2)
#define X2_EDGE     5  // Means (i-1/2, j, k-1/2)
#define X3_EDGE     6  // Means (i-1/2, j-1/2, k)

// Reconstruction labels 

#define FLAT             10
#define LINEAR           20
#define PARABOLIC        30
#define LimO3            40
#define WENO3            50
#define MP5              60
#define WENOZ            70

#define MULTID           3

// Limiter labels.
// These labels must be in the range [21,29] in order to allow 
// different limiters to be used during the transverse reconstruction.

#define FLAT_LIM          21
#define MINMOD_LIM        22
#define VANALBADA_LIM     23
#define OSPRE_LIM         24
#define UMIST_LIM         25
#define VANLEER_LIM       26
#define MC_LIM            27

/*! \name Physical constants in c.g.s units.
     The following set of macros express some useful physical constants
     in c.g.s units (erg, cm and sec). Values have been taken from
     http://physic.nist.gov/cuu/Constants/index.html
*/
/**@{ */
#define CONST_AH      1.008              /*!< Atomic weight of Hydrogen  */
#define CONST_AHe     4.004              /**< Atomic weight of Helium  */
#define CONST_AZ      30.0               /**< Mean atomic weight of heavy elements */
#define CONST_amu     1.66053886e-24     /**< Atomic mass unit.          */
#define CONST_au      1.49597892e13      /**< Astronomical unit.         */
#define CONST_c       2.99792458e10      /**< Speed of Light.            */
#define CONST_e       4.80320425e-10     /**< Elementary (proton) charge */
#define CONST_eV      1.602176463158e-12 /**< Electron Volt in erg.      */
#define CONST_G       6.6726e-8          /**< Gravitational Constant.    */
#define CONST_h       6.62606876e-27     /**< Planck Constant.           */
#define CONST_kB      1.3806505e-16      /**< Boltzmann constant.        */
#define CONST_ly      0.9461e18          /**< Light year.                */
#define CONST_mp      1.67262171e-24     /**< Proton mass.               */
#define CONST_mn      1.67492728e-24     /**< Neutron mass.              */
#define CONST_me      9.1093826e-28      /**< Electron mass.             */
#define CONST_mH      1.6733e-24         /**< Hydrogen atom mass.        */
#define CONST_Msun    2.e33              /**< Solar Mass.                */
#define CONST_Mearth  5.9736e27          /**< Earth Mass.                */
#define CONST_NA      6.0221367e23       /**< Avogadro Contant.          */
#define CONST_pc      3.0856775807e18    /**< Parsec.                    */
#define CONST_PI      3.14159265358979   /**< \f$ \pi \f$.               */
#define CONST_Rearth  6.378136e8         /**< Earth Radius.              */
#define CONST_Rgas    8.3144598e7        /**< Perfect gas constant       */
#define CONST_Rsun    6.96e10            /**< Solar Radius.              */
#define CONST_sigma   5.67051e-5         /**< Stephan Boltmann constant. */
#define CONST_sigmaT  6.6524e-25         /**< Thomson Cross section.    */
/**@} */

////////////////////////////////////////////////////////////////////////
/* Including header files here                                        */
////////////////////////////////////////////////////////////////////////

#include "definitions.hpp"  // Problem-dependent header file  

////////////////////////////////////////////////////////////////////////
/* Set default values of fine-tuning macro-define
    constants.
    This section of the code is for general-purpose macros
    although other may exists elsewhere.                            */
////////////////////////////////////////////////////////////////////////

#ifndef ACC_INIT
  #define ACC_INIT   YES
#endif

#if PHYSICS == HD || PHYSICS == RHD
  #undef ASSIGN_VECTOR_POTENTIAL
#endif
    
#ifndef ASSIGN_VECTOR_POTENTIAL
  #define ASSIGN_VECTOR_POTENTIAL   NO
#endif

#ifndef BACKGROUND_FIELD
 #define BACKGROUND_FIELD NO
#endif

#ifndef CHAR_LIMITING
 #define CHAR_LIMITING  NO
#endif

#ifndef CHECK_DIVB_CONDITION
  #define CHECK_DIVB_CONDITION   NO
#endif

#ifndef CONS_DBL  
  #define CONS_DBL     NO   // Save conservative var in dbl files
#endif

#ifndef COOLING  
  #define COOLING       NO
#endif

#ifndef DEBUG_MIN_MAX
 #define DEBUG_MIN_MAX  YES
#endif

#ifndef DUST_FLUID
 #define DUST_FLUID   NO
 #define IF_DUST_FLUID(a)   
#endif

#ifndef EIGHT_WAVES_REDUCED
 #define EIGHT_WAVES_REDUCED   NO
#endif

#ifndef ENTROPY_SWITCH
 #define ENTROPY_SWITCH  NO
#endif

#ifndef EOS
 #define EOS  -1
#endif

#ifndef FALLBACK_MODE
 #define FALLBACK_MODE   NO
#endif

#if DIVB_CONTROL == DIV_CLEANING
 #define GLM_MHD
#endif

#ifndef HALL_MHD
 #define HALL_MHD          NO
#endif

#ifndef INTERNAL_BOUNDARY
 #define INTERNAL_BOUNDARY   NO
#endif

#ifndef INTERNAL_BOUNDARY_CFL
 #define INTERNAL_BOUNDARY_CFL       YES
#endif

#ifndef INTERNAL_BOUNDARY_REFLECT
 #define INTERNAL_BOUNDARY_REFLECT   NO
#endif

#ifndef LIMITER
 #define LIMITER         VANLEER_LIM
#endif

#ifndef MAPPED_GRID
 #define MAPPED_GRID          NO
#endif

#ifndef MAPPED_GRID_FORCE_PERIODS
 #define MAPPED_GRID_FORCE_PERIODS  NO
#endif

#ifndef MULTIPLE_LOG_FILES
 #define MULTIPLE_LOG_FILES   NO
#endif

#ifndef USE_AL
 #define USE_AL  YES
#endif

#ifndef PARALLEL_COMM_ASYNC
 #define PARALLEL_COMM_ASYNC  YES
#endif

#ifndef RECONSTRUCT_4VEL
 #define RECONSTRUCT_4VEL     NO  // When set to YES, reconstruct 4-velocity rather
                                  // than 3-velocity (only for RHD and RMHD physics
                                  // modules)  
#endif


#ifndef REDUCED_G
  #define REDUCED_G (CONST_G)*((UNIT_DENSITY)*(UNIT_LENGTH)*(UNIT_LENGTH))/((UNIT_VELOCITY)*(UNIT_VELOCITY))   /**< scaled value of gravitational constant >*/
#endif

#ifndef RESISTIVITY
 #define RESISTIVITY           NO
#endif

#ifndef RING_AVERAGE
 #define RING_AVERAGE          NO
#endif

#if RING_AVERAGE != NO
  #ifndef RING_AVERAGE_REC
    #define RING_AVERAGE_REC  MP5
  #endif

  #ifndef RING_AVERAGE_REC_PRIM
    #define RING_AVERAGE_REC_PRIM   NO
  #endif

  #ifndef RING_AVERAGE_CT_REC
    #define RING_AVERAGE_CT_REC FLAT
  #endif

  #ifndef RING_AVERAGE_FOURIER
    #define RING_AVERAGE_FOURIER NO
  #endif
#endif

#ifndef ROTATING_FRAME
 #define ROTATING_FRAME        NO
#endif

#ifndef SHOCK_FLATTENING
  #define SHOCK_FLATTENING     NO
#endif

#ifndef UNIFORM_CARTESIAN_GRID
 #if GEOMETRY == CARTESIAN
  #define UNIFORM_CARTESIAN_GRID   YES
 #else
  #define UNIFORM_CARTESIAN_GRID   NO
 #endif
#endif

#ifndef UPDATE_VECTOR_POTENTIAL
 #define UPDATE_VECTOR_POTENTIAL  NO
#endif

#ifndef VISCOSITY
 #define VISCOSITY NO
#endif

#ifndef RADIATION
 #define RADIATION              NO
 #define RADIATION_NR           NO
 #define RADIATION_RELATIVISTIC NO
 #define IRRADIATION            NO
 #define IRRADIATION_UPDATE     NO
#endif

#ifndef VLBI_PLUTO
 #define VLBI_PLUTO   NO
#endif

#ifndef WARNING_MESSAGES
 #define WARNING_MESSAGES    1
#endif

#ifdef _OPENACC
#define ACC_MEASURE_MEM
#endif

#ifdef PARALLEL
#if USE_AL == NO
// Codes for the decomposition mode 
#define AUTO_DECOMP ((int)  9)  // Internal decomposition 
#define MPI_DECOMP  ((int) 10)  // MPI decomposition 
#define USER_DECOMP ((int) 11)  // User defined decomposition 
#define MAX_DIM     ((int)  5)
#endif
#endif

////////////////////////////////////////////////////////////////////////
/* Set HAVE_ENERGY to YES if an energy equation exists                */
////////////////////////////////////////////////////////////////////////

#if (EOS == IDEAL) || (EOS == PVTE_LAW) || (EOS == TAUB)
 #define HAVE_ENERGY       YES
#else
 #define HAVE_ENERGY       NO
#endif

// #ifdef _USE_NVTX
#define BLUE   (0x000000FF)
#define GREEN  (0x0000FF00)
#define YELLOW (0x00FFFF00)
#define SILVER (0x00C0C0C0)
#define RED (0x00FF0000)
#define PURPLE (0x00FF00FF)
#define LIGHTBLUE   (0x0066CCFF)
#define ORANGE (0x00FF6600)
#define BROWN (0x00993300)
#define DEEPBLUE (0x00000099)
// #endif

#ifdef PARALLEL      // Only for parallel computations on static grid 
 #if USE_AL == YES
 #include <al.hpp>
 #else
 #include <mpi.h>
 #endif
#endif
#pragma omp declare target
extern long int g_usedMemory;
#pragma omp end declare target

extern long int g_usedMemory;

#ifdef PARALLEL
  typedef MPI_Datatype  FileDataType;
  typedef MPI_Offset    FileOffset;
  typedef MPI_Comm      ParallelComm;
  #define FILE_DBL MPI_DOUBLE
  #define FILE_FLT MPI_FLOAT
#else
  typedef int           FileDataType;
  typedef long int      FileOffset;
  typedef int           ParallelComm;
  #define FILE_DBL 0
  #define FILE_FLT 0
#endif

#include "Array.hpp"
#include "structs.hpp"      // Structure declaration header file 

////////////////////////////////////////////////////////////////////////
/* Cooling should be included as soon as possible as it
    may change the number of variables.
    The header "cooling.hpp" defines the number of ions 
    NIONS and the array index corresponding to the 1st 
    ion IONS_BEG                                                      */
////////////////////////////////////////////////////////////////////////

#if COOLING != NO   
  #include "Cooling/cooling.hpp" 
#endif

#ifndef H_MASS_FRAC  // Set default values  
  #if COOLING == SNEq
    #define H_MASS_FRAC       0.737743792120194
  #else
    #define H_MASS_FRAC       0.7110
  #endif
#endif

#ifndef He_MASS_FRAC
  #if (EOS == PVTE_LAW) && (COOLING == NO)
    #define  He_MASS_FRAC  (1 - H_MASS_FRAC) // Effective Y and not 0.2741 Baraffe (2008) 
  #elif COOLING == SNEq
    #define He_MASS_FRAC     (0.082*CONST_AHe*H_MASS_FRAC/CONST_AH)
  #else
    #define He_MASS_FRAC      0.2741
  #endif
#endif

#define Z_MASS_FRAC (1.0 - H_MASS_FRAC - He_MASS_FRAC)
#define FRAC_He     (He_MASS_FRAC/CONST_AHe*CONST_AH/H_MASS_FRAC)
#define FRAC_Z      (Z_MASS_FRAC /CONST_AZ *CONST_AH/H_MASS_FRAC)

////////////////////////////////////////////////////////////////////////
/*! Define the total number of variables (NVAR) that
    are evolved in time:

    - \c [0:NFLX-1] number of equations defining the system
      of conservation laws. For example, for the HD module,
      it consists of density, momentum and energy.

    - \c [TRC:TRC+NTRACER-1]: number of user-defined tracers; defined
      in the problem directory header file definitions.hpp

    - \c [ENTR]: entropy (when defined)

    - \c [IONS_BEG:IONS_BEG+NIONS-1]: number of chemical species
         defined in the cooling modules cooling.h, if present.

    - \c [ENR:ENR+NRAD-1]: number of radiation fields.

    In total, there are <tt>NSCL = NTRACER+ENTROPY+NIONS</tt>
    passive scalar to be advected.

    When radiation is enabled, radiation variables are stored after
    the passive scalar block so that fluid Riemann solvers can use
    NFLX loops.

    \verbatim
     NFLX   NTRACER   ENTR     NIONS      NRAD
            <----------------------->
                       NSCL
            |TRC              |IONS_BEG
     <---------------------------------------------->
                         NVAR

    \endverbatim                                                      */
////////////////////////////////////////////////////////////////////////

#ifndef NIONS
  #define NIONS 0
#endif

#if RADIATION
  #define NRAD 4
#else
  #define NRAD 0
#endif

#define TRC      (NFLX)                                    // First tracer index
#define NSCL     (NTRACER + (ENTROPY_SWITCH != 0) + NIONS) // Total number of scalars

#if RADIATION
  #define ENR (TRC + NSCL)
  #define FR1 (ENR + 1)
  #define FR2 (ENR + 2)
  #define FR3 (ENR + 3)
  #define FIR (FR3 + 1)  // irradiation flux divergence (only used when IRRADIATION)

  #if (PHYSICS == HD) || (PHYSICS == MHD)
    #define RADIATION_NR YES
  #else
    #define RADIATION_NR NO
  #endif

  #if (PHYSICS == RHD) || (PHYSICS == RMHD)
    #define RADIATION_RELATIVISTIC YES
  #else
    #define RADIATION_RELATIVISTIC NO
  #endif
#endif

#ifndef IRRADIATION
  #define IRRADIATION NO
#endif

#define NVAR     (NFLX + NSCL + NRAD + IRRADIATION)        // Total number of variables

#if ENTROPY_SWITCH
  #define ENTR  (TRC + NTRACER)
#endif

////////////////////////////////////////////////////////////////////////
/*! These macro definition should take place afer NVAR  
    has been defined.
    Now RHOE is immediately after ions, so cooling loop can be done 
    continuously.                                                     */
////////////////////////////////////////////////////////////////////////

#if COOLING != NO
  #define NVAR_COOLING  (NVAR+1)
  #define RHOE           NVAR
  #define NVAR_COOLING_LOOP(nv) for (nv = IONS_BEG; nv < NVAR_COOLING; nv++)
#endif

#if HAVE_ENERGY
  #define IF_ENERGY(a)  a
#else 
  #define IF_ENERGY(a)  
#endif

typedef void RiemannSolver (Data *, timeStep *, Grid *, RBox *, int);

#include "macros.hpp"       // Function-like macro header file 
#include "Math_Tools/math_tools.hpp"  // Math tools header file 

#if (defined FARGO)
 #include "Fargo/fargo.hpp"           // FARGO header file 
 #define IF_FARGO(a)  a
#else
 #define IF_FARGO(a) 
#endif

////////////////////////////////////////////////////////////////////////
/* Declare global variables                                           */
////////////////////////////////////////////////////////////////////////

extern int SZ;
extern DomainDecomp new_SZ;  // new replacement
extern int SZ_stagx;
extern int SZ_stagy;
extern int SZ_stagz;
extern int SZ_float;
extern int SZ_char;
extern int SZ_Float_Vect;
extern int SZ_rgb;
extern int SZ_short;

extern int    g_intStage;
extern int    g_maxIMEXIter;
extern int    g_rank; 
#pragma acc declare create(g_rank)

extern int    g_nccl;
extern int    g_stepNumber;
#pragma acc declare create(g_stepNumber)

// #pragma omp declare target
// extern int    g_useGPU;
// #pragma omp end declare target

extern int    g_debug;

#if defined(_OPENACC) ||  defined(_OPENMP)
#ifdef PARALLEL
#pragma omp declare target
extern ncclComm_t nccl_comm;
#pragma omp end declare target
#endif
#endif

#include "prototypes.hpp"
#include "rotate.hpp"
#include <metric.hpp>

////////////////////////////////////////////////////////////////////////
/* Include main module header files                                   */
////////////////////////////////////////////////////////////////////////

#if PHYSICS == HD
  #include "HD/hd.hpp"
#elif PHYSICS == MHD
  #include "MHD/mhd.hpp"
#elif PHYSICS == RHD
  #include "RHD/rhd.hpp"
#elif PHYSICS == RMHD
  #include "RMHD/rmhd.hpp"
#elif PHYSICS == ResRMHD
  #include "ResRMHD/resrmhd.hpp"
#endif

#if PHYSICS == GRMHD
  #include "GRMHD/grmhd.hpp"
  #include "GRMHD/metric.hpp"

  // ---- Define whether the space-time is curved ---- 
  #ifndef GR_METRIC
    #define GR_METRIC  CARTESIAN
  #endif
  
  #define BOYER_LINDQUIST 5
  #define KERR_SCHILD     6
  #define USERDEF_METRIC  7

  #ifndef MINKOWSKI
    #if (GR_METRIC == BOYER_LINDQUIST) || (GR_METRIC == KERR_SCHILD) || (GR_METRIC == USERDEF_METRIC)
      #define MINKOWSKI    NO
    #else
      #define MINKOWSKI    YES
    #endif
  #endif
  #ifndef PRIMITIVE_4VEL
    #define PRIMITIVE_4VEL     NO  // When set to YES, use the 4-velocity rather
                                     // than the 3-velocity as primitive variable 
  #endif
  #ifndef GR_FD
    #define GR_FD        NO     // Finite differences formulation 
  #endif
  #ifndef BARYCENTER
    #define BARYCENTER   NO     // Compute initial condition and source terms in cell's barycenter 
  #endif
#endif


#include "plm_coeffs.hpp"      // PLM header file 

#if (PARTICLES != NO)                   // Particle Header File 
 #include "Particles/particles.hpp"
#endif

#ifdef HIGH_ORDER
 #include "High_Order/high_order.hpp"
#endif

// Decide when HLL switch is to be enabled
 
#define HLL_SWITCH     (SHOCK_FLATTENING == MULTID) \
                    || (FALLBACK_MODE    != NO)

#ifndef SELF_GRAV
  #define SELF_GRAV   NO
#endif

#ifndef POISSON_SOLVER
  #define POISSON_SOLVER   NO
#endif

#if POISSON_SOLVER == YES
	#ifndef POISSON_SOLVER_FREQ
		#define POISSON_SOLVER_FREQ 100
	#endif
	#include "PoissonSolver/poisson_solver.hpp" 
#endif

#if RADIATION
 #include "Radiation/radiation.hpp"
#endif
