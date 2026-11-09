////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Collects global variables definitions.

  This file contains definitions for all global variables
  (visible anywhere in the code) used by PLUTO.
  Global variables names, by convention, are prefixed with a "g_"
  unless they're used as constants throughout the code in which case
  they keep the full-capitalized notation typical of macros.

  For modules, global variables are prefixed with the initial letters
  of the module name, e.g., \c sb_vy.

  In the following "local" means "for the local processor".
  "Interior" means inside the computational domain.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Nov 06, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

int SZ;
DomainDecomp new_SZ;  // new replacement
int SZ_stagx;
int SZ_stagy;
int SZ_stagz;
int SZ_float;
int SZ_char;
int SZ_Float_Vect;
int SZ_rgb;
int SZ_short;

int    g_intStage;
#pragma acc declare create(g_intStage)

int    g_maxIMEXIter;
#pragma acc declare create(g_maxIMEXIter)

int    g_rank;
// #pragma acc declare create(g_rank)

int    g_stepNumber;
#pragma acc update device(g_stepNumber)

int    g_nccl=0;
#pragma acc declare copyin(g_nccl)

// double g_time;
// double g_dt;

long int g_usedMemory = 0;

int    g_debug  = 0;
// #pragma omp declare target  
// int    g_useGPU = 0;
// #pragma acc declare copyin(g_useGPU)
// #pragma omp end declare target

#if defined(_OPENACC) || defined(_OPENMP)
#ifdef PARALLEL
#pragma omp declare target
ncclComm_t nccl_comm;
#pragma omp end declare target
#endif
#endif
