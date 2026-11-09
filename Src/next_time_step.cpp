////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Compute the next time step

  Compute and return the time step for the next time level using the 
  available information from the previous integration 
  (Dts->invDt_hyp and Dts->invDt_par, etc.).

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Nov 24, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
double NextTimeStep (timeStep *Dts, RunConfig *run_config, 
                     Runtime *runtime, Grid *grid)
{
  double dt_new;

  // ------------------------------------------------------
  // 1. Parallel reduction
  // ------------------------------------------------------

  #ifdef PARALLEL

    // 1a. Hyperbolic time step reduction

    // MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // Dts->invDt_hyp = xglob;

    MPI_Allreduce (MPI_IN_PLACE, &Dts->invDt_hyp, 1, 
                   MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // 1b. Parabolic time step reduction

    #if (RESISTIVITY)
    MPI_Allreduce (MPI_IN_PLACE, &Dts->invDt_par, 1, 
                   MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #endif

    // 1c. Cooling time step reduction

    #if COOLING != NO
    // xloc = Dts->dt_cool;
    // MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    // Dts->dt_cool = xglob;

    MPI_Allreduce (MPI_IN_PLACE, &Dts->dt_cool, 1, 
                   MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    #endif

    // 1d. Particle time step reduction

    #if PARTICLES != NO
    // xloc = Dts->invDt_particles;
    // MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // Dts->invDt_particles = xglob;

    MPI_Allreduce (MPI_IN_PLACE, &Dts->invDt_particles, 1, 
                   MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #endif
  #endif

  // ------------------------------------------------------
  // 2. Time step determination
  // ------------------------------------------------------

  #if PARTICLES_DT_CONTROL == YES
  Dts->invDt_hyp = std::max(Dts->invDt_hyp, Dts->invDt_particles);
  #endif
  
  #if (RESISTIVITY)
  dt_new  = runtime->cfl/(Dts->invDt_hyp + 2.0*Dts->invDt_par);  
  #else
  dt_new  = runtime->cfl/Dts->invDt_hyp;  
  #endif

  #if COOLING != NO
  dt_new= MIN( dt_new, Dts->dt_cool);
  #endif

  // ------------------------------------------------------
  // 3. For GLM MHD, recompute glm_ch from dt_new
  // ------------------------------------------------------

  #ifdef GLM_MHD
    int idim;
    double dxmin = 1.e18;

    DIM_LOOP(idim){ //  Min cell length   
      dxmin = MIN(dxmin, grid->dl_min[idim]);
    }

    #if PHYSICS != ResRMHD
    run_config->glm_ch = runtime->cfl*dxmin/dt_new;

    #ifdef PARALLEL
    double xglob;
    double xloc = run_config->glm_ch;
    MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    run_config->glm_ch = xglob;
    #endif

    #pragma acc update device(run_config->glm_ch)
    #pragma omp target update to(run_config->glm_ch)

    #else
    #endif  // PHYSICS != ResRMHD
  #else
  #endif  // GLM_MHD 

  // ------------------------------------------------------
  // 5. Limit time step variations
  // ------------------------------------------------------

  dt_new = MIN(dt_new, runtime->cfl_max_var*Dts->dt);    
  // #pragma acc update device(Dts)
  // #pragma omp target update to(Dts)    

  return dt_new;
}
