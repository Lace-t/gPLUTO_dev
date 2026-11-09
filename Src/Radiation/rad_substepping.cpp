////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Substepping integration routines.

  Integrate the inhomogeneous radiation transport equations
  by carrying IMEX integrations with dt=dt_rad until t+dt_hydro
  is reached.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    Jan 13, 2020
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if RADIATION_NR

////////////////////////////////////////////////////////////////////////
/*!
  Integrate radiation transport equations over the hydro time step dt
  by carrying IMEX sub-steps with adaptive step size dt_rad.

  \param[in,out] data  pointer to Data structure
  \param[in]     dt    hydro time step to be covered
  \param[in,out] Dts   pointer to the hydro time step structure
  \param[in]     grid  pointer to Grid structure
*/
////////////////////////////////////////////////////////////////////////
void RadSubstepping (Data *data, double dt, timeStep *Dts, Grid *grid)
{
  static timeStep Dts_rad;
  static double dt_rad;
  double t, t_new, dt_hydro, dt_rad_next, dt_rad_save;
  const Runtime *runtime         = RuntimeGet();
  const double   cfl             = runtime->rad_cfl;
  const double   rad_first_dt    = runtime->rad_first_dt;
  const double   rad_cfl_max_var = runtime->rad_cfl_max_var;

  if (Dts_rad.cmax == NULL){
    int *ntot = grid->ntot;
    Dts_rad.cmax      = ARRAY_1D(ntot[IDIR]*ntot[JDIR]*ntot[KDIR], double);
    Dts_rad.invDt_hyp = 0.0;
    #pragma acc enter data copyin(Dts_rad)
    #pragma omp target enter data map(to: Dts_rad)
    dt_rad = rad_first_dt;
  }

  dt_hydro = dt;
  t = 0.0;

  while (t < dt_hydro){
    // Update t_new
    t_new = t + dt_rad;

    // Save current time step
    dt_rad_save = dt_rad;

    // Correct time step if needed
    if (t_new > dt_hydro){
      dt_rad = dt_hydro - t;
      t += dt_hydro;
    } else {
      t = t_new;
    }

    // IMEX radiation step
    AdvanceRadStep(data, dt_rad, &Dts_rad, grid);

    #ifdef PARALLEL
    {
      double xloc = Dts_rad.invDt_hyp;
      double xglob;
      MPI_Allreduce(&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      Dts_rad.invDt_hyp = xglob;
    }
    #endif

    // Update time step
    dt_rad_next = cfl / Dts_rad.invDt_hyp;

    // Constant time step if rad_cfl_max_var == 1.0, update otherwise
    if (rad_cfl_max_var != 1.0){
      dt_rad = MIN(dt_rad_next, rad_cfl_max_var*dt_rad_save);
    } else {
      dt_rad = dt_rad_save;
    }

    // Exit if dt_rad_next is too small
    if (dt_rad_next < rad_first_dt*1.e-9){
      print("! RadSubstepping(): dt is too small (%12.6e).\n", dt_rad_next);
      print("! Cannot continue.\n");
      QuitPLUTO("from RadSubstepping()", 1);
    }

    // Reset hyperbolic time step
    Dts_rad.invDt_hyp = 0.0;
  }
}

#endif // RADIATION_NR
