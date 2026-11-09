////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Relativistic magnetized Kelvin-Helmholtz problem.

  Implements the initial condition for a doubly periodic
  KH-unstable shear layer in 2D/3D Cartesian coordinates.
  The setup is the same as in [Mizuno_2013] [Mig_etal2019]

  \b References
     - [Mizuno2013]   Mizuno, ApJS (2013) 205, 7
     - [Mig_etal2012] Mignone, Mattia, Bodo, Del Zanna, MNRAS (2019) 486,4252

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Jan 31, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*!
  Assign Iinitial conditions.
*/
////////////////////////////////////////////////////////////////////////
{
  double arg;
  double *inputParam = run_config.inputParam;
  double A0    = 0.1;
  double a     = 0.01;
  double alpha = 0.1;
  double vsh   = 0.5;
  
  run_config.gamma = (4.0/3.0);
  run_config.small_prs = (1.e-4);

  arg    = (fabs(y) - 0.5)/a;
  v[VX1] = vsh*tanh(arg);
  v[RHO] = 0.5*(tanh(arg) + 1.0) + 0.5*(1.0 - tanh(arg))*0.01;

  arg    = (fabs(y) - 0.5)/alpha;

  v[VX2] = DSIGN(y)*A0*vsh*sin(2.0*M_PI*x)*exp(-arg*arg);
  #if DIMENSIONS == 2
  v[VX3] = 0.0;
  #elif DIMENSIONS == 3
  double rnd;

// IMPORTANT: 
// We avoid using the fixed-seed technique to spawn always the same sequence 
// of pseudo-random numbers.
// This may depend on how many times Startup() calls the Init() function
// Safer to use a unique generator using triplets of (x,y,z) and a 
// hashing function.

  double dx = (run_config.domEnd[IDIR] - run_config.domBeg[IDIR])/run_config.Nx1;
  double dy = (run_config.domEnd[JDIR] - run_config.domBeg[JDIR])/run_config.Nx2;
  double dz = (run_config.domEnd[KDIR] - run_config.domBeg[KDIR])/run_config.Nx3;

// Map floating-point coordinates to integers
  int64_t ix = floor((x - run_config.domBeg[IDIR])/dx);
  int64_t iy = floor((y - run_config.domBeg[JDIR])/dy);
  int64_t iz = floor((z - run_config.domBeg[KDIR])/dz);

  rnd = GridRandomNumber(ix, iy, iz);  // Generate prn in [0,1]
  rnd = (2.0*rnd-1); // Map into [-1,1]

  v[VX3] = A0*vsh*rnd*exp(-arg*arg);
  #endif

  v[PRS] = 1.0;
  v[TRC] = (1.0 - DSIGN(arg))/2.0;

  v[BX1] = sqrt(2.0*v[PRS]*inputParam[SIGMA_POL]);
  v[BX2] = 0.0;
  v[BX3] = sqrt(2.0*v[PRS]*inputParam[SIGMA_TOR]);

  #if PHYSICS == ResRMHD
  v[EX1] = -(v[VX2]*v[BX3] - v[VX3]*v[BX2]);
  v[EX2] = -(v[VX3]*v[BX1] - v[VX1]*v[BX3]);
  v[EX3] = -(v[VX1]*v[BX2] - v[VX2]*v[BX1]);
  #endif

  #ifdef GLM_MHD
  v[PSI_GLM] = 0.0;
  #ifdef PHI_GLM
  v[PHI_GLM] = 0.0;
  #endif
  #endif

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*!
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
  
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  int i,j,k, nv;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  double v[NVAR], vy_min, vy_max, Bp_tot, Bp2, dvy, dV;
  DIM_EXPAND(double Lx  = grid->xend_glob[IDIR] - grid->xbeg_glob[IDIR];,
             double Ly  = grid->xend_glob[JDIR] - grid->xbeg_glob[JDIR];,
             double Lz  = grid->xend_glob[KDIR] - grid->xbeg_glob[KDIR]; )
  double Bm2, Em2, Bm2_tot, Em2_tot, Ek, lor, Ek_tot;
  DIM_EXPAND(double *dx = grid->dx[IDIR];,
             double *dy = grid->dx[JDIR];,
             double *dz = grid->dx[KDIR]; )
  double time = data->Dts->time;
  FILE *fp;
  
  #ifdef HIGH_ORDER
  RBox box;
  BoundaryConservative (data, grid);
  RBoxDefine (grid->lbeg, grid->lend, &box);
  PointValue(data, grid, &box);
  #endif

  vy_max = vy_min = Bp_tot = Bm2_tot = Em2_tot = Ek_tot = 0.0;
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    #ifdef HIGH_ORDER
    NVAR_LOOP(nv) v[nv] = data->Vpc(i,j,k,nv);
    #else
    NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k,nv);
    #endif

    dV      = DIM_EXPAND(dx[i], *dy[j], *dz[k]);
    Bp2     = v[BX1]*v[BX1] + v[BX2]*v[BX2];
    Bm2     = v[BX1]*v[BX1] + v[BX2]*v[BX2]  + v[BX3]*v[BX3];
    Em2     = v[EX1]*v[EX1] + v[EX2]*v[EX2]  + v[EX3]*v[EX3];
    lor     = v[VX1]*v[VX1] + v[VX2]*v[VX2]  + v[VX3]*v[VX3];
    lor     = 1.0/sqrt(1.0 - lor);  // Lorentz factor
    Ek      = v[RHO]*lor*(lor - 1.0);

    Bp_tot  += sqrt(Bp2)*dV;
    Bm2_tot += 0.5*Bm2*dV;
    Em2_tot += 0.5*Em2*dV;
    Ek_tot  += Ek*dV;
    
    vy_max  = MAX(v[VX2], vy_max);
    vy_min  = MIN(v[VX2], vy_min);
  }}}

  #ifdef PARALLEL
  double gscrh;
  
  MPI_Allreduce (&vy_max, &gscrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  vy_max = gscrh;

  MPI_Allreduce (&vy_min, &gscrh, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  vy_min = gscrh;

  MPI_Allreduce (&Bp_tot, &gscrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Bp_tot = gscrh;

  MPI_Allreduce (&Bm2_tot, &gscrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Bm2_tot = gscrh;

  MPI_Allreduce (&Em2_tot, &gscrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Em2_tot = gscrh;

  MPI_Allreduce (&Ek_tot, &gscrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Ek_tot = gscrh;
  #endif

  dvy      = 0.5*(vy_max - vy_min);
  Bp_tot  /= DIM_EXPAND(Lx, *Ly, *Lz); 
  Bm2_tot /= DIM_EXPAND(Lx, *Ly, *Lz); 
  Em2_tot /= DIM_EXPAND(Lx, *Ly, *Lz); 
  Ek_tot  /= DIM_EXPAND(Lx, *Ly, *Lz); 
  
  if (g_rank == 0){
    if (g_stepNumber == 0) {
      fp = fopen ("kh_growth.dat","w");
      fprintf (fp,"# ----------------------------------------------------");
      fprintf (fp,"------------------------------------------------------\n");
      fprintf (fp,"# Kelvin Helmholtz instability\n");
      fprintf (fp,"# Resolution = %d X %d X %d\n", RuntimeGet()->npoint[IDIR],
                                                   RuntimeGet()->npoint[JDIR],
                                                   RuntimeGet()->npoint[KDIR]);
      fprintf (fp,"# Solver     = %s\n", RuntimeGet()->solv_type);
      #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
      fprintf (fp,"# divB       = CT\n");
      #elif DIVB_CONTROL == DIV_CLEANING
      fprintf (fp,"# divB       = GLM\n");
      #endif
      fprintf (fp, "#\n");
      fprintf (fp,"# %4s  %12s  %12s  %12s  %12s  %12s\n", "time","  dvy  ",
                  " <Bp> "," <Bm2> "," <Em2> ", " <Ek> ");
      fprintf (fp,"# ----------------------------------------------------");
      fprintf (fp,"------------------------------------------------------\n");
    } else {
      fp = fopen ("kh_growth.dat","a");
    }  
    
    fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", 
                  time, dvy, Bp_tot, Bm2_tot, Em2_tot, Ek_tot);
    fclose(fp);
  }  
}
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
}

