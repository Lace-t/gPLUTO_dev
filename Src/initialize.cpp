////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Initialize PLUTO.

  Initialize() performs a number of initialization tasks before 
  starting the main computation loop.\n
  More precisely, it completes the following sequence of steps:

  - runtime initialization (i.e. pluto.ini)
  - parallel domain decomposition
  - grid generation
  - memory allocation
  - initial conditions
  - set output attributes

  The function GetDecompMode() sets the parallel domain decomposition 
  mode which can be equal to
  
  - AL_AUTO_DECOMP   default;
  - AL_USER_DECOMP   if the -dec n1 [n2] [n3] command line 
                     argument has been given;
                     In this case only, procs[] contains the number
                     of processors in the three directions;
  - AL_MPI_DECOMP   [todo]

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 19, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

int GetDecompMode (cmdLine *cmd_line, int procs[]);
void CheckConfig(void);
#ifdef PARALLEL
#if USE_AL == NO
void DecompDomain(MPI_Comm, MPI_Datatype,
                  int , int *, int *,
                  int *, int *, int *,
                  int , DomainDecomp *);
#endif
#endif

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Initialize(Data *data, timeStep *Dts, cmdLine *cmd_line,
                Runtime *runtime, RunConfig *run_config, Grid *grid)
{
  int *ntot = grid->ntot;
  int nv;
  int  procs[3]   = {1,1,1};

  #ifdef PARALLEL
  int  nprocs = 1, decomp_mode;
  int  idim;
  int  nghost;
  int  gsize[3] = {1,1,1};
  int  lsize[3] = {1,1,1};
  int  beg[3] = {0,0,0};
  int  end[3] = {0,0,0};
  int  gbeg[3] = {0,0,0};
  int  gend[3] = {0,0,0};
  int  lbeg[3] = {0,0,0};
  int  lend[3] = {0,0,0};
  int  is_gbeg[3], is_gend[3];
  int  ghosts[3] = {0,0,0};
  int  periods[3] = {0,0,0};
  int  pardim[3]  = {0,0,0};
  #endif

  // ------------------------------------------------------
  // 1. Parallel domain decomposition
  // ------------------------------------------------------

#ifdef PARALLEL

  // 1a. Get number of processors 

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  print ("> nprocs: %d\n",nprocs),

  // 1b. Get number of ghost zones and set periodic boundaries 

  nghost = GetNghost();
  MPI_Allreduce (&nghost, &idim, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  nghost = idim;

  for (idim = 0; idim < DIMENSIONS; idim++){
    gsize[idim]   = runtime->npoint[idim];
    ghosts[idim]  = nghost;
    periods[idim] =    (runtime->left_bound[idim] == PERIODIC ? 1:0)
                    || (runtime->left_bound[idim] == SHEARING ? 1:0);
    pardim[idim]  = cmd_line->parallel_dim[idim];
  }

  #if !INCLUDE_JDIR
  ghosts[JDIR] = 0;
  #endif
  #if !INCLUDE_KDIR
  ghosts[KDIR] = 0;
  #endif

  // 1c. Find parallel decomposition mode and number of processors 

  decomp_mode = GetDecompMode(cmd_line, procs);

  // 1d. Set double distributed array descriptor

#if USE_AL == NO
// 1d. Set double distributed array descriptor
  DecompDomain(MPI_COMM_WORLD, MPI_DOUBLE, DIMENSIONS,
               gsize, ghosts, periods, pardim, procs,
               decomp_mode, &new_SZ);
  // Read results back into local variables (used by SetGrid etc.)
  for (int i = 0; i < DIMENSIONS; i++) {
    lsize[i] = new_SZ.lsize[i];
    beg[i]   = new_SZ.beg[i];
    end[i]   = new_SZ.end[i];
    lbeg[i]  = new_SZ.lbeg[i];
    lend[i]  = new_SZ.lend[i];
    gbeg[i]  = new_SZ.gbeg[i];
    gend[i]  = new_SZ.gend[i];
    is_gbeg[i] = new_SZ.is_gbeg[i];
    is_gend[i] = new_SZ.is_gend[i];
    procs[i] = new_SZ.cart_dims[i];
  }
#else
  AL_Sz_init (MPI_COMM_WORLD, &SZ);
  AL_Set_type (MPI_DOUBLE, 1, SZ);
  AL_Set_dimensions (DIMENSIONS, SZ);
  AL_Set_global_dim (gsize, SZ);
  AL_Set_ghosts (ghosts, SZ);
  AL_Set_periodic_dim (periods, SZ);
  AL_Set_parallel_dim (pardim, SZ);

  AL_Decompose (SZ, procs, decomp_mode);
  AL_Get_local_dim (SZ, lsize);
  AL_Get_bounds (SZ, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ, is_gbeg, is_gend);

  // 1e. Set float distributed array descriptor

  AL_Sz_init (MPI_COMM_WORLD, &SZ_float);
  AL_Set_type (MPI_FLOAT, 1, SZ_float);
  AL_Set_dimensions (DIMENSIONS, SZ_float);
  AL_Set_global_dim (gsize, SZ_float);
  AL_Set_ghosts (ghosts, SZ_float);
  AL_Set_periodic_dim (periods, SZ_float);
  AL_Set_parallel_dim (pardim, SZ_float);

  AL_Decompose (SZ_float, procs, decomp_mode);
  AL_Get_local_dim (SZ_float, lsize);
  AL_Get_bounds (SZ_float, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_float, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_float, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_float, is_gbeg, is_gend);
#endif

#if PHYSICS==GRMHD
  SetGrid(runtime, procs, grid, 1, *run_config);
#endif

  #ifdef STAGGERED_MHD
#if USE_AL == NO
// // 1g. x-staggered array descriptor (double)
// #if INCLUDE_IDIR
// for (idim = 0; idim < DIMENSIONS; idim++) gsize[idim] = runtime->npoint[idim];
// gsize[IDIR] += 1;

// #ifdef SHEARINGBOX
// periods[IDIR] = 0;
// #endif

// DecompDomain(MPI_COMM_WORLD, MPI_DOUBLE, DIMENSIONS,
//               gsize, ghosts, periods, pardim, procs,
//               decomp_mode, &new_SZ_stagx);
// new_SZ_stagx.isstaggered[IDIR] = 1;

// #ifdef SHEARINGBOX
// periods[IDIR] = 1;
// #endif
// #endif  // INCLUDE_IDIR

// // 1h. y-staggered array descriptor (double)
// #if INCLUDE_JDIR
// for (idim = 0; idim < DIMENSIONS; idim++) gsize[idim] = runtime->npoint[idim];
// gsize[JDIR] += 1;

// DecompDomain(MPI_COMM_WORLD, MPI_DOUBLE, DIMENSIONS,
//               gsize, ghosts, periods, pardim, procs,
//               decomp_mode, &new_SZ_stagy);
// new_SZ_stagy.isstaggered[JDIR] = 1;
// #endif  // INCLUDE_JDIR

// // 1i. z-staggered array descriptor (double)
// #if INCLUDE_KDIR
// for (idim = 0; idim < DIMENSIONS; idim++) gsize[idim] = runtime->npoint[idim];
// gsize[KDIR] += 1;

// DecompDomain(MPI_COMM_WORLD, MPI_DOUBLE, DIMENSIONS,
//               gsize, ghosts, periods, pardim, procs,
//               decomp_mode, &new_SZ_stagz);
// new_SZ_stagz.isstaggered[KDIR] = 1;
// #endif  // INCLUDE_KDIR
#else
  int  stagdim[3] = {0,0,0};
  
  // 1g. x-staggered array descriptor (double)

  #if INCLUDE_IDIR
  stagdim[IDIR] = AL_TRUE;
  stagdim[JDIR] = AL_FALSE;
  stagdim[KDIR] = AL_FALSE;

  for (idim = 0; idim < DIMENSIONS; idim++) {
    gsize[idim] = runtime->npoint[idim];
  }
  gsize[IDIR] += 1;

  #ifdef SHEARINGBOX
  periods[IDIR] = 0;
  #endif

  AL_Sz_init (MPI_COMM_WORLD, &SZ_stagx);
  AL_Set_type (MPI_DOUBLE, 1, SZ_stagx);
  AL_Set_dimensions (DIMENSIONS, SZ_stagx);
  AL_Set_global_dim (gsize, SZ_stagx);
  AL_Set_ghosts (ghosts, SZ_stagx);
  AL_Set_staggered_dim(stagdim, SZ_stagx);
  AL_Set_periodic_dim (periods, SZ_stagx);
  AL_Set_parallel_dim (pardim, SZ_stagx);

  AL_Decompose (SZ_stagx, procs, decomp_mode);
  AL_Get_local_dim (SZ_stagx, lsize);
  AL_Get_bounds (SZ_stagx, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_stagx, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_stagx, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_stagx, is_gbeg, is_gend);

  #ifdef SHEARINGBOX
  periods[IDIR] = 1;
  #endif
  #endif  // INCLUDE_IDIR

  #if INCLUDE_JDIR
  stagdim[IDIR] = AL_FALSE;
  stagdim[JDIR] = AL_TRUE;
  stagdim[KDIR] = AL_FALSE;

  for (idim = 0; idim < DIMENSIONS; idim++) {
    gsize[idim] = runtime->npoint[idim];
  }

  gsize[JDIR] += 1;

  AL_Sz_init (MPI_COMM_WORLD, &SZ_stagy);
  AL_Set_type (MPI_DOUBLE, 1, SZ_stagy);
  AL_Set_dimensions (DIMENSIONS, SZ_stagy);
  AL_Set_global_dim (gsize, SZ_stagy);
  AL_Set_ghosts (ghosts, SZ_stagy);
  AL_Set_staggered_dim(stagdim, SZ_stagy);
  AL_Set_periodic_dim (periods, SZ_stagy);
  AL_Set_parallel_dim (pardim, SZ_stagy);

  AL_Decompose (SZ_stagy, procs, decomp_mode);
  AL_Get_local_dim (SZ_stagy, lsize);
  AL_Get_bounds  (SZ_stagy, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_stagy, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_stagy, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_stagy, is_gbeg, is_gend);
  #endif // INCLUDE_JDIR

  #if INCLUDE_KDIR
  stagdim[IDIR] = AL_FALSE;
  stagdim[JDIR] = AL_FALSE;
  stagdim[KDIR] = AL_TRUE;

  for (idim = 0; idim < DIMENSIONS; idim++) {
    gsize[idim] = runtime->npoint[idim];
  }
  gsize[KDIR] += 1;

  AL_Sz_init (MPI_COMM_WORLD, &SZ_stagz);
  AL_Set_type (MPI_DOUBLE, 1, SZ_stagz);
  AL_Set_dimensions (DIMENSIONS, SZ_stagz);
  AL_Set_global_dim (gsize, SZ_stagz);
  AL_Set_ghosts (ghosts, SZ_stagz);
  AL_Set_staggered_dim(stagdim, SZ_stagz);
  AL_Set_periodic_dim (periods, SZ_stagz);
  AL_Set_parallel_dim (pardim, SZ_stagz);

  AL_Decompose (SZ_stagz, procs, decomp_mode);
  AL_Get_local_dim (SZ_stagz, lsize);
  AL_Get_bounds  (SZ_stagz, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_stagz, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_stagz, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_stagz, is_gbeg, is_gend);
  #endif // INCLUDE_KDIR 
  #endif // STAGGERED_MHD 
#endif // USE_AL
#endif

  // ------------------------------------------------------
  // 2. Initialize grid
  // ------------------------------------------------------

  int shrink = 1;
#if PHYSICS!=GRMHD
  SetGrid (runtime, procs, grid, shrink);
#else
  SetGrid(runtime, procs, grid, shrink, *run_config);
#endif
  OutputGrid(grid,runtime);
  
  // 2a. Find processors coordinates in a Cartesian topology 

  #ifdef PARALLEL
  {
    int coords[3] = {0,0,0};
    int rank;
    MPI_Comm cartcomm;
  #if USE_AL == NO
  cartcomm = new_SZ.cart_comm;
  #else
    AL_Get_cart_comm(SZ, &cartcomm);
  #endif
    MPI_Cart_get(cartcomm, DIMENSIONS, procs, periods, coords);
    MPI_Cart_rank(cartcomm, coords, &rank);
    if (rank != g_rank) {
      QuitPLUTO ("! Initialize: rank and g_rank are different",1);
      // printf ("! Initialize: rank and g_rank are different\n");
      // QUIT_PLUTO(1);
    }
    for (idim = 0; idim < DIMENSIONS; idim++) {
      grid->rank_coord[idim] = coords[idim];
    }
  }
  #endif


  // ------------------------------------------------------
  // 3. Make sure all processors have the same domain size
  // ------------------------------------------------------

  #ifdef PARALLEL
  int ntot_max[3];
  int ntot_min[3];

  MPI_Allreduce (ntot, ntot_max, 3, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (ntot, ntot_min, 3, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  DIM_LOOP(nv) {
    if (ntot_max[nv] != ntot_min[nv]){
      QuitPLUTO("! Initialize(): processors not equally decomposed",1);
    }
  }
  #endif

  // ------------------------------------------------------
  // 4. Initialize time step info
  // ------------------------------------------------------

  Dts->invDt_hyp = 0.0;
  Dts->invDt_par = 0.0;
  Dts->dt        = runtime->first_dt;  // Initial time step 
  Dts->cmax = ARRAY_1D(  ntot[IDIR]*ntot[JDIR]*ntot[KDIR], double);
  Dts->Dmax = ARRAY_1D(3*ntot[IDIR]*ntot[JDIR]*ntot[KDIR], double);

  // ------------------------------------------------------
  // 5. Compute chunk size for the RING_AVERAGE 
  //    technique.
  // ------------------------------------------------------

  #if RING_AVERAGE > 1
    #if GEOMETRY == CYLINDRICAL
    grid->ring_av_csize = ARRAY_1D(grid->ntot[IDIR] , int);
    #elif GEOMETRY == SPHERICAL
    grid->ring_av_csize = ARRAY_1D(grid->ntot[JDIR] , int);
    #endif
    RingAverageSize(grid);
  #endif

  // ------------------------------------------------------
  // 6. Memory allocation
  // ------------------------------------------------------

  print ("\n> Memory allocation\n");

  // 5a. Main arrays

  data->Uc.create(ntot[IDIR],ntot[JDIR],ntot[KDIR],NVAR);
  data->Uc0.create(ntot[IDIR],ntot[JDIR],ntot[KDIR],NVAR);
  data->Vc.create(ntot[IDIR],ntot[JDIR],ntot[KDIR],NVAR);

  data->integration_mode = DEFAULT;
  data->err.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  data->sweep = (Sweep*)malloc(sizeof(Sweep));
  AllocateSweepArrays (data->sweep, grid);

  for (nv = 0; nv < MAX_RESTART_PARAMS; nv++) data->restartParam[nv] = 0.0;

  // 5b. High-order additional arrays

  #ifdef HIGH_ORDER
  data->Upc.create(ntot[IDIR],ntot[JDIR],ntot[KDIR],NVAR);
  data->Vpc.create(ntot[IDIR],ntot[JDIR],ntot[KDIR],NVAR);
  data->ho_lim.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #ifdef STAGGERED_MHD
  DIM_EXPAND(data->Vps[BX1s].create((ntot[IDIR]+1),ntot[JDIR],ntot[KDIR]); ,
             data->Vps[BX2s].create(ntot[IDIR],(ntot[JDIR]+1),ntot[KDIR]); ,
             data->Vps[BX3s].create(ntot[IDIR],ntot[JDIR],(ntot[KDIR]+1));)
  #endif
  #endif

  // 5c. Constrained transport arrays

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  DIM_EXPAND(data->Vs[BX1s].create((ntot[IDIR]+1),ntot[JDIR],ntot[KDIR]); 
             data->Vs[BX1s].location_ = X1_FACE;                          ,

             data->Vs[BX2s].create(ntot[IDIR],(ntot[JDIR]+1),ntot[KDIR]); 
             data->Vs[BX2s].location_ = X2_FACE;                          ,

             data->Vs[BX3s].create(ntot[IDIR],ntot[JDIR],(ntot[KDIR]+1));
             data->Vs[BX3s].location_ = X3_FACE;)

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  DIM_EXPAND(data->Vs[EX1s].create((ntot[IDIR]+1),ntot[JDIR],ntot[KDIR]); ,
             data->Vs[EX2s].create(ntot[IDIR],(ntot[JDIR]+1),ntot[KDIR]); ,
             data->Vs[EX3s].create(ntot[IDIR],ntot[JDIR],(ntot[KDIR]+1));)
  #endif
  #if (PARTICLES_LP_MREC_SPECTRA == YES) || ((PHYSICS == ResRMHD) \
                      && (DIVE_CONTROL == CONSTRAINED_TRANSPORT))
  data->q.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif
  data->emf = (EMF*)malloc(sizeof(EMF));
  CT_AllocateEMF (data->emf, grid);
  #endif  // DIVB_CONTROL == CONSTRAINED_TRANSPORT

  // Electric field & Vector potential allocation 

  #if    (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) \
      || (PHYSICS == GRMHD)
  data->Ex1.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]); 
  data->Ex2.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->Ex3.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  data->A[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->A[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->A[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  data->Ex1.location_ = data->A[IDIR].location_ = X1_EDGE;
  data->Ex2.location_ = data->A[JDIR].location_ = X2_EDGE;
  data->Ex3.location_ = data->A[KDIR].location_ = X3_EDGE;
  #endif

  #if PARTICLES != NO
    // #if PARTICLES == PARTICLES_CR
    // data->Fcr = ARRAY_4D(4, Nntot[IDIR],ntot[JDIR],ntot[KDIR], double);
    // data->Jcr = ARRAY_4D(3, ntot[IDIR],ntot[JDIR],ntot[KDIR], double);
    // data->qcr = ARRAY_3D(ntot[IDIR],ntot[JDIR],ntot[KDIR], double);
    // #endif
    #if PARTICLES == PARTICLES_DUST
    data->Fdust = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  #endif

  // 5d. Other arrays 

  data->flag.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->C_dt.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  data->C_dtpar[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);      
  data->C_dtpar[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);      
  data->C_dtpar[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);  

  data->phi_m.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  // 5e. Self-gravity / Poisson solver related arrays
  
  #if (SELF_GRAV == YES)
    data->acc_grav.create(3,ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif

  #if (POISSON_SOLVER)
    data->poisson_sol.create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
    data->poisson_sol.setValue(0.);
  #endif

  // 5f. Non-ideal MHD related arrays
  
  data->J[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->J[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->J[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);

  #if RESISTIVITY 
  data->eta[IDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->eta[JDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  data->eta[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]);
  #endif
  
  #if RING_AVERAGE != NO 
  data->Phi[IDIR].create(ntot[IDIR]+INCLUDE_IDIR,ntot[JDIR],ntot[KDIR]);
  data->Phi[JDIR].create(ntot[IDIR],ntot[JDIR]+INCLUDE_JDIR,ntot[KDIR]);
  data->Phi[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]+INCLUDE_KDIR);

  data->PhiM[IDIR].create(ntot[IDIR]+INCLUDE_IDIR,ntot[JDIR],ntot[KDIR]);
  data->PhiM[JDIR].create(ntot[IDIR],ntot[JDIR]+INCLUDE_JDIR,ntot[KDIR]);
  data->PhiM[KDIR].create(ntot[IDIR],ntot[JDIR],ntot[KDIR]+INCLUDE_KDIR);

  #if RING_AVERAGE_FOURIER == YES
  #if GEOMETRY == CYLINDRICAL
  data->A0[IDIR].create (ntot[IDIR], ntot[KDIR]+INCLUDE_KDIR);
  data->A1[IDIR].create (ntot[IDIR], ntot[KDIR]+INCLUDE_KDIR);
  data->B1[IDIR].create (ntot[IDIR], ntot[KDIR]+INCLUDE_KDIR);
  #if INCLUDE_KDIR == YES
  data->A0[KDIR].create (ntot[IDIR], ntot[KDIR]+INCLUDE_KDIR);
  data->A1[KDIR].create (ntot[IDIR], ntot[KDIR]+INCLUDE_KDIR);
  data->B1[KDIR].create (ntot[IDIR], ntot[KDIR]+INCLUDE_KDIR);
  #endif
  #elif GEOMETRY == SPHERICAL
  data->A0[JDIR].create (ntot[IDIR]+INCLUDE_IDIR, ntot[JDIR]);
  data->A1[JDIR].create (ntot[IDIR]+INCLUDE_IDIR, ntot[JDIR]);
  data->B1[JDIR].create (ntot[IDIR]+INCLUDE_IDIR, ntot[JDIR]);

  data->A0[IDIR].create (ntot[IDIR]+INCLUDE_IDIR, ntot[JDIR]);
  data->A1[IDIR].create (ntot[IDIR]+INCLUDE_IDIR, ntot[JDIR]);
  data->B1[IDIR].create (ntot[IDIR]+INCLUDE_IDIR, ntot[JDIR]);
  #endif
  #endif
  #endif

  // ------------------------------------------------------
  // 7. Set Run_Config structure
  // ------------------------------------------------------

  #if EOS == IDEAL 
  run_config->gamma = 5.0/3.0;  // Default Values 
  #endif
  
  #if RADIATION
  run_config->rad_constant           = 1.0;
  run_config->rad_ideal_gas_constant = 1.0;
  run_config->rad_absorption_coeff   = 0.0;
  run_config->rad_scattering_coeff   = 0.0;
  run_config->rad_c                  = 1.0;
  run_config->rad_reduced_c          = 1.0;
  run_config->rad_max_iter           = 100;
  run_config->rad_tolerance          = 1e-7;
  run_config->small_enr           = 1e-16;
  #endif

  run_config->time = 0.0;
  for (nv = 0; nv < USER_DEF_PARAMETERS; nv++) {
    run_config->inputParam[nv] = runtime->aux[nv];
  }

  for (int dir = 0; dir < 3; dir++) {
    run_config->domBeg[dir] = grid->xbeg_glob[dir];
    run_config->domEnd[dir] = grid->xend_glob[dir];
    run_config->x[dir]      = grid->xC[dir];
    run_config->dx[dir]     = grid->dx[dir];
  }
  run_config->Nx1 = grid->nint_glob[IDIR];
  run_config->Nx2 = grid->nint_glob[JDIR];
  run_config->Nx3 = grid->nint_glob[KDIR];

  // ------------------------------------------------------
  // 8. Set pointers to other structures;
  //    initialize sweep structure;
  //    initialize output
  // ------------------------------------------------------

  data->run_config = run_config;
  data->output     = runtime->output;
  data->Dts        = Dts;

  // Need to use run_config for userdef metrics 
  #if PHYSICS==GRMHD
  SetGrid(runtime, procs, grid, shrink, *run_config);
  #endif

  SetOutput (data, runtime, grid);
  #if (PARTICLES != NO)
  Particles_SetOutput(data, runtime);
  #endif

  data->fluidRiemannSolver = SetSolver (runtime->solv_type);
  #if RADIATION
  data->radiationRiemannSolver = RadSetSolver (runtime->rad_solv_type);
  #endif

  // ------------------------------------------------------
  // 9. Assign initial conditions
  // ------------------------------------------------------

  Startup (data, grid);

  // Compute initial irradiation flux divergence (set by the user in irradiation.cpp).
  #if IRRADIATION
  RadIrradiationFlux(data, grid);
  #endif

  // #pragma acc enter data copyin(run_config)
  // #pragma acc enter data copyin(grid)
  #pragma acc enter data copyin(Dts)
  #pragma omp target enter data map(to:Dts[:1])
  // #pragma acc enter data copyin(data)
  //#pragma omp target enter data map(to:grid[:1], Dts[:1], data[:1], run_config[:1])

  #if (PARTICLES != NO)
    if (cmd_line->prestart == NO) Particles_Set(data, grid); 
  #endif

  // ------------------------------------------------------
  // 10. Print normalization units
  // ------------------------------------------------------

  ShowUnits(*run_config);  

  // --------------------------------------------
  // 11. Enable writing of grid information to be 
  //     read by Gnuplot. This is done for 
  //      educational purposes only.
  // --------------------------------------------

  #ifdef GNUPLOT_HEADER
  GnuplotSetting(runtime, grid);
  #endif

  CheckConfig();
}

#ifdef PARALLEL
////////////////////////////////////////////////////////////////////////
/*!
  Returns the parallel domain decomposition mode.
  
  \param [in]  cmd_line  pointer to the cmdLine structure
  \param [out] procs     an array of integers giving the number
  of processors in each direction only if
  the -dec command line option has been given
  
  \return  The decomposition mode:
  
  - AL_AUTO_DECOMP   defaults
  - AL_USER_DECOMP   if the -dec n1 [n2] [n3] command line
  argument has been given;
  In this case only, procs[] contains the number
  of processors in the three directions;
  -  AL_MPI_DECOMP   if [todo]
*/
////////////////////////////////////////////////////////////////////////
int GetDecompMode (cmdLine *cmd_line, int procs[])
{
  int  nprocs;
  long int npx = 1, npy = 1, npz = 1;

  DIM_EXPAND(npx = cmd_line->nproc[IDIR];  ,
             npy = cmd_line->nproc[JDIR];  ,
             npz = cmd_line->nproc[KDIR];)

  // 1. if -dec has not been given
  //  decomposition will be set to be AL_AUTO_DECOMP

  if (npx == -1 || npy == -1 || npz == -1){
    #if USE_AL == YES
    return AL_AUTO_DECOMP;
    #else
    return AUTO_DECOMP;
    #endif
  }

  // 2. Enter in user decomp mode if the number of
  //  processors along all directions has been given

  if (npx > 0 && npy > 0 && npz > 0){
    procs[IDIR] = npx;
    procs[JDIR] = npy;
    procs[KDIR] = npz;

  // Check if decomposition is correct 

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (procs[IDIR]*procs[JDIR]*procs[KDIR] != nprocs){
      printf ("! The specified parallel decomposition (%d, %d, %d) is not\n",
               procs[IDIR],procs[JDIR],procs[KDIR]);
      printf ("! consistent with the number of processors (%d).\n",nprocs);
      QuitPLUTO(" ", 1);
    }
    #if USE_AL == YES
    return AL_USER_DECOMP;
    #else
    return USER_DECOMP;
    #endif
  }

  QuitPLUTO("! GetDecompMode: invalid decomposition mode",1);
  return -1;
}
#endif

////////////////////////////////////////////////////////////////////////
void CheckConfig(void)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if (PHYSICS == MHD) && (BACKGROUND_FIELD  == YES) 
  #if DIVB_CONTROL == EIGHT_WAVES
  QuitPLUTO("! CheckConfig(): Background field works with CT or GLM ONLY",1);
  #endif

  #if (DIVB_CONTROL == DIV_CLEANING) && (GLM_EXTENDED == YES)
  QuitPLUTO("! CheckConfig(): Background field not supported with EGLM",1);
  #endif

  #if CHAR_LIMITING == YES
  QuitPLUTO("! CheckConfig(): Background field does not support characteristic limiting",1);
  #endif
  
  #endif
}

#ifdef PARALLEL
#if USE_AL == NO
void DecompDomain(MPI_Comm base_comm, MPI_Datatype dtype,
                  int ndim, int *gsize, int *ghosts,
                  int *periods, int *pardim, int *procs,
                  int decomp_mode, DomainDecomp *sz)
{
  int i;
  int nprocs;
  MPI_Comm_size(base_comm, &nprocs);

  // 0. Store basic info
  sz->comm  = base_comm;
  sz->type  = dtype;
  sz->ndim  = ndim;
  sz->size  = nprocs;
  for (i = 0; i < ndim; i++) {
    sz->gsize[i]      = gsize[i];
    sz->ghosts[i]     = ghosts[i];
    sz->periods[i]    = periods[i];
    sz->isparallel[i] = pardim[i];
    sz->isstaggered[i] = 0;
    sz->cart_dims[i]   = 0;  // will be filled by MPI_Dims_create
  }

  // 1. Determine number of processes per direction --

  if (decomp_mode == AUTO_DECOMP) {
    // 1a. Let MPI decide the base decomposition
    for (i = 0; i < ndim; i++) {
      sz->cart_dims[i] = (pardim[i] ? 0 : 1);
    }
    MPI_Dims_create(nprocs, ndim, sz->cart_dims);

    // 1b. ALIGNMENT LOGIC: Match largest process counts to largest grid dimensions
    int p_dims[MAX_DIM], g_dims[MAX_DIM];
    int p_idx[MAX_DIM], g_idx[MAX_DIM];
    
    // Copy current dimensions and initialize indices
    for (i = 0; i < ndim; i++) {
      p_dims[i] = sz->cart_dims[i];
      p_idx[i]  = i;
      g_dims[i] = gsize[i];
      g_idx[i]  = i;
    }

    // Bubble sort processor counts and grid sizes in descending order
    for (int a = 0; a < ndim; a++) {
      for (int b = a + 1; b < ndim; b++) {
        if (p_dims[a] < p_dims[b]) {
          int tmp = p_dims[a]; p_dims[a] = p_dims[b]; p_dims[b] = tmp;
          tmp = p_idx[a]; p_idx[a] = p_idx[b]; p_idx[b] = tmp;
        }
        if (g_dims[a] < g_dims[b]) {
          int tmp = g_dims[a]; g_dims[a] = g_dims[b]; g_dims[b] = tmp;
          tmp = g_idx[a]; g_idx[a] = g_idx[b]; g_idx[b] = tmp;
        }
      }
    }

    // Map the sorted process dimensions to the sorted grid dimensions
    for (i = 0; i < ndim; i++) {
      sz->cart_dims[ g_idx[i] ] = p_dims[i];
    }
  } else {
    // USER_DECOMP: procs[] has been set by the user via -dec
    for (i = 0; i < ndim; i++) {
      sz->cart_dims[i] = procs[i];
    }
  }

  // 2. Create Cartesian communicator --

  MPI_Cart_create(base_comm, ndim, sz->cart_dims, periods,
                  /*reorder=*/1, &sz->cart_comm);

  // 3. Get this rank's coordinates in the Cartesian grid --

  MPI_Comm_rank(sz->cart_comm, &sz->rank);
  MPI_Cart_coords(sz->cart_comm, sz->rank, ndim, sz->coords);

  // 4. Find left/right neighbors in each direction --

  for (i = 0; i < ndim; i++) {
    MPI_Cart_shift(sz->cart_comm, i, 1, &sz->left[i], &sz->right[i]);
  }

  // 5. Compute local domain size and bounds per direction --

  for (i = 0; i < ndim; i++) {
    int np   = sz->cart_dims[i];
    int ci   = sz->coords[i];
    int N    = gsize[i];
    int ng   = ghosts[i];

    int base  = N / np;
    int extra = N % np;

    int local_n = base + (ci < extra ? 1 : 0);
    int offset  = ci * base + (ci < extra ? ci : extra);

    sz->lsize[i] = local_n;

    // global position of this rank's data in the global array
    sz->beg[i]  = ng + offset;
    sz->end[i]  = ng + offset + local_n - 1;

    // local indices within this rank's array (always starts at ng)
    sz->lbeg[i] = ng;
    sz->lend[i] = ng + local_n - 1;

    // global extent of the full domain (same for all ranks)
    sz->gbeg[i] = ng;
    sz->gend[i] = ng + N - 1;

    sz->is_gbeg[i] = (ci == 0);
    sz->is_gend[i] = (ci == np - 1);
  }

  // 6. Create 1D sub-communicators along each direction

  for (i = 0; i < ndim; i++) {
    int remain[3] = {0, 0, 0};
    remain[i] = 1;
    MPI_Cart_sub(sz->cart_comm, remain, &sz->oned_comm[i]);
  }
}
#endif
#endif
