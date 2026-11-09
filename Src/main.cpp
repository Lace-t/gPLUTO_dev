////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief PLUTO main function.

  It contains the PLUTO main function and several other top-level routines.
  main() provides basic code initialization, handles the the principal
  integration loop and calls the output driver write_data.cpp.
  Other useful functions contained in this file are NextTimeStep() responsible for computing the
  next time step based on the information available at the last time
  level.

  The standard integration loop consists of the following steps:

   - Check for last step & adjust dt if necessary
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Check output/analysis:  t(n) < tout < t(n)+dt(n)
   - write to disk/call analysis using {U(n), t(n), dt(n)}
   - Advance solution using dt(n): U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Show dominant time step (n)
   - Get next time step dt(n+1)
   - [MPI] reduction operations (n)
   - Increment n --> n+1

  \author  A. Mignone  (andrea.mignone@unito.it)
           M. Rossazza (marco.rossazza@unito.it)
           V. Berta    (vittoria.berta@unito.it)

  \date   Nov 23, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "time.h"
#if defined(_OPENACC)
#include <cuda_runtime.h>
#endif

#ifndef DATAMINMAX
#define DATAMINMAX  YES
#endif

static void CheckForOutput (Data *, timeStep *, Runtime *, time_t t0, Grid *);
static void CheckForAnalysis (Data *, timeStep *, Runtime *, Grid *);
static void DataMinMax (Data *, Grid *);
static string TotalExecutionTime(double);

int main(int argc, char *argv[])
{
  auto clock_beg_whole = chrono::high_resolution_clock::now();

  int       i, j, k;
  int       nstep = 0, last_step = 0;
  Data      data;
  Grid      grid;
  timeStep  Dts;
  cmdLine   cmd_line;
  Runtime   runtime;
  RunConfig run_config;

  #ifdef PARALLEL
  int      nRanks;
  #pragma omp target enter data map(always,alloc:nccl_comm)
  #endif
  #if PARTICLES != NO
  pcSet(data);
  #endif
  char first_step = 1;
  RBox box_tot;

  #pragma omp target enter data map(always,alloc:g_usedMemory)
  // #pragma omp target enter data map(always,alloc:g_useGPU)

  // ------------------------------------------------------
  // 0. Initialize parallel environment
  // ------------------------------------------------------

#ifdef PARALLEL
#if USE_AL == NO
MPI_Init(&argc, &argv);
#else
  AL_Init (&argc, &argv);
#endif
  MPI_Comm_rank (MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  // 0a. Set the device for each rank
  
  #if defined(_OPENACC) || defined(_OPENMP)
  ncclUniqueId id;

  // Device num is the local rank mod number 
  // of devices on the node
  int device_type;
  int num_devices;
  int gpuId;

  #if defined(_OPENACC)
  
  device_type = acc_get_device_type();
  num_devices = acc_get_num_devices(acc_get_device_type());
  gpuId = g_rank % num_devices;

  acc_set_device_num(gpuId, acc_get_device_type());
  acc_init(acc_get_device_type());

  #elif defined(_OPENMP)

  device_type = omp_get_default_device();
  num_devices = omp_get_num_devices();
  gpuId = g_rank % num_devices;
  omp_set_default_device(gpuId);

  #endif

  printf ("[Rank %d]\n", g_rank);
  printf ("  device number = %d\n", gpuId);

  if (num_devices != 1) g_nccl = YES;

  if (g_nccl) {
  // get NCCL unique ID at rank 0 and broadcast it to all others
  if (g_rank == 0) {
    ncclResult_t res = ncclGetUniqueId(&id);
    if(res != ncclSuccess) {
      printf("[%d] res = %d (ncclGetUniqueId)\n",g_rank, res);
    }
  }
  MPI_Bcast((void *) &id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD);

  // initializing NCCL
  ncclResult_t res =  ncclCommInitRank(&nccl_comm, nRanks, id, g_rank);
  if(res != ncclSuccess) {
    printf("[%d] res = %d (ncclCommInitRank)\n",res, g_rank);
  }

  } // Multiple ranks may share one GPU; MPI handles their exchanges.
  #endif // _OPENACC and _OPENMP 
  
#endif // PARALLEL 

  // ----------------------------------------------------------
  // 1. Parse command line, read parameter file, 
  //     set default structure members, initialize grid.
  // ----------------------------------------------------------

  ParseCmdLineArgs (argc, argv, &cmd_line);
  if (g_rank == 0) RuntimeSetup (&runtime, &cmd_line, cmd_line.input_file);

  #ifdef PARALLEL
  MPI_Bcast (&runtime,  sizeof (Runtime), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  RuntimeSet(&runtime);

  // --------------------------------------------------------
  // 2. Initialize members of timeStep structure
  // -------------------------------------------------------- 

  #if COOLING != NO
  Dts.dt_cool   = 1.e38;
  #endif 
  #if PARTICLES != NO
  #if PARTICLES == PARTICLES_CR
  Dts.Nsub_particles  = MAX(1, -PARTICLES_CR_NSUB);
  #else
  Dts.Nsub_particles = 1;
  #endif
  Dts.invDt_particles  = 1.e-38;
  Dts.omega_particles  = 1.e-38;
  Dts.particles_tstart = runtime.particles_tstart;
  #endif

  // ------------------------------------------------------
  // 3. Open log file for the first time
  // ----------------------------------------------

  if (cmd_line.restart == NO && cmd_line.h5restart == NO){
    LogFileOpen (runtime.log_dir, "w");
  }else{
    LogFileOpen (runtime.log_dir, "a");
  }

  // -------------------------------------------------
  // 4. Initialize parallel environment, grid,
  //     memory allocation.
  // -------------------------------------------------

  ShowConfig (cmd_line.input_file, run_config);
  Initialize(&data, &Dts, &cmd_line, &runtime, &run_config, &grid);
  
  #if COOLING != NO
    #if COOLING ==  TABULATED
      // read tabulated cooling file outside the ACC zone
      ReadTabFile(run_config);
    #endif
  #endif
  
  Dts.time = 0.;
  time_t tbeg = time(0);
  g_stepNumber = nstep = 0;
  // #pragma acc update device(g_stepNumber)
  // #pragma omp target update to(g_stepNumber)


  
	#if 0
  // -------------------------------------------------------
  // 5a. The array where the solution is going to be stored is 
  // &ta.poisson_sol. 
  //     This array has the full size of the Grid, 
  //     but the points that are going to be modified are
  //     the local grid points, plus one additional for the 
  //     Poisson boundary condition. 
  //
  //     i.e., from IBEG to IEND is the active domain, 
  //     using IBEG-1 and IEND+1 as b.c. for the Poisson equation: 
  //                          ∇φ = 4πS. 
  //
  //     SingleGridPoissonSolver or MultiGridPoissonSolver takes 
  //     3 arguments:
  //     1. Ary3D storing the source term S for the Poisson equation;
  //     2. Ary3D storing the solution φ of the Poisson equation;
  //     3. The Grid structure. 
  //
  //     Example of a SELF_GRAV correction routine
  // ------------------------------------------------------

    #if (SELF_GRAV == YES)
      Ary3D density;
      CopyInData(density, &data, &grid, nstep);
      #if (POISSON_SOLVER_METHOD == MULTIGRID)
        MultiGridPoissonSolver(density, data.poisson_sol, &grid); 
      #else
        SingleGridPoissonSolver(density, data.poisson_sol, &grid);
      #endif
      CalculateGravAcceleration(&data, &grid, data.poisson_sol);
    #endif

    //  Example of a Divergence correction routine 
    #if (DIVERGENCE_FREE == YES)
      #if (POISSON_SOLVER_METHOD == MULTIGRID)
        MultiGridPoissonSolver(data.divB, data.poisson_sol, &grid); 
      #else
        SingleGridPoissonSolver(data.divB, data.poisson_sol, &grid);
      #endif
      DivergenceFreeCorrect(&data, &grid);
    #endif

    //  Test Problems for PoissonSolver 
    #if (POISSON_TEST == YES)
      #if (POISSON_SOLVER_METHOD == MULTIGRID)
        MultiGridPoissonSolver(data.divB, data.poisson_sol, &grid); 
      #else
        SingleGridPoissonSolver(data.divB, data.poisson_sol, &grid);
      #endif

      // PoissonTestErrorCheck(data.poisson_sol, &grid);
      // PoissonTestWrite(data.poisson_sol, &grid, g_rank);
    #endif
  #endif

  // -------------------------------------------------
  // 0g. Check if restart is necessary.
  //     If not, write initial condition to disk.
  //     Copy structures on GPU.
  // -------------------------------------------------

  if (cmd_line.restart == YES) {
    RestartFromFile (&data, &runtime, &Dts, cmd_line.nrestart, DBL_OUTPUT, &grid);
    #if PARTICLES
    Particles_Restart(&runtime, &data, &Dts, cmd_line.nrestart, &grid);
    #endif
    #pragma acc update device(data)
    #pragma omp target update to(data)
  }else if (cmd_line.h5restart == YES){
    RestartFromFile (&data, &runtime, &Dts, cmd_line.nrestart, DBL_H5_OUTPUT, &grid);
    #pragma acc update device(data)
    #pragma omp target update to(data)
  }
  
  if (cmd_line.write && (cmd_line.restart == NO && cmd_line.h5restart == NO)){
    CheckForAnalysis (&data, &Dts, &runtime, &grid);
    CheckForOutput   (&data, &Dts, &runtime, tbeg, &grid);
  }
  nstep = g_stepNumber;

  RBoxDefine(0, grid.ntot[IDIR] - 1,  
             0, grid.ntot[JDIR] - 1,  
             0, grid.ntot[KDIR] - 1,  &box_tot);

  #pragma acc enter data copyin(box_tot)
  #pragma omp target enter data map(to: box_tot)
  //#pragma omp target enter data map(to: data[:1] , grid[:1], box_tot[:1], Dts[:1])

  #if VLBI_PLUTO == YES
    // user parameters:  angle ? Freq ? file beg/end ? 
    VLBI_Pluto(data, runtime, grid);
    QuitPLUTO(1);
  #endif

  // ===================================================================
  // 6.      M A I N     L O O P     S T A R T S     H E R E
  // ===================================================================
  
  #pragma acc update device(Dts,data) 
  print ("\n> Starting computation\n");
  LogFileFlush();
  auto clock_beg = chrono::high_resolution_clock::now();

  mynvtxstart_("MAIN_WHILE",SILVER);

  while (last_step!=1){


		// --------------------------------------------------------
  	// 7. Poisson Solver
  	// -------------------------------------------------------
		#if (POISSON_SOLVER)
			if(g_stepNumber % POISSON_SOLVER_FREQ == 0)
			PoissonMain(data,grid);
		#endif 
		

    #if COOLING != NO
    Dts.dt_cool   = 1.e38;   
    #endif

    g_stepNumber = Dts.stepNumber = nstep;
    // #pragma acc update device(g_stepNumber)
    // #pragma omp target update to(g_stepNumber)

    // ----------------------------------------------------
    // 6a. Example of a SELF_GRAV correction routine:
    //     Solve the gravitational potential at 
    //     the beginning of every time step.
    // ----------------------------------------------------

    #if (SELF_GRAV == YES)
      CopyInData(density, &data, &grid, nstep+1);
      #if (POISSON_SOLVER_METHOD == MULTIGRID)
        MultiGridPoissonSolver(density, data.poisson_sol, &grid); 
      #else
        SingleGridPoissonSolver(density, data.poisson_sol, &grid);
      #endif
      CalculateGravAcceleration(&data, &grid, data.poisson_sol);
    #endif

    // ----------------------------------------------------
    // 6b. Check if this is the last integration step:
    //     - final tstop has been reached: adjust time step
    //     - or max number of steps has been reached
    // ----------------------------------------------------

    if ((Dts.time + Dts.dt) >= runtime.tstop*(1.0 - 1.e-8)) {
      Dts.dt    = (runtime.tstop - Dts.time);
      // #pragma acc update device(Dts)
      // #pragma omp target update to(Dts)
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps >= 0) {
      last_step = 1;
    }

    // ----------------------------------------------------
    // 6c. Check if it's time to write or perform analysis
    // ----------------------------------------------------

    #if defined(_OPENACC)
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    if(static_cast<double>(free_mem) < 0.08*static_cast<double>(total_mem)){
     QuitPLUTO("Memory Sturated!\n", 1);
    }
    #endif

    if (!first_step && cmd_line.write) {
      CheckForAnalysis (&data, &Dts, &runtime, &grid);
      if (!last_step) CheckForOutput  (&data, &Dts, &runtime, tbeg, &grid);
    }

    if (g_stepNumber%runtime.log_freq == 0) {
      print ("> nstep = %d; t = %0e; dt = %0e; %6.1f %%\n",
                nstep, Dts.time, Dts.dt, 100.0*Dts.time/runtime.tstop);
 
      #if PARTICLES != NO
      {
        uint64_t nptot = data.pc->Size();
        #ifdef PARALLEL
        MPI_Allreduce (MPI_IN_PLACE, &nptot, 1, MPI_UNSIGNED_LONG_LONG, 
                       MPI_SUM, MPI_COMM_WORLD);
        #endif
        print("> nParticles = %llu\n", nptot);
      }
      #endif
    }

    // ----------------------------------------------------
    // 6d. Reset flag
    // ----------------------------------------------------
 
    #pragma omp target teams distribute parallel for collapse(3) 
    #pragma acc parallel loop collapse(3) present(data, box_tot)
    BOX_LOOP(box_tot,i,j,k) {
      data.flag(i,j,k) = 0;
    }
    data.err_count        = 0;
    data.integration_mode = DEFAULT;
    #pragma acc update device (data.integration_mode)
    #pragma omp target update to(data.integration_mode)
    
    // ----------------------------------------------------
    // 6e. Advance solution array by a single time step
    //     Dts.dt. After this step U^n -> U^{n+1}
    // ----------------------------------------------------

    #ifdef GLM_MHD  // initialize glm_ch
    GLM_Init (&data, &grid);
    #endif

    if (cmd_line.jet != -1) SetJetDomain (&data, cmd_line.jet, 
                                          runtime.log_freq, &grid);

    #ifdef FARGO
    FARGO_AverageVelocity(&data, &grid);
    #endif

    #if FALLBACK_MODE != FALLBACK_RECOMPUTE
    SplitSource (&data, 0.5*(Dts.dt), &Dts, &grid);
    AdvanceStep (&data, &grid, &Dts);
    SplitSource (&data, 0.5*(Dts.dt), &Dts, &grid);
    #else
    Fallback(&data, &grid, &Dts);
    #endif

    if (cmd_line.jet != -1) UnsetJetDomain (cmd_line.jet, &grid);

    // ----------------------------------------------------
    // 6f. Increment time & compute next time step
    // ----------------------------------------------------
    
    Dts.time             += Dts.dt;
    data.run_config->time = Dts.time;
    #pragma acc update device(data.run_config->time)
    #pragma omp target update to(data.run_config->time)
    
    Dts.dt = NextTimeStep(&Dts, &run_config, &runtime, &grid);
    // #pragma acc update device(Dts)
    // #pragma omp target update to(Dts)

    #ifdef GLM_MHD
    if (g_stepNumber%runtime.log_freq == 0) {
      print ("  glm_ch = %f\n",run_config.glm_ch);
    }
    #endif

    // ----------------------------------------------------
    // 6g. Compute Max Mach number & reduce var,
    //     Reduce variables through min/max 
    // ----------------------------------------------------
   
    if (g_stepNumber%runtime.log_freq == 0) {
      DataMinMax(&data, &grid);
    }

    nstep++;
    first_step = 0;
    LogFileFlush();
  }

  // ===================================================================
  //          M A I N       L O O P      E N D S       H E R E
  // ===================================================================

  mynvtxstop_();

  auto clock_end = std::chrono::high_resolution_clock::now();

  // 7a. Check for output

  if ((cmd_line.write) && !(cmd_line.maxsteps == 0)) {
    CheckForOutput (&data, &Dts, &runtime, tbeg, &grid);
    CheckForAnalysis (&data, &Dts, &runtime, &grid);
  }

  // 7b. Print diagnostic

  std::chrono::duration<double> clock_tot = clock_end - clock_beg;
  double clock_diff = clock_tot.count();
  string tel = TotalExecutionTime(clock_diff);

  #ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);

    print(
      "\n> Total allocated memory (fluid) %6.2f Mb (proc #%d)\n",
      (float)g_usedMemory / 1.e6, g_rank
    );

    #if PARTICLES != NO
      print(
        "\n> Particles allocated memory %6.2f Mb (proc #%d)\n",
        data.pc->TotalSize() / 1.e6, g_rank
      );
    #endif

    MPI_Barrier(MPI_COMM_WORLD);
  #else
    print(
      "\n> Total allocated memory (fluid) %6.2f Mb\n",
      (float)g_usedMemory / 1.204e6
    );

    #if PARTICLES != NO
      print(
        "> Particles allocated memory     %6.2f Gb\n",
        static_cast<double>(data.pc->TotalSize()) / 1.204e9
      );
    #endif
  #endif

  if (clock_diff < 1.e3) {
    print ("> Elapsed time (compute)            %8.4f (sec) [%s]\n", 
            clock_diff, tel.c_str());
  }else{
    print ("> Elapsed time (compute)            %10.4e (sec) [%s]\n", 
            clock_diff, tel.c_str());
  }
  print ("> Average time/step        %10.2e (sec)\n", clock_diff/nstep);

  // 7c. Close file, finalize environment 

  #pragma acc exit data delete(data)
  #pragma acc exit data delete(Dts)
  #pragma acc exit data delete(grid)
  #pragma acc exit data delete(box_tot)
  #pragma omp target exit data map(delete: data, Dts, grid, box_tot)

  auto clock_end_whole = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> clock_tot_whole = clock_end_whole - clock_beg_whole;
  double clock_diff_whole = clock_tot_whole.count();
  string tel_whole = TotalExecutionTime(clock_diff_whole);

  if (clock_diff_whole < 1.e3) {
    print (
      "> Elapsed time (total)            %8.4f (sec) [%s]\n", 
            clock_diff_whole, tel_whole.c_str()
          );
  }else{
    print (
      "> Elapsed time (total)            %10.4e (sec) [%s]\n", 
            clock_diff_whole, tel_whole.c_str()
          );
  }
  print ("> Done.\n");

  // For density array for SELF_GRAV 
 
 #if (SELF_GRAV == YES)
    density.destroy();
  #endif  
 
  #ifdef PARALLEL
  LogFileClose();
  #ifdef _OPENACC
  if (g_nccl) ncclCommDestroy(nccl_comm);
  #endif
  MPI_Barrier (MPI_COMM_WORLD);
#if USE_AL == NO
MPI_Finalize();
#else
  AL_Finalize ();
#endif  
  #endif

  return 0;
}

////////////////////////////////////////////////////////////////////////
/*!
  Check if file output has to be performed.
*/
////////////////////////////////////////////////////////////////////////
void CheckForOutput (Data *data, timeStep *Dts, Runtime *runtime,
                     time_t t0, Grid *grid)
{
  static int first_call = 1;
  int    n, check_dt, check_dn, check_dclock;
  int    first_step     = (g_stepNumber == 0 ? 1:0);
  int    last_step      = (Dts->time >= runtime->tstop);
  int    restart_update = 0;
  double t  = Dts->time;
  double dt = Dts->dt;
  double dtime[MAX_OUTPUT_TYPES];
  static time_t tbeg[MAX_OUTPUT_TYPES], tend;
  Output *output;

  // 0. Initialize

  if (first_call) for (n = 0; n < MAX_OUTPUT_TYPES; n++) tbeg[n] = t0;

  time (&tend);

  for (n = 0; n < MAX_OUTPUT_TYPES; n++) dtime[n] = difftime(tend, tbeg[n]);
  #ifdef PARALLEL
  MPI_Bcast(dtime, MAX_OUTPUT_TYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif

  // ------------------------------------------------------
  // 1. Start main loop on outputs
  // ------------------------------------------------------

  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    output = runtime->output + n;
    check_dt = check_dn = check_dclock = 0;

    // 1a. Enable writing if  t(n) < t(write) < t(n+1)

    if (output->dt > 0.0){
      check_dt = (int) ((t + dt)/output->dt) - (int)(t/output->dt);
      check_dt = check_dt || first_step || last_step;
    }

    // 1b. Enable writing if current step is an integer
    //     multiple of output interval

    if (output->dn > 0){
      check_dn = (g_stepNumber%output->dn) == 0;
      check_dn = check_dn || first_step || last_step;
    }

    // 1c. Enable writing if clock time has been reached

    if (output->dclock > 0.0){
      if (dtime[n] >= output->dclock) {
        check_dclock = 1;
        tbeg[n] = tend;
      }else{
        check_dclock = 0;
      }
      check_dclock = check_dclock || first_step || last_step;
    }

    // 1d. If any of the previous is true dump data to disk 
    if (check_dt || check_dn || check_dclock) {

      #if PARTICLES
      if (  output->type == PARTICLES_DBL_OUTPUT
         || output->type == PARTICLES_FLT_OUTPUT
         || output->type == PARTICLES_VTK_OUTPUT
         || output->type == PARTICLES_TAB_OUTPUT
         || output->type == PARTICLES_HDF5_OUTPUT){
        Particles_WriteData(data, Dts, output, grid);
      } else
      #endif

      WriteData(data, Dts, output, grid);

      // 1e. Save the file number of the dbl and dbl.h5 output format
      //     for writing restart.out once we exit the loop.

      if ((output->type == DBL_OUTPUT) ||
          (output->type == DBL_H5_OUTPUT)) restart_update = 1;
    }
  }

  // 2. Dump restart information if required.
  //    Note that if both dbl and dbl.h5 formats are used,
  //    bookeeping is done using dbl format.

  if (restart_update) RestartDump (runtime, data->restartParam, Dts);

  first_call = 0;
}

////////////////////////////////////////////////////////////////////////
/*!
  Check if Analysis needs to be called
*/
////////////////////////////////////////////////////////////////////////
void CheckForAnalysis (Data *data, timeStep *Dts, Runtime *runtime, Grid *grid)
{
  int check_dt = 0, check_dn = 0;
  double t, tnext;

  t     = Dts->time;
  tnext = t + Dts->dt;
  check_dt = (int) (tnext/runtime->anl_dt) - (int)(t/runtime->anl_dt);
  check_dt = check_dt || g_stepNumber == 0 || fabs(t - runtime->tstop) < 1.e-9;
  check_dt = check_dt && (runtime->anl_dt > 0.0);

  check_dn = (g_stepNumber%runtime->anl_dn) == 0;
  check_dn = check_dn && (runtime->anl_dn > 0);

  if (check_dt || check_dn) Analysis (data, grid);

  #if PARTICLES != NO
  check_dt = (int) (tnext/runtime->Panl_dt) - (int)(t/runtime->Panl_dt);
  #endif
  check_dt = check_dt || g_stepNumber == 0 || fabs(t - runtime->tstop) < 1.e-9;

  #if PARTICLES != NO
  check_dt = check_dt && (runtime->Panl_dt > 0.0);
  check_dn = (g_stepNumber%runtime->Panl_dn) == 0;
  check_dn = check_dn && (runtime->Panl_dn > 0);
  if (check_dt || check_dn) ParticlesAnalysis (data, grid);
  #endif // PARTICLES != NO 
}

////////////////////////////////////////////////////////////////////////
/*!
  Converta floating-point variable (dt, in seconds) to a string
  displaying days:hours:minutes:seconds
*/
////////////////////////////////////////////////////////////////////////
string TotalExecutionTime(double dt) 
{
    int days = static_cast<int>(dt / 86400.0);
    int hours = static_cast<int>((dt - 86400.0 * days) / 3600.0);
    int mins = static_cast<int>((dt - 86400.0 * days - 3600.0 * hours) / 60.0);
    int secs = static_cast<int>(dt - 86400.0 * days - 3600.0 * hours - 60.0 * mins);

    return std::to_string(days) + "d:" +
           std::to_string(hours) + "h:" +
           std::to_string(mins) + "m:" +
           std::to_string(secs) + "s";
}


////////////////////////////////////////////////////////////////////////
/*!
  Print max & min value of variabels
*/
////////////////////////////////////////////////////////////////////////
void DataMinMax (Data *data, Grid *grid)
{
  int i,j,k,nv;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  double cs, mach_max;
  double vc_max[NVAR], vc_min[NVAR];
  #if !HAVE_ENERGY
  RunConfig &run_config = *(data->run_config);
  double iso_sound_speed = run_config.isoSoundSpeed;
  #endif

  // ------------------------------------------------------
  // 1. Print total number of errors 
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0
  if (data->err_count > 0) {
    print ("  ! %llu errors detected [proc #%d]\n", data->err_count, g_rank);
  }
  #endif

#if DEBUG_MIN_MAX == YES

  mach_max = 0.0;

  #pragma omp target teams distribute parallel for collapse(3) reduction(max: mach_max)
  #pragma acc parallel loop collapse(3) reduction(max:mach_max) \
              present(grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    #if HAVE_ENERGY
    cs = data->Vc(i,j,k,PRS) / data->Vc(i,j,k,RHO);
    cs = sqrt(cs);
    #else
    cs = iso_sound_speed;
    #endif

    double mach_IDIR = fabs(data->Vc(i,j,k,VX1)) / cs;
    double mach_JDIR = fabs(data->Vc(i,j,k,VX2)) / cs;
    double mach_KDIR = fabs(data->Vc(i,j,k,VX3)) / cs;
    
    // Compute local maximum
    double local_max = MAX(mach_IDIR, mach_JDIR);
    local_max = MAX(local_max, mach_KDIR);

    // Single reduction update
    mach_max = MAX(mach_max, local_max);
  }}}

  #ifdef PARALLEL
  MPI_Allreduce (MPI_IN_PLACE, &mach_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  #endif
  print ("  Max(MACH) =  %f\n",mach_max);

  if(isnan(mach_max)){
    QuitPLUTO("! Main(): mach_max = nan",1);
  }

  // ------------------------------------------------------
  // 2. Take max / min of primitive variables and reduce them
  // ------------------------------------------------------

  // 1. Initialize arrays on the host before the parallel region
  NVAR_LOOP(nv) {
      vc_max[nv] = -1e18;
      vc_min[nv] = +1e18;
  }

  // 2. Perform a single parallel launch with array reductions
  #pragma omp target teams distribute parallel for collapse(3) \
                reduction(max: vc_max[:NVAR]) reduction(min: vc_min[:NVAR])
  #pragma acc parallel loop collapse(3) present(grid, data) \
              reduction(max: vc_max[:NVAR]) reduction(min: vc_min[:NVAR])
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    
    // 3. Evaluate all variables within the inner loop
    NVAR_LOOP(nv){
      double val = data->Vc(i,j,k,nv);
      
      vc_max[nv] = MAX(vc_max[nv], val);
      vc_min[nv] = MIN(vc_min[nv], val);
    }

  }}}

  #ifdef PARALLEL
  MPI_Allreduce (MPI_IN_PLACE, vc_min, NVAR, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, vc_max, NVAR, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  #endif

  // ------------------------------------------------------
  // 3. Print to log
  // ------------------------------------------------------


/*  Alternative construction using C++ strings ...
    **********************************************
string format_string_max="";
string format_string_min="";
#if (PHYSICS == HD) || (PHYSICS == RHD)
intList var_list = {5, {RHO, VX1, VX2, VX3, PRS}};
#elif (PHYSICS == MHD) || (PHYSICS == RMHD)
  #if HAVE_ENERGY
  intList var_list = {8, {RHO, VX1, VX2, VX3, BX1, BX2, BX3, PRS}};
  #else
  intList var_list = {7, {RHO, VX1, VX2, VX3, BX1, BX2, BX3}};
  #endif
  #ifdef GLM_MHD
  var_list.nvar++;
  var_list.
  #endif

#endif

FOR_EACH(nv, &var_list){
  char val[32];

  sprintf (val,"% 0e",vc_max[nv]);
  format_string_max.append(string(val)+"  ");

  sprintf (val,"% 0e",vc_min[nv]);
  format_string_min.append(string(val)+"  ");
}
cout << "max  = " << format_string_max << endl;
cout << "min  = " << format_string_min << endl;
exit(1);
*/

  #if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == GRMHD)
  #if HAVE_ENERGY
  #ifdef GLM_MHD
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e % 0e\n",
            vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3],
                         vc_max[BX1], vc_max[BX2], vc_max[BX3],
                         vc_max[PRS], vc_max[PSI_GLM]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e % 0e\n",
            vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3],
                         vc_min[BX1], vc_min[BX2], vc_min[BX3],
                         vc_min[PRS], vc_min[PSI_GLM]);
  #else
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e\n",
            vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3],
                         vc_max[BX1], vc_max[BX2], vc_max[BX3], vc_max[PRS]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e\n",
            vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3],
                         vc_min[BX1], vc_min[BX2], vc_min[BX3], vc_min[PRS]);
  #endif  // GLM_MHD 
  #else
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e\n",
            vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3],
                         vc_max[BX1], vc_max[BX2], vc_max[BX3]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e\n",
            vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3],
                         vc_min[BX1], vc_min[BX2], vc_min[BX3]);
  #endif
  #endif


  #if PHYSICS==ResRMHD
  #if DIVB_CONTROL  == DIV_CLEANING
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e % 0e % 0e % 0e % 0e % 0e\n",
                        vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3],
                                     vc_max[BX1], vc_max[BX2], vc_max[BX3],
                                     vc_max[EX1], vc_max[EX2], vc_max[EX3],
                        vc_max[PRS], vc_max[PSI_GLM]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e % 0e % 0e % 0e % 0e % 0e % 0e  % 0e\n",
                        vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3],
                                     vc_min[BX1], vc_min[BX2], vc_min[BX3],
                                     vc_min[EX1], vc_min[EX2], vc_min[EX3],
                        vc_min[PRS], vc_min[PSI_GLM]);
  #else
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e % 0e  % 0e % 0e\n",
                        vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3],
                                     vc_max[BX1], vc_max[BX2], vc_max[BX3], 
                                     vc_max[EX1], vc_max[EX2], vc_max[EX3],
                        vc_max[PRS]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e  % 0e  % 0e  % 0e % 0e % 0e  % 0e % 0e\n",
                        vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3],
                                     vc_min[BX1], vc_min[BX2], vc_min[BX3],
                                     vc_min[EX1], vc_min[EX2], vc_min[EX3], 
                        vc_min[PRS]);
  #endif
  #endif


  #if (PHYSICS == HD) || (PHYSICS == RHD)
  #ifdef ENG
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e % 0e\n",
             vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3], vc_max[PRS]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e % 0e\n",
             vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3], vc_min[PRS]);
  #else
  print ("  Max(vc)   = % 0e  % 0e  % 0e  % 0e\n",
             vc_max[RHO], vc_max[VX1], vc_max[VX2], vc_max[VX3]);
  print ("  Min(vc)   = % 0e  % 0e  % 0e  % 0e\n",
             vc_min[RHO], vc_min[VX1], vc_min[VX2], vc_min[VX3]);
  #endif 
  #endif

  // ------------------------------------------------------
  // 4. Particle output
  // ------------------------------------------------------

  #if PARTICLES == PARTICLES_LP
  constexpr int particle_size =( PARTICLES == PARTICLES_LP) + (PARTICLES_LP_SPECTRA == YES);
  #if PARTICLES_LP_SPECTRA == YES
  char var[particle_size][12] = {"Pos", "Chi"};
  #else
  char var[particle_size][12] = {"Pos"};
  #endif
  double avg[particle_size];

  uint64_t size = data->pc->Size();
  for(i = 0; i < particle_size; i++) avg[i] = 0.0;

  #pragma acc data copy(avg[:particle_size])
  {
  #pragma acc parallel loop present(avg[:particle_size], data)
  for(uint64_t np = 0; np < size; np++){
    double value;
    for(i = 0; i < DIMENSIONS; i++) {
      value = data->pc->pos(np)[i];
      #pragma acc atomic update
      avg[0] += value;
    }
    #if  PARTICLES_LP_SPECTRA == YES
    for(i = 0; i < PARTICLES_LP_NEBINS; i++) {
        value = data->pc->chi(np)[i];
        #pragma acc atomic update
        avg[1] += value;
    }
    #endif
  }  // end particles loop 
  }


  #ifdef PARALLEL
  int nRanks;
  MPI_Allreduce (MPI_IN_PLACE,  avg, particle_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
  #endif

  for(i = 0; i < particle_size; i++) print ("  %s       = %0e\n", var[i], avg[i]);
  #endif

#endif // DEBUG_MIN_MAX 

}
