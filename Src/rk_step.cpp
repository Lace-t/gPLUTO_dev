////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Advance equations with Runge Kutta time integrators.

  Main driver for RK time-stepping (RK2, RK3 and RK4).

  \authors A. Mignone         (andrea.mignone@unito.it)\n
           M. Rossazza        (marco.rossazza@unito.it)\n
           V. Berta           (vittoria.berta@unito.it)\n
           D. Melon Fuksman   (fuksman@mpia.de)\n

  \date    May 5, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// RK-SSP Coefficients for the 1st stage

#if TIME_STEPPING != RK4
  #define W10   1.0
  #define W11   1.0
  #define W1R   1.0
#else
  #define W10   1.0
  #define W11   1.0
  #define W1R   0.391752226571890
#endif

// RK-SSP Coefficients for the 2nd stage

#if TIME_STEPPING == RK2
  #define W20   0.5
  #define W21   0.5
  #define W2R   0.5
#elif TIME_STEPPING == RK3
  #define W20   0.75
  #define W21   0.25
  #define W2R   0.25
#elif TIME_STEPPING == RK4
  #define W20   0.444370493651235
  #define W21   0.555629506348765
  #define W2R   0.368410593050371
#endif

// RK-SSP Coefficients for the 3rd stage

#if TIME_STEPPING == RK3
  #define W30   (1.0/3.0)
  #define W31   (2.0/3.0)
  #define W3R   (2.0/3.0)

  // hexadecimal floating-point literals option gives you the exact IEEE-754
  // #define W30  0x1.5555555555555p-2  
  // #define W31  0x1.5555555555555p-1
  // #define W3R  0x1.5555555555555p-1

#elif TIME_STEPPING == RK4
  #define W30   0.620101851488403
  #define W31   0.379898148511597
  #define W3R   0.251891774271694
#endif

// RK-SSP Coefficients for the 4th and 5th stages

#if TIME_STEPPING == RK4
  #define W40   0.178079954393132
  #define W43   0.821920045606868
  #define W4R   0.544974750228521

  #define W52   0.517231671970585
  #define W53   0.096059710526147
  #define W5R4  0.063692468666290
  #define W54   0.386708617503268
  #define W5R   0.226007483236906
#endif

////////////////////////////////////////////////////////////////////////
/*!
   Advance the equations by a single time step using unsplit
   integrators based on the method of lines.

   \param[in,out]   data  pointer to Data structure
   \param[in]       grid  pointer to array of Grid structures
   \param[in,out]    Dts  pointer to time step structure
*/
////////////////////////////////////////////////////////////////////////
uint32_t AdvanceStep(Data *data, Grid *grid, timeStep *Dts)
{
  mynvtxstart_("RKSTEP", GREEN);
  int nv, i, j, k;
  double dt = Dts->dt;
  #ifdef FARGO
  RunConfig &run_config = *(data->run_config);
  #endif

  if (Dts->time >= RuntimeGet()->tfreeze) {
    #if PARTICLES != NO
    Particles_Inject (data, grid);
    #if PARTICLES == PARTICLES_CR
    Particles_CR_Update (data, grid);
    #elif PARTICLES == PARTICLES_DUST
    Particles_Dust_Update (data, grid);
    #elif PARTICLES == PARTICLES_LP
    #if TIME_STEPPING == RK2
    g_intStage = 1;
    Particles_LP_Update (data, grid);

    g_intStage = 2;
    Particles_LP_Update (data, grid);

    #elif TIME_STEPPING == RK3
    g_intStage = 1;
    Particles_LP_Update (data, grid);

    g_intStage = 2;
    Particles_LP_Update (data, grid);


    g_intStage = 3;
    Particles_LP_Update (data, grid);

    #endif // TIME_STEPPING == RK2 
    #if PARTICLES_LP_SPECTRA == YES
    Boundary (data, 0, grid); // Spectra update requires interpolating
                                  // * fluid quantities at particle position. 
    Particles_LP_UpdateSpectra (data, Dts,grid);
    #endif // PARTICLES_LP_SPECTRA == YES 
    #endif // PARTICLES == PARTICLES_CR 
    #endif // PARTICLES != NO 
    mynvtxstop_();
    return 0;
  }
  
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int nx1_tot=grid->ntot[IDIR];
  int nx2_tot=grid->ntot[JDIR];
  int nx3_tot=grid->ntot[KDIR];

  RBox box, boxt;

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine(0, nx1_tot - 1, 0, nx2_tot - 1, 0, nx3_tot - 1, &boxt);
  #pragma acc enter data copyin(box, boxt)
  #pragma omp target enter data map(to: box, boxt)

  // 0a. Allocate memory 

  const  Ary4D &U0 = data->Uc0;
  static Ary4D Rc;
  static Ary3D Bs0[6], Rs[6];
  static Ary4D U1, U2, U3, U4, R4;
  static Ary3D Bs1[6], Bs2[6], Bs3[6], Bs4[6], Rs4[6];

  #if RADIATION_RELATIVISTIC
  #if TIME_STEPPING == RK4
  QuitPLUTO("! AdvanceStep(): RADIATION is not yet supported with TIME_STEPPING == RK4.\n",1);
  #endif
  #ifdef HIGH_ORDER
  QuitPLUTO("! AdvanceStep(): RADIATION is not yet supported with HIGH_ORDER.\n",1);
  #endif
  #if RADIATION_IMEX_SSP2
  #if TIME_STEPPING != RK2
  QuitPLUTO("! AdvanceStep(): RADIATION_IMEX_SSP2 requires TIME_STEPPING == RK2\n",1);
  #endif
  static Ary4D Srad1;
  static Ary4D Srad2;
  const double rk_gamma = 0.29289321881;
  const double dt_rk1 = dt*rk_gamma ;
  const double dt_rk2 = dt*(1.0-2.0*rk_gamma);
  const double dt_rk3 = dt*0.5*(1.0-rk_gamma);
  #endif
  #endif

  if (Rc.isAllocated() == 0) {
    Rc.create (nx1_tot,nx2_tot,nx3_tot,NVAR);

    #if TIME_STEPPING == RK4
    U1.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    U2.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    U3.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    U4.create (nx1_tot,nx2_tot,nx3_tot,NVAR);    
    R4.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    #endif

    #pragma acc enter data copyin(U1, U2, U3, U4, Rc, R4)
    #pragma omp target enter data map(to: U1, U2, U3, U4, Rc, R4)

    #if RADIATION_RELATIVISTIC && RADIATION_IMEX_SSP2
     Srad1.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
     Srad2.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
     #pragma acc enter data copyin(Srad1,Srad2)
     #pragma omp target enter data map(to: Srad1,Srad2)
    #endif

    #ifdef STAGGERED_MHD
    int (*ns)[3];
    ns = grid->ntot_stag;

    CT_EXPAND(nv, Bs0[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                   Rs[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
    )
    #if TIME_STEPPING == RK4
    CT_EXPAND(nv, Bs1[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Bs2[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Bs3[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Bs4[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs4[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
    )
    #endif  // TIME_STEPPING == RK4

    #pragma acc enter data copyin(Bs0, Bs1, Bs2, Bs3, Bs4, Rs, Rs4)
    #pragma omp target enter data map(to: Bs0, Bs1, Bs2, Bs3, Bs4, Rs, Rs4)

    #endif
  }

  // 0b. - save initial stage
  //     - compute Vc = V(Uc) (necessary for fallback)

  #if FALLBACK_MODE == FALLBACK_BLENDING
  static Ary4D Ucn, Vcn;

  if (Ucn.isAllocated() == 0) {
    Ucn.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    Vcn.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
  }

  NFLUID_LOOP(nv) BOX_LOOP(boxt,i,j,k){
    data->Uc0(i,j,k,nv) = data->Uc(i,j,k,nv);
    Ucn(i,j,k,nv) = data->Uc(i,j,k,nv);
    Vcn(i,j,k,nv) = data->Vc(i,j,k,nv);
  }
  
double dt0 = Dts->dt;
Dts->dt = dt0*0.5;

  AdvanceStepEuler(data, grid, Dts);
  NFLUID_LOOP(nv) BOX_LOOP(boxt,i,j,k) U0(i,j,k,nv) = data->Uc(i,j,k,nv);

  AdvanceStepEuler(data, grid, Dts);
  NFLUID_LOOP(nv) BOX_LOOP(boxt,i,j,k) U0(i,j,k,nv) = data->Uc(i,j,k,nv);

  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Usafe(i,j,k,nv) = data->Uc(i,j,k,nv);  // Save initial stage
    data->Uc(i,j,k,nv)    = U0(i,j,k,nv);        // Restore initial stage
  }

  // Restore
  Dts->dt = dt0;
  NFLUID_LOOP(nv) BOX_LOOP(boxt,i,j,k){
    data->Uc(i,j,k,nv) = data->Uc0(i,j,k,nv) = Ucn(i,j,k,nv);
    data->Vc(i,j,k,nv) = Vcn(i,j,k,nv);
  }

  #endif

  // 0c. Apply boundary conditions on primitive vars. 

  g_intStage = 1;
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  RMHD_PCP_SetStageTime(data, Dts, 1);
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i, j, k) {
    data->flag(i, j, k) &= ~FLAG_PCP_CT_FAIL;
  }
  #endif
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  #if RING_AVERAGE > 1
  RingAverageCons(data, Rs, grid);
  #endif

  Boundary (data, 0, grid);   // Assign b.c. on data->Vc, data->Vs
  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);   // Set b.c. on point-values (primitive) 
  #endif

  // 0d. Check Field divergence (if required)

  #if CHECK_DIVB_CONDITION == 1
  CT_CheckDivB (data->Vs, grid, &box);    
  #elif CHECK_DIVB_CONDITION == 2
  CT_CheckDivB (data->Vs, grid, &boxt);    
  #endif

  // Inject particles
  #if PARTICLES != NO
  mynvtxstart_("Particles Inject",ORANGE);
  Particles_Inject(data, grid);
  mynvtxstop_();
  #endif

  // 0e. Flag zones lying in a shock or discontinuity 
  //     (Useful for shock flattening or entropy/energy selective update)

  mynvtxstart_("Shock",BLUE);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
  FlagShock (data, grid, &box);
  #endif
  #if (defined HIGH_ORDER) && ((HO_LAP_LIMITER != NO) || (HO_ORDER_REDUCTION != NO))
  DiscontinuityDetect(data, grid, &box);
  #endif
  mynvtxstop_();

  // 0f. Advance particles with pred. step 

  #if PARTICLES == PARTICLES_LP
  mynvtxstart_("Particles_LP_Update#1",SILVER);
  Particles_LP_Update (data, grid);
  mynvtxstop_();
  #endif
  
  // 0g. Store Vc copy for later CR particles update
  #if PARTICLES == PARTICLES_CR
  Ary4D V0;
  V0.create(nx1_tot, nx2_tot, nx3_tot, NVAR);
  BOX_LOOP(box,i,j,k) NVAR_LOOP(nv) V0(i,j,k,nv) = data->Vc(i,j,k,nv);
  #endif


  // ------------------------------------------------------
  // 1. First stage (EULER, RK2, RK3, RK4)
  // ------------------------------------------------------

  // 1a. First radiation step if IMEX_SSP2 is enabled

  #if RADIATION_RELATIVISTIC && RADIATION_IMEX_SSP2
  #pragma acc parallel loop collapse(4) present(data, box, U0)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U0)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    U0(i,j,k,nv) = data->Uc(i,j,k,nv);
  }
  ConsToPrim3D (data, 0, grid, &box);
  RadStep3D (data, Srad1, grid, dt_rk1);
  Boundary (data, 0, grid);
  #endif
  
  // 1b. Compute right hand side.

  #if PHYSICS == RMHD && PCP_FOFC == YES
  uint32_t stage_error_1 = RMHD_PCP_RunStage(
      data, Dts, dt, U0, Bs0, Rc, Rs, 0.0, 1.0, &box, grid);
  if (stage_error_1) {
    RMHD_PCP_SetStageTime(data, Dts, 0);
    #pragma acc exit data delete(box, boxt)
    #pragma omp target exit data map(delete: box, boxt)
    mynvtxstop_();
    return stage_error_1;
  }
  #else
  UpdateStage(data, Rc, Rs, grid, Dts, &box, dt);

  // 1b'. Update irradiation flux
  #if RADIATION_RELATIVISTIC && IRRADIATION_UPDATE
  RadIrradiationFlux(data, grid);
  #endif

  // 1c. Advance zone-centered variables

  #pragma acc parallel loop collapse(4) present(data, box, U0, U1, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U1, Rc) 
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    #if !(RADIATION_RELATIVISTIC && RADIATION_IMEX_SSP2)
    U0(i,j,k,nv) = data->Uc(i,j,k,nv);  // Save initial stage
    #endif
    data->Uc(i,j,k,nv) = U0(i,j,k,nv) + W1R*Rc(i,j,k,nv);
  }

  // 1d. Advance staggered fields

  #ifdef STAGGERED_MHD
  #pragma acc parallel loop collapse(3) present(data, box, Bs0, Rs)
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to:data[:1], box, Bs0, Rs)
  BOX_STAG_LOOP(box,i,j,k) {
    CT_EXPAND(nv, Bs0[nv](i,j,k) = data->Vs[nv](i,j,k);) // Save initial stage
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Bs0[nv](i,j,k) + W1R*Rs[nv](i,j,k);)
  }

  // 1d2. Replace cell-centered B with staggered average. 
  //      Apply energy correction if enabled.

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif  // DIVB_CONTROL == CONSTRAINED_TRANSPORT 

  // 1e. Convert to primitive variables

  #if RING_AVERAGE > 1
  RingAverageCons(data, Rs, grid);
  #endif

  // 1f. Radiation implicit step (second if IMEX_SSP2 is enabled)

  #if RADIATION_RELATIVISTIC
  #if RADIATION_IMEX_SSP2
  AddRadSource1(data, grid, Srad1, &box, dt_rk2);
  ConsToPrim3D (data, 0, grid, &box);
  RadStep3D (data, Srad2, grid, dt_rk1);
  #else
  ConsToPrim3D (data, 0, grid, &box);
  RadStep3D (data, grid, dt);
  #endif 
  #endif

  // 1g. Convert to primitive variables

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc 
  #else
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  uint32_t check_error_1 = RMHD_PCP_FinishStage(data, grid, &box);
  if (check_error_1) {
    RMHD_PCP_SetStageTime(data, Dts, 0);
    #pragma acc exit data delete(box, boxt)
    #pragma omp target exit data map(delete: box, boxt)
    mynvtxstop_();
    return check_error_1;
  }
  #else
  ConsToPrim3D (data, 0, grid, &box);
  #endif
  #endif

  #endif // PCP fluid stage / standard fluid stage

  // ------------------------------------------------------
  // 2. Second stage (RK2, RK3, RK4)
  // ------------------------------------------------------

#if TIME_STEPPING != EULER

  // 2a. Apply boundary conditions 

  g_intStage = 2;
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  RMHD_PCP_SetStageTime(data, Dts, 2);
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i, j, k) {
    data->flag(i, j, k) &= ~FLAG_PCP_CT_FAIL;
  }
  #endif
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);       // Set b.c. on point-values (primitive) 
  #else
  Boundary (data, 0, grid);       // Set b.c. on data->Vc 
  #endif

  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  // Preserve PCP's stage-wise flattening flags without replaying particles.
  #if SHOCK_FLATTENING == MULTID || ENTROPY_SWITCH
  FlagShock(data, grid, &box);
  #endif
  #endif

  // 2b. Advance particles 

  #if PARTICLES == PARTICLES_LP
  mynvtxstart_("Particles_LP_Update#2",SILVER);
  Particles_LP_Update (data, grid);
  mynvtxstop_();
  #endif

  // 2c. Compute right hand side.

  #if PHYSICS == RMHD && PCP_FOFC == YES
  uint32_t stage_error_2 = RMHD_PCP_RunStage(
      data, Dts, dt, U0, Bs0, Rc, Rs, W20, W21, &box, grid);
  if (stage_error_2) {
    RMHD_PCP_SetStageTime(data, Dts, 0);
    #pragma acc exit data delete(box, boxt)
    #pragma omp target exit data map(delete: box, boxt)
    mynvtxstop_();
    return stage_error_2;
  }
  #else
  UpdateStage(data, Rc, Rs, grid, Dts, &box, dt);

  #if RADIATION_RELATIVISTIC && IRRADIATION_UPDATE
  RadIrradiationFlux(data, grid);
  #endif

  // 2d. Radiation step if needed
  #if RADIATION_RELATIVISTIC && !RADIATION_IMEX_SSP2
  #pragma acc parallel loop collapse(4) present(data, box, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, Rc)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) += Rc(i,j,k,nv); // We need to compute this intermediate
                                        // stage of Uc before RadStep3D
  }
  ConsToPrim3D (data, 0, grid, &box);
  RadStep3D (data, grid, dt);
  #endif

  // 2e. Advance zone-centered variables

  #if !RADIATION || RADIATION_NR
  #pragma acc parallel loop collapse(4) present(data, box, U0, U1, U2, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U1, U2, Rc)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    #if TIME_STEPPING == RK4
    U1(i,j,k,nv) = data->Uc(i,j,k,nv);   // Save intermediate stage 
    #endif
    data->Uc(i,j,k,nv) = W20*U0(i,j,k,nv) + W21*data->Uc(i,j,k,nv) 
                                          + W2R*Rc(i,j,k,nv);
  }
  #else
  #if RADIATION_IMEX_SSP2
  #pragma acc parallel loop collapse(4) present(data, box, U0, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U0, Rc)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = 0.5*U0(i,j,k,nv) + 0.5*(data->Uc(i,j,k,nv) + Rc(i,j,k,nv));
  }
  AddRadSource2(data, grid, Srad1, Srad2, &box, dt_rk1, dt_rk3);
  #else
  #pragma acc parallel loop collapse(4) present(data, box, U0)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U0)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = W20*U0(i,j,k,nv) + W21*data->Uc(i,j,k,nv);
  }
  #endif
  #endif

  // 2f. Advance staggered fields

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #pragma acc parallel loop collapse(3) present(data, box, Bs0, Rs)
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box, Bs0, Rs) 
  BOX_STAG_LOOP(box,i,j,k) {
    #if TIME_STEPPING == RK4
    CT_EXPAND(nv, Bs1[nv](i,j,k) = data->Vs[nv](i,j,k)); // Save intermediate stage 
    #endif

    CT_EXPAND(nv, data->Vs[nv](i,j,k) =   W20*Bs0[nv](i,j,k) 
                                        + W21*data->Vs[nv](i,j,k) 
                                        + W2R*Rs[nv](i,j,k); )
  }

  // 2g. Replace cell-centered B with staggered average. 
  //     Apply energy correction if enabled.

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // DIVB_CONTROL == CONSTRAINED_TRANSPORT 

  #if RING_AVERAGE > 1
  RingAverageCons(data, Rs, grid);
  #endif

  // 2h. Apply FARGO orbital shift 

  #if (defined FARGO) && (TIME_STEPPING == RK2)
  FARGO_ShiftSolution (data, data->Uc, data->Vs, grid, run_config, Dts->dt);
  #endif

  // 2i. Convert to primitive variables

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost) 
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc 
  #else
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  uint32_t check_error_2 = RMHD_PCP_FinishStage(data, grid, &box);
  if (check_error_2) {
    RMHD_PCP_SetStageTime(data, Dts, 0);
    #pragma acc exit data delete(box, boxt)
    #pragma omp target exit data map(delete: box, boxt)
    mynvtxstop_();
    return check_error_2;
  }
  #else
  ConsToPrim3D (data, 0, grid, &box);
  #endif
  #endif


  #endif // PCP fluid stage / standard fluid stage

#endif //  TIME_STEPPING != EULER 

  // ------------------------------------------------------
  // 3. Third stage (RK3, RK4) 
  // ------------------------------------------------------

#if (TIME_STEPPING == RK3) || (TIME_STEPPING == RK4)
  
  // 3a. Apply Boundary conditions 
    
  g_intStage = 3;
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  RMHD_PCP_SetStageTime(data, Dts, 3);
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i, j, k) {
    data->flag(i, j, k) &= ~FLAG_PCP_CT_FAIL;
  }
  #endif
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);      // Set b.c. on point-values (primitive)
  #else
  Boundary (data, 0, grid);      // Set b.c. on data->Vc 
  #endif

  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  // Preserve PCP's stage-wise flattening flags without replaying particles.
  #if SHOCK_FLATTENING == MULTID || ENTROPY_SWITCH
  FlagShock(data, grid, &box);
  #endif
  #endif

  // 3b. Update particles 

  #if PARTICLES == PARTICLES_LP
  mynvtxstart_("Particles_LP_Update#3",SILVER);
  Particles_LP_Update (data, grid);
  mynvtxstop_();
  #endif

  // 3c. Compute right hand side.

  #if PHYSICS == RMHD && PCP_FOFC == YES
  uint32_t stage_error_3 = RMHD_PCP_RunStage(
      data, Dts, dt, U0, Bs0, Rc, Rs, W30, W31, &box, grid);
  if (stage_error_3) {
    RMHD_PCP_SetStageTime(data, Dts, 0);
    #pragma acc exit data delete(box, boxt)
    #pragma omp target exit data map(delete: box, boxt)
    mynvtxstop_();
    return stage_error_3;
  }
  #else
  UpdateStage(data, Rc, Rs, grid, Dts, &box, dt);

  #if RADIATION_RELATIVISTIC && IRRADIATION_UPDATE
  RadIrradiationFlux(data, grid);
  #endif


  // 3d. Radiation step if needed
  #if RADIATION_RELATIVISTIC && !RADIATION_IMEX_SSP2
  #pragma acc parallel loop collapse(4) present(data, box, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, Rc)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) += Rc(i,j,k,nv);
  }
  ConsToPrim3D (data, 0, grid, &box);
  RadStep3D (data, grid, dt);
  #endif

  // 3e. Advance zone-centered variables

  #if !RADIATION || RADIATION_NR
  #pragma acc parallel loop collapse(4) present(data, box, Rc, U0, U2, U3)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, Rc, U2, U3)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    #if TIME_STEPPING == RK4
    U2(i,j,k,nv) = data->Uc(i,j,k,nv);   // Save intermediate stage 
    #endif
    data->Uc(i,j,k,nv) =  W30*U0(i,j,k,nv) + W31*data->Uc(i,j,k,nv) 
                                           + W3R*Rc(i,j,k,nv);
  }
  #else
  // Already integrated
  #pragma acc parallel loop collapse(4) present(data, box, U0)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U0)
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = W30*U0(i,j,k,nv) + W31*data->Uc(i,j,k,nv);
  }
  #endif

  // 3g. Advance staggered variables

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #pragma acc parallel loop collapse(3) present(data, box, Bs0, Rs)
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box, Bs0, Rs) 
  BOX_STAG_LOOP(box,i,j,k) {
    #if TIME_STEPPING == RK4
    CT_EXPAND(nv, Bs2[nv](i,j,k) = data->Vs[nv](i,j,k)); // Save intermediate stage 
    #endif
    CT_EXPAND(nv, data->Vs[nv](i,j,k) =   W30*Bs0[nv](i,j,k) 
                                        + W31*data->Vs[nv](i,j,k) 
                                        + W3R*Rs[nv](i,j,k);)
  }

  // 3h. Replace cell-centered B with staggered average. 
  //     Apply energy correction if enabled.

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // DIVB_CONTROL == CONSTRAINED_TRANSPORT

  #if RING_AVERAGE > 1
  RingAverageCons(data, Rs, grid);
  #endif

  // 3h. Apply FARGO orbital shift 

  #if (defined FARGO) && (TIME_STEPPING == RK3)
  FARGO_ShiftSolution (data, data->Uc, data->Vs, grid, run_config, Dts->dt);
  #endif

  // 3i. Convert to primitive variables

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost) 
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc 
  #else
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  uint32_t check_error_3 = RMHD_PCP_FinishStage(data, grid, &box);
  if (check_error_3) {
    RMHD_PCP_SetStageTime(data, Dts, 0);
    #pragma acc exit data delete(box, boxt)
    #pragma omp target exit data map(delete: box, boxt)
    mynvtxstop_();
    return check_error_3;
  }
  #else
  ConsToPrim3D (data, 0, grid, &box);
  #endif
  #endif

  #endif // PCP fluid stage / standard fluid stage

#endif // TIME_STEPPING == RK3/RK4 

  // ------------------------------------------------------
  // 4. Fourth Stage (RK4) - High Order only
  // ------------------------------------------------------

#if TIME_STEPPING == RK4

  // 4a. Apply Boundary conditions 

  g_intStage = 4;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  // BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost) 
  // PointValue(data, grid, &box);       // Obtain point values Upc, Vpc 
  Boundary (data, 1, grid);         // Set b.c. on point-values (primitive) */

  // 4b. Compute right hand side.

  UpdateStage(data, R4, Rs4, grid, Dts, &box, dt);

  // 4c. Advance zone-centered variables

  #pragma acc parallel loop collapse(4) present(data, box, U0, U4, U3, R4)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U4, U3, R4) 
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    U3(i,j,k,nv) = data->Uc(i,j,k,nv);   // Save intermediate stage 
    data->Uc(i,j,k,nv) = W40*U0(i,j,k,nv) + W43*U3(i,j,k,nv) 
                                          + W4R*R4(i,j,k,nv);
  }

  // 4d. Advance staggered variables

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #pragma acc parallel loop collapse(3) present(data, box, Bs0, Bs3, Rs4)
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box, Bs0, Bs3, Rs4)
  BOX_STAG_LOOP(box,i,j,k) {
    CT_EXPAND(nv, Bs3[nv](i,j,k) = data->Vs[nv](i,j,k)); // Save intermediate stage 
    CT_EXPAND(nv, data->Vs[nv](i,j,k) =   W40*Bs0[nv](i,j,k) 
                                        + W43*Bs3[nv](i,j,k) 
                                        + W4R*Rs4[nv](i,j,k);)
  }

  // 4e. Replace cell-centered B with staggered average. 
  //     Apply energy correction if enabled.

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif  // DIVB_CONTROL == CONSTRAINED_TRANSPORT 

  #if RING_AVERAGE > 1
  RingAverageCons(data, Rs, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost) 
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc 
  #else
  ConsToPrim3D (data, 0, grid, &box);
  #endif

  // ------------------------------------------------------
  // 5. Fifth stage (RK4) - High Order only
  // ------------------------------------------------------

  // 5a. Apply Boundary conditions 

  g_intStage = 5;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  Boundary (data, 1, grid);           // Set b.c. on point-values (primitive)

  // 5b. Compute right hand side.

  UpdateStage(data, Rc, Rs, grid, Dts, &box, dt);

  // 5c. Advance zone-centered variables

  #pragma acc parallel loop collapse(4) present(data, box, U2, U3, U4, R4, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U2, U3, U4, R4, Rc) 
  NFLUID_LOOP(nv) BOX_LOOP(box,i,j,k){
    U4(i,j,k,nv) = data->Uc(i,j,k,nv);   // Save intermediate stage 
    data->Uc(i,j,k,nv) = W52*U2(i,j,k,nv)  + W53*U3(i,j,k,nv) + 
                         W5R4*R4(i,j,k,nv) + W54*U4(i,j,k,nv) +
                         W5R*Rc(i,j,k,nv);
  }

  // 5d. Advance staggered fields

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #pragma acc parallel loop collapse(3) present(data, box, Bs2, Bs3, Bs4, Rs4, Rs)
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box, Bs2, Bs3, Bs4, Rs4, Rs)
  BOX_STAG_LOOP(box,i,j,k) {

    CT_EXPAND(nv, Bs4[nv](i,j,k) = data->Vs[nv](i,j,k)); // Save intermediate stage 
    CT_EXPAND(nv, data->Vs[nv](i,j,k) =   W52 *Bs2[nv](i,j,k) 
                                        + W53 *Bs3[nv](i,j,k) 
                                        + W5R4*Rs4[nv](i,j,k) 
                                        + W54 *Bs4[nv](i,j,k) 
                                        + W5R *Rs[nv](i,j,k);)
  }

  // 5e. Replace cell-centered B with staggered average. 
  //     Apply energy correction if enabled.

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // DIVB_CONTROL == CONSTRAINED_TRANSPORT 
  
  #if RING_AVERAGE > 1
  RingAverageCons(data, Rs, grid);
  #endif

  // 5f. Apply FARGO orbital shift

  //   #ifdef FARGO
  //   FARGO_ShiftSolution (data->Uc, data->Vs, grid, run_config);
  //   #endif

  // 5g. Convert to primitive variables

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost) 
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc 
  #else
  ConsToPrim3D (data, 0, grid, &box);
  #endif

#endif // TIME_STEPPING == RK4 

  
  // ------------------------------------------------------
  // 6. Evolve particles.
  // ------------------------------------------------------

#if PARTICLES != NO
  // 6a. Evolve LP particle spectra.
  #if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == YES)
  mynvtxstart_("Particles_LP_Spectra",ORANGE);
  Boundary (data, 0, grid);
  Particles_LP_UpdateSpectra (data, Dts, grid);
  mynvtxstop_();
  #endif

  //6b. Evolve CR particles
  #if PARTICLES == PARTICLES_CR
  double temp;

  BOX_LOOP(box,i,j,k) NVAR_LOOP(nv) {
    temp = data->Vc(i,j,k,nv);
    data->Vc(i,j,k,nv) = data->Vc(i,j,k,nv)*0.5 + V0(i,j,k,nv)*0.5;
    V0(i,j,k,nv) = temp;
  }
      
  Particles_CR_Update(data, grid);
 
  BOX_LOOP(box,i,j,k) NVAR_LOOP(nv) {
    data->Vc(i,j,k,nv) = V0(i,j,k,nv);
  }
  #endif
#endif // PARTICLES != NO

  // ------------------------------------------------------
  // 7. Update irradiation flux (HD/MHD)
  // ------------------------------------------------------
  #if RADIATION_NR && IRRADIATION_UPDATE
  RadIrradiationFlux(data, grid);
  #endif

  // ------------------------------------------------------
  // 8. Conversion necessary for High Order, for 2 reasons:
  //    i)  reproduce the min/max values in the log file, 
  //        which are computed from data->Vc also for high 
  //        order methods.
  //    ii) allow next step to proceed correctly when 
  //        shock flattening / detection is enabled.
  // ------------------------------------------------------

  #ifdef HIGH_ORDER
  ConsToPrim3D (data, 0, grid, &box);
  #endif

  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  RMHD_PCP_SetStageTime(data, Dts, 0);
  #endif
  #pragma acc exit data delete (box, boxt)
  #pragma omp target exit data map(delete: box, boxt)

  mynvtxstop_();

  return 0;
}
