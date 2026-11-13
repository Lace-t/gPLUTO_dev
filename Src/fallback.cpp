////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Fallback integration function.

  When FALLBACK is activated (triggered) from other functions, 
  track zones marked with catastrophic errors and retry step by 
  lowering spatial accuracy.

      - Store solution at U(n):   
        U0 = data->Uc(), data->flag() = 0

      - Advance solution at the next time level with standard method,
        keeping track of previously activated flags (stored in flag_keep()).

        if (ttl(*) < 1) flag_keep(*) = 0;
        else            data->flag(*) = flag_keep(*);

        U(n) --> U(n+1)

      - On failure, collect current root flags separately from the
        cumulative fallback policy. Add HLL at each root, FLAT in its
        inner neighbourhood and MINMOD in the outer layer.

      - Restore U(n), V(n) and staggered B(n), then retry at the same dt.
        Newly failed roots extend the cumulative policy. Stop same-dt
        retries when no new root region appears or FALLBACK_LOCAL_PASSES
        is reached, then halve dt up to FALLBACK_MAX_RETRIES times.

      - Keep accumulated policy across retries and retain its actions in
        flag_keep() for the existing TTL lifetime across accepted steps.

  Here: 

      - flag_keep(): array recording persistent fallback policy across 
                    timesteps (NOT root-failure indicator). It basically stores, 
                    for each cell, persistent actions from step to step (as long
                    as its time to live > 0)
      - ttl():    Time-to-live counter controlling how long the fallback 
                  policy is kept (flag_keep() remains active for ttl() steps)

  \author  A. Mignone  (andrea.mignone@unito.it)
  \date    Jan 16, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifndef FALLBACK_BUFFER
  #define FALLBACK_BUFFER   1
#endif

#ifndef FALLBACK_TTL_MAX
  #define FALLBACK_TTL_MAX   2
#endif

#ifndef FALLBACK_SUBSTEPS 
  #define FALLBACK_SUBSTEPS 1
#endif

// --------------------------------------------------------
// Check keywords bounds 
// --------------------------------------------------------

#if FALLBACK_BUFFER < 0
  #error FALLBACK_BUFFER cannot be less than 0
#endif

#if FALLBACK_TTL_MAX < 1
  #error FALLBACK_TTL_MAX cannot be less than 1
#endif

#if FALLBACK_SUBSTEPS < 1
  #error FALLBACK_SUBSTEPS cannot be less than 1
#endif


#ifndef FALLBACK_LOCAL_PASSES
  #define FALLBACK_LOCAL_PASSES 4
#endif
#if FALLBACK_LOCAL_PASSES < 1
  #error FALLBACK_LOCAL_PASSES must be positive
#endif

#ifndef FALLBACK_MAX_RETRIES
  #define FALLBACK_MAX_RETRIES 1
#endif
#if FALLBACK_MAX_RETRIES < 0
  #error FALLBACK_MAX_RETRIES cannot be negative
#endif

#if FALLBACK_MODE == FALLBACK_RECOMPUTE
// Both the initial attempt and every recomputation use the same source/RK
// sequence. All ranks must agree before entering another attempt.
static uint32_t FallbackAdvance(Data *data, Grid *grid, timeStep *Dts)
{
  SplitSource(data, 0.5 * Dts->dt, Dts, grid);
  uint32_t err = AdvanceStep(data, grid, Dts);
#ifdef PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_UNSIGNED, MPI_BOR, MPI_COMM_WORLD);
#endif
  if (!err) {
    SplitSource(data, 0.5 * Dts->dt, Dts, grid);
  }
  return err;
}

// Restore the same complete fluid snapshot for every retry and on exhaustion.
static void FallbackRestore(Data *data, Grid *grid, Ary4D &Uc0, Ary4D &Vc0,
                            Ary3D *Vs0)
{
  int nx = grid->ntot[IDIR];
  int ny = grid->ntot[JDIR];
  int nz = grid->ntot[KDIR];
  #pragma acc parallel loop collapse(3) present(data, Uc0, Vc0)
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        for (int n = 0; n < NVAR; ++n) {
          data->Uc(i,j,k,n) = Uc0(i,j,k,n);
          #ifdef HIGH_ORDER
          data->Vpc(i,j,k,n) = Vc0(i,j,k,n);
          #else
          data->Vc(i,j,k,n) = Vc0(i,j,k,n);
          #endif
        }
      }
    }
  }
  #ifdef STAGGERED_MHD
  for (int n = 0; n < DIMENSIONS; ++n) {
    int ex = nx + (n == IDIR);
    int ey = ny + (n == JDIR);
    int ez = nz + (n == KDIR);
    #pragma acc parallel loop collapse(3) present(data, Vs0[:6])
    for (int k = 0; k < ez; ++k) {
      for (int j = 0; j < ey; ++j) {
        for (int i = 0; i < ex; ++i) {
          #ifdef HIGH_ORDER
          data->Vps[n](i,j,k) = Vs0[n](i,j,k);
          #else
          data->Vs[n](i,j,k) = Vs0[n](i,j,k);
          #endif
        }
      }
    }
  }
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Fallback(Data *data, Grid *grid, timeStep *Dts)
{
  int  i, j, k, nv;

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  int zone_count;
  uint32_t err;
  int retry = 0;
  int local_pass = 0;
  double inv_hyp0 = Dts->invDt_hyp;
  double inv_par0 = Dts->invDt_par;
  double accepted_dt;

  double t0, dt0;
  static Ary3Duint32_t policy; // Cumulative policy for this complete step
  static Ary3Duint32_t flag1; 
  static Ary3Duint32_t flag_keep;  // Bookkeep flagged zones from step to step   
  static Ary3Duint32_t ttl;        // "Time to live": flag lifetime (in steps)
  static Ary4D Uc0, Vc0;
  static Ary3D Vs0[6];
  RBox box, box_tot;

  int rank=0;

  // ------------------------------------------------------
  // 0. Define boxes and allocate memory, reset flag arrays 
  // ------------------------------------------------------

  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &box_tot);
  #pragma acc enter data copyin(box, box_tot)

  if (Uc0.isAllocated() == 0){
    Uc0.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    Vc0.create (nx1_tot,nx2_tot,nx3_tot,NVAR);

    policy.create(nx1_tot,nx2_tot,nx3_tot);
    flag1.create (nx1_tot,nx2_tot,nx3_tot);
    flag_keep.create (nx1_tot,nx2_tot,nx3_tot);
    ttl.create (nx1_tot,nx2_tot,nx3_tot);

    #pragma acc enter data copyin(Uc0, Vc0, policy, flag1, flag_keep, ttl)

    #pragma acc parallel loop collapse(3) present(flag1, flag_keep, ttl, box_tot)
    BOX_LOOP(box_tot, i,j,k) {
      flag1(i,j,k) = flag_keep(i,j,k) = 0;
      ttl(i,j,k)   = FALLBACK_TTL_MAX;
    }

    #ifdef STAGGERED_MHD 
    DIM_EXPAND (Vs0[BX1s].create((nx1_tot+1),nx2_tot,nx3_tot); ,
                Vs0[BX2s].create(nx1_tot,(nx2_tot+1),nx3_tot); ,
                Vs0[BX3s].create(nx1_tot,nx2_tot,(nx3_tot+1));)
    #pragma acc enter data copyin(Vs0)

    #endif
  }

  // ------------------------------------------------------
  // 1a. Save a copy of solution array at t^n (Uc0, Vs0)
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(data, box_tot, Uc0, Vc0)  
  BOX_LOOP(box_tot, i,j,k){
    NVAR_LOOP(nv) {
      Uc0(i,j,k,nv) = data->Uc(i,j,k,nv);
      #ifdef HIGH_ORDER
      Vc0(i,j,k,nv) = data->Vpc(i,j,k,nv);
      #else
      Vc0(i,j,k,nv) = data->Vc(i,j,k,nv);
      #endif
    }

  }

  #ifdef STAGGERED_MHD
  for (int n = 0; n < DIMENSIONS; ++n) {
    int ex = nx1_tot + (n == IDIR);
    int ey = nx2_tot + (n == JDIR);
    int ez = nx3_tot + (n == KDIR);
    #pragma acc parallel loop collapse(3) present(data, Vs0)
    for (int k = 0; k < ez; ++k) {
      for (int j = 0; j < ey; ++j) {
        for (int i = 0; i < ex; ++i) {
          #ifdef HIGH_ORDER
          Vs0[n](i,j,k) = data->Vps[n](i,j,k);
          #else
          Vs0[n](i,j,k) = data->Vs[n](i,j,k);
          #endif
        }
      }
    }
  }
  #endif

  t0  = Dts->time;
  dt0 = Dts->dt;
  accepted_dt = dt0;

  data->err_count        = 0;
  data->integration_mode = DEFAULT;
  #pragma acc update device(data->integration_mode)
  #pragma omp target update to(data->integration_mode)

  // ------------------------------------------------------
  // 1b. Enable/disable previously activated flags
  //     Count number of zones that have reached life limit
  // ------------------------------------------------------

  {
    int dcount = 0;
    int kcount = 0;
    #pragma acc parallel loop collapse(3) present(data,flag_keep,ttl,box_tot) \
                         reduction(+:dcount, kcount)
    BOX_LOOP(box_tot,i,j,k) {
      if (flag_keep(i,j,k) != 0) ttl(i,j,k)--;

      if (ttl(i,j,k) < 1){
        flag_keep(i,j,k) = 0;
        ttl(i,j,k)       = FALLBACK_TTL_MAX;
        dcount++; 
      }else{
        data->flag(i,j,k) = flag_keep(i,j,k);
        kcount += (flag_keep(i,j,k) != 0);
      }
    }

    #ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &dcount, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &kcount, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif

    // Print some diagnostic (note: here serial and parallel print
    // may give different total since the loop includes the ghost zones)
    // if (dcount > 0 || kcount > 0) {
    //   printLog ("> Fallback(): %d active cells [%d disabled]\n", 
    //              kcount, dcount);
    // }
  }

  // ------------------------------------------------------
  // 1c. Advance solution array w/ standard method
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(data, policy, box_tot)
  BOX_LOOP(box_tot, i,j,k) {
    policy(i,j,k) = data->flag(i,j,k);
  }

  err = FallbackAdvance(data, grid, Dts); 

  // ------------------------------------------------------
  // 2. Collect current failures and extend the cumulative local policy.
  //    All ranks agree on same-dt retries, halvings and acceptance.
  // ------------------------------------------------------

  while (true) {
    int progress = 0;
    zone_count = 0;
    // Keep only current failure seeds. In particular, do not reuse stale
    // ghost seeds from an earlier candidate or confuse policy with failure.
    #pragma acc parallel loop collapse(3) present(data, flag1, box_tot)
    BOX_LOOP(box_tot, i,j,k) {
      flag1(i,j,k) = 0;
    }
    #pragma acc parallel loop collapse(3) present(data, flag1, policy, box) \
                         reduction(+:zone_count) reduction(max:progress)
    BOX_LOOP(box, i,j,k) {
      if (data->flag(i,j,k) & FLAG_FALLBACK) {
        flag1(i,j,k) = FLAG_FALLBACK;
        ++zone_count;
        if (!(policy(i,j,k) & FLAG_HLL)) {
          progress = 1;
        }
      }
    }
    #ifdef PARALLEL
    int counts[2] = {zone_count, progress};
    MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    zone_count = counts[0];
    progress = (counts[1] != 0);
    #endif

    if (!err && zone_count == 0) {
      break;
    }

    // Match PLUTO: try local recovery at least once at each dt. Subsequent
    // passes at that dt are useful only if newly failed roots add a region.
    if (local_pass > 0 && (!progress || local_pass >= FALLBACK_LOCAL_PASSES)) {
      if (retry >= FALLBACK_MAX_RETRIES) {
        FallbackRestore(data, grid, Uc0, Vc0, Vs0);
        Dts->time = t0;
        data->run_config->time = t0;
        Dts->invDt_hyp = inv_hyp0;
        Dts->invDt_par = inv_par0;
        Dts->dt = dt0;
        data->integration_mode = DEFAULT;
        #pragma acc update device(data->run_config->time, data->integration_mode)
        #pragma acc exit data delete(box, box_tot)
        QuitPLUTO("Fallback: retries exhausted; state rolled back", 1);
      }
      ++retry;
      local_pass = 0;
      accepted_dt = ldexp(dt0, -retry);
    }

    // Exchange the roots BEFORE dilation so rank-boundary failures produce
    // the same neighbourhood on both sides. OR then merges overlap/corners.
    FlagBoundary(flag1, FLAG_COPY, grid);
    #pragma acc parallel loop collapse(3) present(policy, flag1, box_tot)
    BOX_LOOP(box_tot, i,j,k) {
      const int flat_radius = FALLBACK_BUFFER;
      const int outer_radius = flat_radius + 1;
      uint32_t added = 0;
      int ib = MAX(0, i - outer_radius*INCLUDE_IDIR);
      int ie = MIN(nx1_tot - 1, i + outer_radius*INCLUDE_IDIR);
      int jb = MAX(0, j - outer_radius*INCLUDE_JDIR);
      int je = MIN(nx2_tot - 1, j + outer_radius*INCLUDE_JDIR);
      int kb = MAX(0, k - outer_radius*INCLUDE_KDIR);
      int ke = MIN(nx3_tot - 1, k + outer_radius*INCLUDE_KDIR);
      // Gather from an immutable seed map: each cell has one writer.
      for (int kk = kb; kk <= ke; ++kk) {
        for (int jj = jb; jj <= je; ++jj) {
          for (int ii = ib; ii <= ie; ++ii) {
            if (!(flag1(ii,jj,kk) & FLAG_FALLBACK)) {
              continue;
            }
            if (ii == i && jj == j && kk == k) {
              added |= FLAG_HLL;
            }
            if (abs(ii-i) <= flat_radius && abs(jj-j) <= flat_radius &&
                abs(kk-k) <= flat_radius) {
              added |= FLAG_FLAT;
            } else {
              added |= FLAG_MINMOD;
            }
          }
        }
      }
      #ifdef HIGH_ORDER
      if (added) {
        added |= FLAG_HO_LAP_LIMITER;
      }
      #endif
      policy(i,j,k) |= added;
    }
    FlagBoundary(policy, FLAG_OR, grid);

    int zone_count_buf = 0;
    #pragma acc parallel loop collapse(3) present(data, policy, flag_keep, ttl, box_tot) \
                         reduction(+:zone_count_buf)
    BOX_LOOP(box_tot, i,j,k) {
      uint32_t actions = policy(i,j,k) & (FLAG_HLL | FLAG_FLAT | FLAG_MINMOD);
      #ifdef HIGH_ORDER
      actions |= policy(i,j,k) & FLAG_HO_LAP_LIMITER;
      #endif
      data->flag(i,j,k) = policy(i,j,k);
      if (actions) {
        flag_keep(i,j,k) |= actions;
        ttl(i,j,k) = FALLBACK_TTL_MAX;
        ++zone_count_buf;
      }
    }
    #ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &zone_count_buf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif

    FallbackRestore(data, grid, Uc0, Vc0, Vs0);
    Dts->time = t0;
    Dts->invDt_hyp = inv_hyp0;
    Dts->invDt_par = inv_par0;
    Dts->dt = accepted_dt / FALLBACK_SUBSTEPS;
    data->run_config->time = t0;
    data->integration_mode = FALLBACK_RECOMPUTE;
    data->err_count = 0;
    #pragma acc update device(data->run_config->time, data->integration_mode)
    ++local_pass;
    // printLog("> Fallback: halvings=%d local_pass=%d roots=%d buffered=%d new_region=%s\n",
    //          retry, local_pass, zone_count, zone_count_buf, progress ? "YES" : "NO");
    for (int substep = 0; substep < FALLBACK_SUBSTEPS; ++substep) {
      // printLog("> [Fallback] retry=%d local_pass=%d substep=%d/%d t=%.16e dt=%.16e\n",
      //          retry, local_pass, substep + 1, FALLBACK_SUBSTEPS, Dts->time, Dts->dt);
      err = FallbackAdvance(data, grid, Dts);
      if (err) {
        break;
      }
      Dts->time += Dts->dt;
      data->run_config->time = Dts->time;
      #pragma acc update device(data->run_config->time)
    }
  }

  // main advances the clock by the interval actually accepted.
  Dts->time = t0;
  Dts->dt = accepted_dt;

#ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == 0 &&  (retry > 0 || local_pass > 0)) {
    printLog("> Fallback success: t=%.4e dt=%.4e halvings=%d local_pass=%d\n",
             t0, accepted_dt, retry, local_pass);
  }

  // ------------------------------------------------------
  // 4. Quit check
  // ------------------------------------------------------

  data->integration_mode = DEFAULT;
  data->run_config->time = t0;
  #pragma acc update device(data->integration_mode, data->run_config->time)

  #pragma acc exit data delete(box, box_tot)

}
#endif //  FALLBACK_MODE == FALLBACK_RECOMPUTE

////////////////////////////////////////////////////////////////////////
void Fallback_Trigger(Ary3Duint32_t& err, Ary3Duint32_t& flag, RBox *box)
/*!
  Set which errors are going to trigger the fallback condition.
*/
////////////////////////////////////////////////////////////////////////
{
  int i,j,k;

  BOX_LOOP(*box,i,j,k){
    if (   (err(i,j,k) & ERR_CONS2PRIM_FAIL) 
        || (err(i,j,k) & ERR_ENG_SOLVE_NO_SOLUTION) 
        || (err(i,j,k) & ERR_NAN)        
        || (err(i,j,k) & ERR_NEGATIVE_RHO) 
        || (err(i,j,k) & ERR_NEGATIVE_ENG) 
        || (err(i,j,k) & ERR_NEGATIVE_PRS) 
        || (err(i,j,k) & ERR_POINT_VALUE_FAIL)
        || (err(i,j,k) & ERR_GAUSSIAN_SOLVE)
//        || (err(i,j,k) & ERR_RIEMANN_SOLVE_FAIL) 
        || (err(i,j,k) & ERR_QUIT)) {
      flag(i,j,k) |= FLAG_FALLBACK;
    }
  }
}
