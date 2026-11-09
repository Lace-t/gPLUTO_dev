#include "gpluto.hpp"

#pragma acc routine seq
int RMHD_PCP_Recover(double *u, double *v, RunConfig &run_config, int allow_fix)
{
  double work[NVAR], primitive[NVAR];
  for (int n = 0; n < NVAR; ++n) {
    work[n] = u[n];
    primitive[n] = 0.0;
    if (!isfinite(work[n])) {
      return 1;
    }
  }
  // EnergySolve recovers fluid variables; magnetic fields map directly.
  primitive[BX1] = work[BX1];
  primitive[BX2] = work[BX2];
  primitive[BX3] = work[BX3];

  int admissible = RMHD_PCP_Check(work);
  uint32_t err = 0;
  if (!admissible) {
    if (!allow_fix || work[RHO] <= 0.0) {
      return 1;
    }
    primitive[BX1] = work[BX1];
    primitive[BX2] = work[BX2];
    primitive[BX3] = work[BX3];
    err = RMHD_PressureFix(work, primitive, run_config);
  } else {
    State state;
    STATE_SET(state, primitive, work, NULL, 0.0, 0.0, 0.0);
    uint32_t flag = 0;
    state.flag = &flag;
    if (allow_fix) {
      err = ConsToPrim(state, run_config);
    } else {
      // FOFC and limiter trials must not commit a pressure-floor repair.
      err = RMHD_EnergySolve(work, primitive, run_config);
#if PCP_NR == YES
      if (err && !RMHD_PCP_NR(work, primitive, run_config)) {
        err = 0;
      }
#endif
      if (err) {
        return 1;
      }
    }
  }
  if (err & (ERR_CONS2PRIM_FAIL | ERR_PRESSURE_FIX)) {
    return 1;
  }
  if (!allow_fix) {
    // The strict path cannot alter U: identical input retains its admission.
    // Reuse the entrance criterion instead of evaluating it a second time.
    for (int n = 0; n < NVAR; ++n) {
      if (work[n] != u[n]) {
        return 1;
      }
    }
  } else if (!RMHD_PCP_Check(work)) {
    // PressureFix can change energy, so repaired states need a fresh check.
    return 1;
  }
#if NSCL > 0
  for (int n = NFLX; n < NVAR; ++n) {
    primitive[n] = work[n] / work[RHO];
  }
#endif
#ifdef GLM_MHD
  primitive[PSI_GLM] = work[PSI_GLM];
#endif
  double v2 = primitive[VX1] * primitive[VX1] + primitive[VX2] * primitive[VX2] +
              primitive[VX3] * primitive[VX3];
  if (!(v2 < 1.0) || primitive[RHO] <= 0.0 || primitive[PRS] <= 0.0) {
    return 1;
  }
  for (int n = 0; n < NVAR; ++n) {
    if (!isfinite(primitive[n])) {
      return 1;
    }
  }
  for (int n = 0; n < NVAR; ++n) {
    u[n] = work[n];
    v[n] = primitive[n];
  }
  return 0;
}

// The face at sweep index i lies between cells i-1 and i.
void RMHD_PCP_LimitFaces(Data *data, Grid *grid, RBox *box, int dir)
{
  int ib, ie, jb, je, kb, ke, ni, nj;
  ARRAY_STRIDES(grid, dir, ni, nj);
  if (dir == IDIR) {
    RBoxPermuteIndices<IDIR>(box, ib, ie, jb, je, kb, ke);
  } else if (dir == JDIR) {
    RBoxPermuteIndices<JDIR>(box, ib, ie, jb, je, kb, ke);
  } else {
    RBoxPermuteIndices<KDIR>(box, ib, ie, jb, je, kb, ke);
  }
  RunConfig &run_config = *data->run_config;
#pragma acc parallel loop collapse(3) present(data, grid, run_config)
  for (int k = kb; k <= ke; ++k) {
    for (int j = jb; j <= je; ++j) {
      for (int i = ib; i <= ie + 1; ++i) {
        long off = ni * (j + nj * k);
        Ary2D vc(&data->sweep->vC[NVAR * off], ni, NVAR);
        Ary2D vL(&data->sweep->vL[NVAR * off], ni, NVAR);
        Ary2D vR(&data->sweep->vR[NVAR * off], ni, NVAR);
        double trial_v[2][NVAR];
        bool changed[2] = {false, false};
        bool failed = false;
        for (int side = 0; side < 2; ++side) {
          double face_u[NVAR], center_u[NVAR], center_v[NVAR];
          for (int n = 0; n < NVAR; ++n) {
            trial_v[side][n] = side == 0 ? vL(i, n) : vR(i, n);
          }
          PrimToCons(trial_v[side], face_u, run_config);
          if (RMHD_PCP_Check(face_u)) {
            continue;
          }
          for (int n = 0; n < NVAR; ++n) {
            center_v[n] = vc(i - 1 + side, n);
          }
          PrimToCons(center_v, center_u, run_config);
          int limited = RMHD_PCP_Limit(center_u, face_u);
          if (limited < 0 || RMHD_PCP_Recover(face_u, trial_v[side], run_config, 0)) {
            failed = true;
            break;
          }
          changed[side] = (limited > 0);
        }
        // Keep interface updates atomic: discard both trials if either fails.
        if (!failed) {
          if (changed[0]) {
            for (int n = 0; n < NVAR; ++n) {
              vL(i, n) = trial_v[0][n];
            }
          }
          if (changed[1]) {
            for (int n = 0; n < NVAR; ++n) {
              vR(i, n) = trial_v[1][n];
            }
          }
        }
      }
    }
  }
}

#if PCP_FOFC == YES
#include "pcp_flux.inc"
#endif

#if (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
#if PCP_FOFC == YES
// Device snapshots include the upper face in each component's own direction.
struct PCPSnapshot {
  Ary4D u, v;
  Ary3D b[6];
  void create(Grid *g)
  {
    int nx = g->ntot[0], ny = g->ntot[1], nz = g->ntot[2];
    u.create(nx, ny, nz, NVAR);
    v.create(nx, ny, nz, NVAR);
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      b[n].create(nx + (n == 0), ny + (n == 1), nz + (n == 2));
    }
#endif
#pragma acc enter data copyin(this[ : 1])
  }
  void transfer(Data *d, Grid *g, int save)
  {
    int nx = g->ntot[0], ny = g->ntot[1], nz = g->ntot[2];
#pragma acc parallel loop collapse(3) present(d, u, v)
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          for (int n = 0; n < NVAR; ++n) {
            if (save) {
              u(i, j, k, n) = d->Uc(i, j, k, n);
              v(i, j, k, n) = d->Vc(i, j, k, n);
            } else {
              d->Uc(i, j, k, n) = u(i, j, k, n);
              d->Vc(i, j, k, n) = v(i, j, k, n);
            }
          }
        }
      }
    }
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      int ex = nx + (n == 0), ey = ny + (n == 1), ez = nz + (n == 2);
#pragma acc parallel loop collapse(3) present(d, b)
      for (int k = 0; k < ez; ++k) {
        for (int j = 0; j < ey; ++j) {
          for (int i = 0; i < ex; ++i) {
            if (save) {
              b[n](i, j, k) = d->Vs[n](i, j, k);
            } else {
              d->Vs[n](i, j, k) = b[n](i, j, k);
            }
          }
        }
      }
    }
#endif
  }
};
#endif

// Check and recover a completed candidate without changing its stage input.
// PressureFix is enabled only after all local FOFC attempts have failed.
static int PCPRecoverStage(Data *data, Grid *grid, RBox &box, RunConfig &run_config,
                           int allow_pressure_fix)
{
  int i, j, k;
  int bad = 0;
#pragma acc parallel loop collapse(3) present(data, grid, run_config, box) reduction(+ : bad)
  BOX_LOOP (box, i, j, k) {
    double u[NVAR], v[NVAR];
    for (int n = 0; n < NVAR; ++n) {
      u[n] = data->Uc(i, j, k, n);
    }
#if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    ConsToOrtho(u, metric);
#endif
    int err = (data->flag(i, j, k) & FLAG_PCP_CT_FAIL)
                  ? 1
                  : RMHD_PCP_Recover(u, v, run_config, allow_pressure_fix);
    if (err) {
      data->flag(i, j, k) |= FLAG_PCP_FOFC | FLAG_FALLBACK;
      ++bad;
    } else {
      data->flag(i, j, k) &= ~FLAG_FALLBACK;
#if GEOMETRY != CARTESIAN
      OrthoToCons(u, metric);
#endif
      for (int n = 0; n < NVAR; ++n) {
        data->Uc(i, j, k, n) = u[n];
        data->Vc(i, j, k, n) = v[n];
      }
    }
  }
#ifdef PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  return bad;
}

void RMHD_PCP_SetStageTime(Data *data, timeStep *Dts, int stage_number)
{
  RunConfig &run_config = *data->run_config;
  double offset = stage_number == 2 ? Dts->dt
                  : stage_number == 3 ? 0.5 * Dts->dt
                                      : 0.0;
  run_config.time = Dts->time + offset;
#pragma acc update device(run_config.time)
}

uint32_t RMHD_PCP_FinishStage(Data *data, Grid *grid, RBox *box)
{
  int bad = PCPRecoverStage(data, grid, *box, *data->run_config, 0);
#if PCP_PRESSURE_FIX == YES
  if (bad) {
    bad = PCPRecoverStage(data, grid, *box, *data->run_config, 1);
    //printLog("> PCP PressureFix: stage=%d rejected=%d\n", g_intStage, bad);
  }
#endif
  if (bad) {
#if FALLBACK_MODE == FALLBACK_RECOMPUTE
    return ERR_CONS2PRIM_FAIL;
#else
    QuitPLUTO("ERR_CONS2PRIM_FAIL", 1);
#endif
  }
  return 0;
}

#if PCP_FOFC == YES
// The RK driver owns boundaries, particles, radiation and source scheduling.
// Ubase/Bbase are its initial RK state; rhs/rs are its reusable work arrays.
uint32_t RMHD_PCP_RunStage(Data *data, timeStep *Dts, double dt,
                           const Ary4D &Ubase, Ary3D (&Bbase)[6],
                           Ary4D &rhs, Ary3D (&rs)[6], double initial_weight,
                           double stage_weight, RBox *box, Grid *grid)
{
  static PCPSnapshot stage;
  static bool allocated = false;
  if (!allocated) {
    stage.create(grid);
    allocated = true;
  }
  stage.transfer(data, grid, 1);
  RunConfig &run_config = *data->run_config;
  int nx = grid->ntot[0], ny = grid->ntot[1], nz = grid->ntot[2];
  const int step_retry = (data->integration_mode == FALLBACK_RECOMPUTE);
  int i, j, k;
  if (g_intStage == 1) {
#pragma acc parallel loop collapse(3) present(data, Ubase, box)
    BOX_LOOP ((*box), i, j, k) {
      for (int n = 0; n < NVAR; ++n) {
        Ubase(i, j, k, n) = data->Uc(i, j, k, n);
      }
    }
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      int ex = box->iend + (n == 0);
      int ey = box->jend + (n == 1);
      int ez = box->kend + (n == 2);
#pragma acc parallel loop collapse(3) present(data, Bbase, box)
      for (int k = box->kbeg; k <= ez; ++k) {
        for (int j = box->jbeg; j <= ey; ++j) {
          for (int i = box->ibeg; i <= ex; ++i) {
            Bbase[n](i, j, k) = data->Vs[n](i, j, k);
          }
        }
      }
    }
#endif
  }
#pragma acc parallel loop collapse(3) present(data)
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        data->flag(i, j, k) &= ~(FLAG_PCP_FOFC | FLAG_PCP_CT_FAIL);
        // Whole-step fallback selects only troubled roots, not every face.
        if (step_retry > 0 && (data->flag(i, j, k) & FLAG_HLL)) {
          data->flag(i, j, k) |= FLAG_PCP_FOFC;
        }
      }
    }
  }
  const int local_passes = PCP_STAGE_LOCAL_PASSES;
  int bad = 0;
  for (int pass = 0; pass < local_passes; ++pass) {
    if (pass > 0) {
      stage.transfer(data, grid, 0);
    }
    // Keep COPY on later RK stages for other flag consumers; retain the
    // established first-stage optimization after the initial FlagShock.
#if SHOCK_FLATTENING == MULTID && ENTROPY_SWITCH == NO
    if (pass > 0 || g_intStage != 1) {
      FlagBoundary(data->flag, FLAG_COPY, grid);
    }
#else
    FlagBoundary(data->flag, FLAG_COPY, grid);
#endif
    const bool pcp_flux_possible = (step_retry > 0 || pass > 0);
    UpdateStage(data, rhs, rs, grid, Dts, box, dt, pcp_flux_possible);
#pragma acc parallel loop collapse(3) present(data, rhs, box, Ubase, stage.u)
    BOX_LOOP ((*box), i, j, k) {
      for (int n = 0; n < NVAR; ++n) {
        data->Uc(i, j, k, n) = initial_weight * Ubase(i, j, k, n) +
                               stage_weight * (stage.u(i, j, k, n) + rhs(i, j, k, n));
      }
    }
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      int ex = box->iend + (n == 0);
      int ey = box->jend + (n == 1);
      int ez = box->kend + (n == 2);
#pragma acc parallel loop collapse(3) present(data, box, rs, Bbase, stage.b)
      for (int k = box->kbeg; k <= ez; ++k) {
        for (int j = box->jbeg; j <= ey; ++j) {
          for (int i = box->ibeg; i <= ex; ++i) {
            data->Vs[n](i, j, k) = initial_weight * Bbase[n](i, j, k) +
                                   stage_weight * (stage.b[n](i, j, k) + rs[n](i, j, k));
          }
        }
      }
    }
#pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP ((*box), i, j, k) {
      data->flag(i, j, k) &= ~FLAG_PCP_CT_FAIL;
    }
    CT_AverageStaggeredFields(data, 1, box, grid);
#endif
    bad = PCPRecoverStage(data, grid, *box, run_config, 0);
    if (!bad) {
      return 0;
    }
    // printLog("> PCP: stage=%d pass=%d retry=%d rejected=%d dt=%.16e\n",
    //          g_intStage, pass, step_retry, bad, dt);
  }
#if PCP_PRESSURE_FIX == YES
  // Exhaust FOFC before pressure repair; a remaining failure goes to main.
  bad = PCPRecoverStage(data, grid, *box, run_config, 1);
  // printLog("> PCP PressureFix: stage=%d retry=%d rejected=%d\n",
  //          g_intStage, step_retry, bad);
#endif
  if (bad) {
#if FALLBACK_MODE == FALLBACK_RECOMPUTE
    return ERR_CONS2PRIM_FAIL;
#else
    QuitPLUTO("ERR_CONS2PRIM_FAIL", 1);
#endif
  }
  return 0;
}
#endif // PCP_FOFC
#endif // stage checks / recovery
