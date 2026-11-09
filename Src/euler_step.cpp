////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Advance equations with 1st-order Euler method.

  Used in fallback approaches

FALLBACK_MODE = 0   -> disabled
FALLBACK_MODE = 1   -> crude average
FALLBACK_MODE = 2   -> + blending with U_safe
FALLBACK_MODE = 3   -> + conservative redistribution
FALLBACK_MODE = 4    + flux recomputation (MOOD)


integration_mode            Error ? 
   DEFAULT            FALLBACK_MODE = FALLBACK_MOOD
   FALLBACK_MOOD      FALLBACK_MODE = FALLBACK_BLENDING

When used in FALLBACK_MODE:

- Make a copy of current level, U0() = Un() 
- Euler step with 1st order scheme --> Usafe()
- Standard RK step, compute stage Uc
- If (! ConsToPrim (Uc, Vc)) Uc() = th*U0() + (1.0 - th)*Usafe()
- Convert Uc()

Uho(n+1) = U(n) - dt/dx*[Fho(+) - Fho(-)] = U(n) + Rho(n)
Ulo(n+1) = U(n) - dt/dx*[Flo(+) - Flo(-)] = U(n) + Rlo(n)

i:
U(n+1) = U(n) - dt/dx*[Flo(+) - Flo(-)]

i-1:
U(n+1) = U(n) - dt/dx*[Flo(-) - Fho(--)]

i+1:
U(n+1) = U(n) - dt/dx*[Fho(++) - Flo(+)]


  \authors A. Mignone (andrea.mignone@unito.it)\n

  \date    Jan 10, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Advance the equations by a single time step using unsplit 
  integrators based on the method of lines.
  
  \param [in,out]   data  pointer to Data structure
  \param [in]       grid  pointer to array of Grid structures
  \param [in,out]    Dts  pointer to time step structure
*/
////////////////////////////////////////////////////////////////////////
uint32_t AdvanceStepEuler(Data *data, Grid *grid, timeStep *Dts)
{
  mynvtxstart_("RKSTEP", GREEN);
  int nv, i, j, k;
  double dt = Dts->dt;

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int nx1_tot=grid->ntot[IDIR];
  int nx2_tot=grid->ntot[JDIR];
  int nx3_tot=grid->ntot[KDIR];

  RBox box, boxt;
  const  Ary4D &U0 = data->Uc0;
  static Ary4D Rc;
  static Ary3D Bs0[6], Rs[6];

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine(0, nx1_tot - 1, 0, nx2_tot - 1, 0, nx3_tot - 1, &boxt);
  #pragma acc enter data copyin(box, boxt)
  #pragma omp target enter data map(to: box, boxt)

  // 0a. Allocate memory 

  if (data->Usafe.isAllocated() == 0){
    data->Usafe.create(nx1_tot, nx2_tot, nx3_tot,NVAR);  
  }

  if (Rc.isAllocated() == 0) {
    Rc.create (nx1_tot,nx2_tot,nx3_tot,NVAR);

    #pragma acc enter data copyin(U0)
    #pragma omp target enter data map(to: U0)

    #ifdef STAGGERED_MHD
    int (*ns)[3];
    ns = grid->ntot_stag;

    CT_EXPAND(nv, Bs0[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                   Rs[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
    )

    #pragma acc enter data copyin(Bs0)
    #pragma omp target enter data map(to: Bs0)

    #endif
  }

  // 0b. Set fallback options 

  RiemannSolver *solver  = data->fluidRiemannSolver;
  #if RADIATION_RELATIVISTIC
  RiemannSolver *rad_solver = data->radiationRiemannSolver;
  #endif
  data->integration_mode = FALLBACK_BLENDING;
  data->use_flat         = YES;
  data->fluidRiemannSolver = SetSolver ("hll");
  #if RADIATION_RELATIVISTIC
  data->radiationRiemannSolver = RadSetSolver ("tvdlf");
  #endif

  // 0c. Apply boundary conditions on primitive vars. 

  g_intStage = 1;
  //#pragma acc update device(g_intStage)
  //#pragma omp target update to(g_intStage)

  Boundary (data, 0, grid);   // Assign b.c. on data->Vc, data->Vs 

  // 0d. Flag zones lying in a shock or discontinuity 
  //     (Useful for shock flattening or entropy/energy selective update)

  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
  FlagShock (data, grid, &box);
  #endif

  // 0e. Advance particles with pred. step 

  // ------------------------------------------------------
  // 1. First stage
  // ------------------------------------------------------

  // 1b. Compute right hand side.

  UpdateStage(data, Rc, Rs, grid, Dts, &box, dt);

  // 1c. Advance zone-centered variables

  #pragma acc parallel loop collapse(4) present(data, box, U0, Rc)
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box, U0, Rc)
  NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
    data->Uc(i,j,k,nv) = U0(i,j,k,nv) + Rc(i,j,k,nv);
  }

  // 1d. Advance staggered fields

  #ifdef STAGGERED_MHD
  #pragma acc parallel loop collapse(3) present(data, box, Bs0, Rs)
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to:data[:1], box, Bs0, Rs)
  BOX_STAG_LOOP(box,i,j,k) {
    CT_EXPAND(nv, Bs0[nv](i,j,k) = data->Vs[nv](i,j,k);) // Save initial stage
  }

  // 1e. Replace cell-centered B with staggered
  //     average. Apply energy correction if enabled.

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif  // DIVB_CONTROL == CONSTRAINED_TRANSPORT 

  #if RADIATION_RELATIVISTIC
  ConsToPrim3D (data, 0, grid, &box);
  RadStep3D (data, grid, dt);
  #endif

  ConsToPrim3D (data, 0, grid, &box);

  #pragma acc exit data delete (box, boxt)
  #pragma omp target exit data map(delete: box, boxt)

  mynvtxstop_();

  // 0?. Restore default integration options

  data->integration_mode   = DEFAULT;
  data->use_flat           = NO;
  data->fluidRiemannSolver = solver;
  #if RADIATION_RELATIVISTIC
  data->radiationRiemannSolver = rad_solver;
  #endif

  return 0;
}
