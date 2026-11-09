/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations with IMEX Runge Kutta time integrators.

  Main driver for SSP IMEX Runge-Kutta time evolution.
  Time stepping methods include 
  - IMEX-SSP2 (2,2,2)
  - IMEX-SSP3 (3,3,2)
  - IMEX-SSP3 (4,3,3)

  - "Implicit-Explicit Runge-Kutta Schemes and Applications to
     Hyperbolic Systems with Relaxation", Journal Sci. Comput (2005)
  - "Beyond ideal MHD: towards a more realistic modelling of relativistic
     astrophysical plasmas"\n
     Palenzuela \& al, RAS (2009) 394, 1727-1740

  Debug: 
// {  ////////////////////////////////////////////////////////////
// RBox tbox;
// Ary4D Qnull;
// Ary3D Vsnull[6];
// RBoxDefine (grid->lbeg, grid->lend, &tbox);
// ShowData (data->Uc, data->Vc, data->Vs, &tbox, "RKStep(stage=1)");
// }  ///////////////////////////////////////////////////////////////

  \authors A. Mignone (andrea.mignone@unito.it)\n
           G. Mattia

  \date    Dec 17, 2025
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "gpluto.hpp"

#define BVEC_LOOP(nv)   for (nv = BX1; nv <= BX3; nv++)
#define EVEC_LOOP(nv)   for (nv = EX1; nv <= EX3; nv++)

#if TIME_STEPPING == RK4
 #error RK4 never tested for this configuration
#endif

/* ********************************************************************* */
/*!
 * Advance the equations by a single time step using unsplit
 * integrators based on the method of lines.
 *
 * \param [in,out]   data  pointer to Data structure
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *
 *********************************************************************** */
uint32_t AdvanceStep (Data *data, Grid *grid, timeStep *Dts)
{
  int  i, j, k, nv;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int nx1_tot=grid->ntot[IDIR];
  int nx2_tot=grid->ntot[JDIR];
  int nx3_tot=grid->ntot[KDIR];

  double dt = Dts->dt;
  double ci;

  const double gIMEX = (1.0 - 1.0/sqrt(2.0));
  const double alpha = 0.24169426078821;
  const double beta  = 0.06042356519705;
  const double eta   = 0.12915286960590;
  const double chi   = 0.5-beta-eta-alpha;

  RBox   box, box_tot;
  RBox   boxp;

  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &box_tot);
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &boxp);
  #else
  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &boxp);
  #endif

  #pragma acc enter data copyin (box, box_tot, boxp)
  #pragma omp target enter data map(to: box, box_tot, boxp)

  const  Ary4D &U0 = data->Uc0;
  static Ary4D S1, S2, S3, S4, Rc, R1, R2, R3, R4;
  static Ary4D Jf1, Jf2, Jf3, Jf4;
  static Ary3D Vs0[6], Rs[6], Rs1[6], Rs2[6], Rs3[6], Rs4[6];

  if (Rc.isAllocated() == 0) {

    Rc.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R1.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R2.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R3.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R4.create (nx1_tot,nx2_tot,nx3_tot,NVAR);

    S1.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S2.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S3.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S4.create (nx1_tot,nx2_tot,nx3_tot,NVAR);

    #pragma omp target enter data map(to: Rc, R1, R2, R3, R4, S1, S2, S3, S4)
    #pragma acc enter data copyin(Rc, R1, R2, R3, R4, S1, S2, S3)

    #ifdef STAGGERED_MHD
    Jf1.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf2.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf3.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf4.create (nx1_tot,nx2_tot,nx3_tot,3);

    int (*ns)[3], dir;
    ns = grid->ntot_stag;
    
    CT_EXPAND(nv, dir = nv%DIMENSIONS;
                  Vs0[nv].create(ns[dir][IDIR], ns[dir][JDIR], ns[dir][KDIR]);
                   Rs[nv].create(ns[dir][IDIR], ns[dir][JDIR], ns[dir][KDIR]);
                  Rs1[nv].create(ns[dir][IDIR], ns[dir][JDIR], ns[dir][KDIR]);
                  Rs2[nv].create(ns[dir][IDIR], ns[dir][JDIR], ns[dir][KDIR]);
                  Rs3[nv].create(ns[dir][IDIR], ns[dir][JDIR], ns[dir][KDIR]);
                  Rs4[nv].create(ns[dir][IDIR], ns[dir][JDIR], ns[dir][KDIR]);
    )

    #pragma acc enter data copyin(Jf1, Jf2, Jf3, Jf4, Vs0, Rs, Rs1, Rs2, Rs3, Rs4)
    #endif // STAGGERED_MHD 
  }

  mynvtxstart_("RKSTEP",GREEN);

  g_intStage = 10;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  Boundary (data, 0, grid);

  // Copy conservative variables at level 0 
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box[:1] , U0[:1])
  #pragma acc parallel loop collapse(4) present(data, box, U0)
  BOX_LOOP(box,i,j,k) NVAR_LOOP(nv) U0(i,j,k,nv) = data->Uc(i,j,k,nv);

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to:data[:1], box_tot[:1])
  #pragma acc parallel loop collapse(3) present(data, box_tot)
  BOX_LOOP(box_tot, i,j,k) {
    CT_EXPAND(nv, Vs0[nv](i,j,k) = data->Vs[nv](i,j,k);)
  }
  #endif // STAGGERED_MHD 

  // ------------------------------------------------------
  // 1. IMEX Stage #1 - Implicit contribution:
  //
  //    U1 = U^n + gamma*dt*S1                 [RK2, RK3]
  //    U1 = U^n + alpha*dt*S1                 [RK4]
  // ------------------------------------------------------

  #if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)
  double dt1 = Dts->dt*gIMEX;
  #elif TIME_STEPPING == RK4  
  double dt1 = Dts->dt*alpha;
  #endif

//  IMEXImplicitUpdate(data, dt1, S1, Jf1, grid, &boxp);   // d->Uc = U^1 
  IMEXExplicitSource (data, S1, grid);

  // ------------------------------------------------------
  // 2a. IMEX stage #2 - Explicit contribution:
  //
  //     U^* = U^n + dt*R1 + dt*(1-2*gamma)*S1     [RK2, RK3]
  //     U^* = U^n - dt*alpha*S1                   [RK4]
  // ------------------------------------------------------

  g_intStage = 1;
  #pragma acc update device(g_intStage)
  #pragma omp target update to(g_intStage)

  Boundary (data, 0, grid);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH)
  FlagShock (data, grid, &box);
  #endif
  UpdateStage (data, R1, Rs1, grid, Dts, &box, dt);

mynvtxstart_("Stage#1",BLUE);
  
  #if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)
  ci = dt*(1.0 - 2.0*gIMEX); // IMEX SSP2 (2,2,2), SSP3 (3,3,2) 
  #elif TIME_STEPPING == RK4  
  ci = -dt*alpha;            // IMEX SSP3 (4,3,3) 
  #endif

  // 2a1. Update zone-centered variables
  #pragma omp target teams distribute parallel for collapse(3) \
                     map (to: data[:1], box[:1], U0, R1, S1[:1] , data->Uc[:1])
  #pragma acc parallel loop collapse(3) present(data, box, U0, R1, S1)
  BOX_LOOP(box, i,j,k){
    #if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + R1(i,j,k,nv);
    #endif

    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*S1(i,j,k,nv);

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*S1(i,j,k,ENTR);
    #endif
  }

  // 2a2. Update staggered variables

  #ifdef STAGGERED_MHD 
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box[:1])
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){

    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + Rs1[nv](i,j,k);)
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT 
    #pragma acc loop seq
    DIM_LOOP(nv) data->Vs[EX1s+nv](i,j,k) -= ci*Jf1(i,j,k,nv);
    #endif
  }
  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // STAGGERED_MHD 

mynvtxstop_();

  // 2a3. Convert to primitive

  ConsToPrim3D (data, 0, grid, &box);
  
  // ------------------------------------------------------
  // 2b. IMEX stage #2 - Implicit contribution:
  //
  //     U2 = U^* + dt*gamma*S2                 [RK2, RK3]
  //     U2 = U^* + dt*alpha*S2                 [RK4]
  // ------------------------------------------------------

  g_intStage = 20;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  #if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)
  double dt2  = Dts->dt*gIMEX;
  #elif TIME_STEPPING == RK4  
  double dt2  = Dts->dt*alpha;
  #endif

//  IMEXImplicitUpdate(data, dt2, S2, Jf2, grid, &boxp);  // d->Uc = U^2 
  IMEXExplicitSource (data, S2, grid);

  // ------------------------------------------------------
  // 3a. IMEX stage #3 - Explicit contribution
  //
  //     U^(n+1) = U^n + dt*(R1 + R2)/2 + dt*(S1 + S2)/2      [RK2]
  //     U^*     = U^n + dt*(R1 + R2)/4 + dt*(1-2*gamma)*S1   [RK3]
  //     U^*     = U^n + dt*R2 + dt*(1-alpha)*S1              [RK4]
  // ------------------------------------------------------

  g_intStage = 2;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  Boundary (data, 0, grid);
  UpdateStage (data, R2, Rs2, grid, Dts, &box, dt);

  #if TIME_STEPPING == RK2    // IMEX SSP2 (2,2,2) 
  ci = 0.5*dt;
  #elif TIME_STEPPING == RK3  // IMEX SSP3 (3,3,2) 
  ci = Dts->dt*(0.5 - gIMEX);
  #elif TIME_STEPPING == RK4  // IMEX SSP3 (4,3,3) 
  ci = dt*(1.0 - alpha);
  #endif

mynvtxstart_("Stage#1",BLUE);
  
  // 3a1. Update zone-centered variables 
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box[:1], U0[:1], R1[:1], R2[:1], S1[:1], S2[:1])
  #pragma acc parallel loop collapse(3) present(data, box, U0, R1, R2, S1, S2)
  BOX_LOOP(box, i,j,k){
    #if TIME_STEPPING == RK2
      #pragma acc loop
      NVAR_LOOP(nv) data->Uc(i,j,k,nv)  =   U0(i,j,k,nv) 
                                          + 0.5*(R1(i,j,k,nv) + R2(i,j,k,nv));
      #pragma acc loop
      EVEC_LOOP(nv) data->Uc(i,j,k,nv) += 0.5*dt*(S1(i,j,k,nv) + S2(i,j,k,nv));

      #if ENTROPY_SWITCH
      data->Uc(i,j,k,ENTR) += 0.5*dt*(S1(i,j,k,ENTR) + S2(i,j,k,ENTR));
      #endif
    #elif TIME_STEPPING == RK3
      #pragma acc loop
      NVAR_LOOP(nv) data->Uc(i,j,k,nv)  =   U0(i,j,k,nv) 
                                          + 0.25*(R1(i,j,k,nv) + R2(i,j,k,nv));
      #pragma acc loop
      EVEC_LOOP(nv) data->Uc(i,j,k,nv) += ci*S1(i,j,k,nv);

      #if ENTROPY_SWITCH
      data->Uc(i,j,k,ENTR) += ci*S1(i,j,k,ENTR);
      #endif
    #elif TIME_STEPPING == RK4
      #pragma acc loop
      NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + R2(i,j,k,nv);
      #pragma acc loop
      EVEC_LOOP(nv) data->Uc(i,j,k,nv) += ci*S2(i,j,k,nv);

      #if ENTROPY_SWITCH
      data->Uc(i,j,k,ENTR) += ci*S2(i,j,k,ENTR);
      #endif
    #endif

  }

  // 3a2. Update staggered variables

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box[:1])
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){

    #if TIME_STEPPING == RK2
    CT_EXPAND(nv, data->Vs[nv](i,j,k) =   Vs0[nv](i,j,k) 
                                        + 0.5*(Rs1[nv](i,j,k) + Rs2[nv](i,j,k));)
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv) data->Vs[EX1s+nv](i,j,k) -= ci*(Jf1(i,j,k,nv) + Jf2(i,j,k,nv));
    #endif

    #elif TIME_STEPPING == RK3

    CT_EXPAND(nv, data->Vs[nv](i,j,k) =   Vs0[nv](i,j,k) 
                                        + 0.25*(Rs1[nv](i,j,k) + Rs2[nv](i,j,k));)
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv) data->Vs[EX1s+nv](i,j,k) -= ci*Jf1(i,j,k,nv); 
    #endif

    #elif TIME_STEPPING == RK4

    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + Rs2[nv](i,j,k);)
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv) data->Vs[EX1s+nv](i,j,k) -= ci*Jf2(i,j,k,nv);
    #endif

    #endif  // TIME_STEPPING
  }
  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // STAGGERED_MHD 

mynvtxstop_();

  // 3a3. Convert to primitive variables

  ConsToPrim3D (data, 0, grid, &box);

  // ------------------------------------------------------
  // 3b. IMEX stage #3 - Implicit contribution:
  //
  //     U3 = U^* + dt*gamma*S3                    [RK3]
  //     U3 = U^* + dt*alpha*S3                    [RK4]
  // ------------------------------------------------------

  #if (TIME_STEPPING == RK3) || (TIME_STEPPING == RK4)

  g_intStage = 30;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  #if   TIME_STEPPING == RK3
  double dt3 = Dts->dt*gIMEX;
  #elif TIME_STEPPING == RK4  
  double dt3 = Dts->dt*alpha;
  #endif

  IMEXImplicitUpdate(data, dt3, S3, Jf3, grid, &boxp);  // d->Uc = U^3 

  // ------------------------------------------------------
  // 4a. IMEX stage #4 - Explicit contribution
  //
  // 
  //    U^(n+1) = U^n + dt*(R1 + R2 + 4*R3)/6
  //                  + dt*(S1 + S2 + 4*S3)/6           [RK3]
  //    U^(n+1) = U^n + dt*(R2 + R3)/4
  //                  + dt*(beta*S1 + eta*S2 + chi*S3)  [RK4]
  // ------------------------------------------------------

  g_intStage = 3;
  // #pragma acc update device(g_intStage)
  // #pragma omp target update to(g_intStage)

  Boundary (data, 0, grid);
  UpdateStage (data, R3, Rs3, grid, Dts, &box, dt);

  // 4a1. Update zone-centered variables 
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box[:1], U0[:1])
  #pragma acc parallel loop collapse(3) present(data, box, U0)
  BOX_LOOP(box, i,j,k){

    #if TIME_STEPPING == RK3
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) = U0(i,j,k,nv) + (      R1(i,j,k,nv) 
                                                       +     R2(i,j,k,nv) 
                                                       + 4.0*R3(i,j,k,nv))/6.0;
    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(      S1(i,j,k,nv) 
                                            +     S2(i,j,k,nv)
                                            + 4.0*S3(i,j,k,nv))/6.0;

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(      S1(i,j,k,ENTR) 
                                +     S2(i,j,k,ENTR)
                                + 4.0*S3(i,j,k,ENTR))/6.0;
    #endif
    
    #elif TIME_STEPPING == RK4
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + 0.25*(R3(i,j,k,nv) + R2(i,j,k,nv));
    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  beta*S1(i,j,k,nv) 
                                            + eta *S2(i,j,k,nv)
                                            + chi *S3(i,j,k,nv));   

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  beta*S1(i,j,k,ENTR)
                                + eta *S2(i,j,k,ENTR)
                                + chi *S3(i,j,k,ENTR));
    #endif 

    #endif
  }

  // 4a2. Update staggered variables 

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box[:1])
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    #pragma acc loop seq
    DIM_LOOP(nv) {
      #if TIME_STEPPING == RK3
      data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + (      Rs1[nv](i,j,k) 
                                              +     Rs2[nv](i,j,k) 
                                              + 4.0*Rs3[nv](i,j,k))/6.0;
      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      int nvE = EX1s+nv;
      double rhs =      (Rs1[nvE](i,j,k) + Rs2[nvE](i,j,k) + 4.0*Rs3[nvE](i,j,k))/6.0
                  - dt*(Jf1(i,j,k,nv) + Jf2(i,j,k,nv) + 4.0*Jf3(i,j,k,nv))/6.0;

      data->Vs[nvE](i,j,k) = Vs0[nvE](i,j,k) + rhs; 
      #endif

      #elif TIME_STEPPING == RK4

      double rhs = 0.25*(Rs2[nv](i,j,k) + Rs3[nv](i,j,k));
      data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + rhs;
      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      int nvE = EX1s+nv;
      rhs = Rs2[nvE](i,j,k) - dt*(  beta*Jf1(i,j,k,nv) 
                                  + eta *Jf2(i,j,k,nv) 
                                  + chi *Jf3(i,j,k,nv))/4.0;
      data->Vs[nvE](i,j,k) = Vs0[nvE](i,j,k) + rhs;
      #endif
      #endif
    }
  }
  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // STAGGERED_MHD 

  // 4a3. Convert to Primitive variables

  ConsToPrim3D (data, 0, grid, &box);

  #endif  // TIME_STEPPING == RK3,RK4 

  // ------------------------------------------------------
  // 4b. IMEX stage #4 - Implicit contribution:
  //
  //     U4 = U^* + dt*alpha*S4                    [RK4]
  // ------------------------------------------------------

  #if TIME_STEPPING == RK4

  g_intStage = 40;
  #pragma acc update device(g_intStage)

  double dt4         = Dts->dt*alpha;
  IMEXImplicitUpdate(data, dt4, S4, Jf4, grid, &boxp);  // d->Uc = U^4 

  // ------------------------------------------------------
  // 5a. IMEX stage #5 - Explicit contribution
  //
  //    U^(n+1) = U^n + dt*(R2 + R3 + 4*R4)/6
  //                  + dt*(S2 + S3 + 4*S4)/6      [RK4]
  // ------------------------------------------------------

  g_intStage = 4;
  #pragma acc update device(g_intStage)

  Boundary (data, 0, grid);
  UpdateStage (data, R4, Rs4, grid, Dts, &box, dt);

  // 5a1. Update zone-centered variables 
  #pragma omp target teams distribute parallel for collapse(4) \
                     map(to: data[:1], box[:1], U0[:1])
  #pragma acc parallel loop collapse(3) present(data, box, U0)
  BOX_LOOP(box, i,j,k){  

    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + (  R2(i,j,k,nv) 
                                                        + R3(i,j,k,nv) 
                                                        + 4.0*R4(i,j,k,nv))/6.0;

    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  S2(i,j,k,nv) 
                                            + S3(i,j,k,nv)
                                            + 4.0*S4(i,j,k,nv))/6.0;   

    #if ENTROPY_SWITCH
    data->Uc(ij,k,ENTR) += dt*(   S2(i,j,k,ENTR)
                                + S3(i,j,k,ENTR)
                                + 4.0*S4(i,j,k,ENTR))/6.0;
    #endif  
  }

  // 5a2. Update staggered variables

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) \
                     map(to: data[:1], box[:1])
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + (      Rs2[nv](i,j,k) 
                                                          +     Rs3[nv](i,j,k) 
                                                          + 4.0*Rs4[nv](i,j,k))/6.0;)
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv){
      data->Vs[EXs+nv](i,j,k) -= dt*(      Jf2(i,j,k,nv) 
                                     +     Jf3(i,j,k,nv)  
                                     + 4.0*Jf4(i,j,k,nv))/6.0;
    }
    #endif
  }
  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif // STAGGERED_MHD 

  // 5a3. Convert to primitive variables

  ConsToPrim3D (data, 0, grid, &box);
  #endif  // TIME_STEPPING == RK4 

  #pragma omp target exit data map (delete: box, box_tot, boxp)
  #pragma acc exit data delete (box, box_tot, boxp)

  mynvtxstop_();
  return 0;
}
