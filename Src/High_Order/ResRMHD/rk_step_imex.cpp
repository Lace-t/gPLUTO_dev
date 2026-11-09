////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Advance equations with IMEX Runge Kutta time integrators.

  Main driver for SSP IMEX Runge-Kutta time evolution.
  Time stepping include:

  RK2:  IMEX-SSP2 (2,2,2) 
  RK3:  IMEX-SSP3 (3,3,2)
  RK4:  IMEX-SSP3 (4,3,3)
  ARK4:
  IERK45:

  - "Implicit-Explicit Runge-Kutta Schemes and Applications to
     Hyperbolic Systems with Relaxation", Journal Sci. Comput (2005)
  - "Beyond ideal MHD: towards a more realistic modelling of relativistic 
     astrophysical plasmas"\n
     Palenzuela \& al, RAS (2009) 394, 1727-1740
  - "A 4th-order scheme for Resistive Relativistic Magnetohydrodynamics"\n
     A. Mignone, V. Berta \& al, MNRAS (2024)
  

  \authors A. Mignone  (andrea.mignone@unito.it)\n
           V. Berta    (vittoria.berta@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it)\n
           
  \date    Jan 14, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define BVEC_LOOP(nv)   for (nv = BX1; nv <= BX3; nv++)
#define EVEC_LOOP(nv)   for (nv = EX1; nv <= EX3; nv++)

#if ENTROPY_SWITCH == YES
  #define ADD_ResRMHD_SOURCES(nv, a)  nv = EX1;  a;  \
                                      nv = EX2;  a;  \
                                      nv = EX3;  a;  \
                                      nv = ENTR; a;
#else
  #define ADD_ResRMHD_SOURCES(nv, a)  nv = EX1;  a;  \
                                      nv = EX2;  a;  \
                                      nv = EX3;  a;
#endif

void IMEXSimpsonSource (Ary4D &, RBox *, Grid *);
void CheckPeriodicity (Data *, Grid *);

////////////////////////////////////////////////////////////////////////
/*!
  Advance the equations by a single time step using unsplit
  integrators based on the method of lines.
  
  \param [in,out]   data  pointer to Data structure
  \param [in,out]    Dts  pointer to time step structure
  \param [in]       grid  pointer to array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
uint32_t AdvanceStep (Data *data, Grid *grid, timeStep *Dts)
{
  int  i, j, k, nv;
  double ci;
  double dt = Dts->dt;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

// --------------------------------------------------------
// Defining Butcher Tableaux for IMEX-SSP2 (2,2,2),
// IMEX-SSP3 (3,3,2), IMEX-SSP3 (4,3,3)  ARK4(3)7L[2]SA1 and
// IERK45. 
// Referring to PR05  a = implicit coeff., b = explicit coeff.       
// --------------------------------------------------------

  const double gIMEX = (1.0 - 1.0/sqrt(2.0));

  #if TIME_STEPPING == RK2
  const double a11 = gIMEX; 
  const double a21 = 1.0 - 2.0*gIMEX, a22 = gIMEX;
  const double a31 = 0.5,             a32 = 0.5;
  const double b21 = 1.0; 
  const double b31 = 0.5,             b32 = 0.5;
  #elif TIME_STEPPING == RK3
  const double a11 = gIMEX;
  const double a21 = 1.0 - 2.0*gIMEX, a22 = gIMEX;
  const double a31 = 0.5 - gIMEX,                    a33 = gIMEX;
  const double a41 = 1.0/6.0,         a42 = 1.0/6.0, a43 = 2.0/3.0;
  const double b21 = 1.0; 
  const double b31 = 0.25,            b32 = 0.25;
  const double b41 = 1.0/6.0,         b42 = 1.0/6.0, b43 = 2.0/3.0;
  #elif TIME_STEPPING == RK4
  const double alpha = 0.24169426078821;
  const double beta  = 0.06042356519705;
  const double eta   = 0.12915286960590;
  const double chi   = 0.5-beta-eta-alpha;
  const double a11 =  alpha;
  const double a21 = -alpha, a22 = alpha;
  const double               a32 = 1.0 - alpha, a33 = alpha;
  const double a41 = beta,   a42 = eta,         a43 = chi, a44 = alpha;
  const double               a52 = 1.0/6.0, a53 = 1.0/6.0, a54 = 2.0/3.0; 
  const double b21 = 1.0;
  const double               b32 = 1.0;
  const double               b42 = 0.25, b43 = 0.25;
  const double               b52 = 1.0/6.0, b53 = 1.0/6.0, b54 = 2.0/3.0; 
  #elif TIME_STEPPING == ARK4
  const double kappa = 0.1235;
  const double a21 = kappa,          a22 = kappa;
  const double a31 = 0.14907768747,  a32 = 0.14907768747,  a33 = kappa;
  const double a41 = 0.12483442871,  a42 = 0.12483442871,  a43 = -0.03816885743;
  const double a44 = kappa;
  const double a51 = -0.0730319403,  a52 = -0.0730319403,  a53 = -0.24343568716;
  const double a54 = 0.34099956776,  a55 = kappa;
  const double a61 = -0.15296500088, a62 = -0.15296500088, a63 = 0.07220562047; 
  const double a64 = 0.40430630248,  a65 = 0.40591807880,  a66 = kappa;
  const double a73 = 0.51611072831,  a74 = -0.14606356393, a75 = 0.23473048589;
  const double a76 = 0.27172234973,  a77 = kappa;

  const double b21 = 0.247;
  const double b31 = 0.06175,        b32 = 0.35990537495;
  const double b41 = 0.05301658458,  b42 = 0.35949264529,  b43 = -0.07750922988;
  const double b51 = 0.05841715944,  b52 = -0.16313824817, b53 = -0.19732090979; 
  const double b54 = 0.37704199852;
  const double b61 = 0.53853032270,  b62 = -0.45497746895, b63 = 1.25629056234; 
  const double b64 = -0.47828452721, b65 = -0.16155888888;
  const double b71 = 0.23221715782,  b72 = 0.23221715782,  b73 = -6.80999437504;
  const double b74 = 7.36185855242,  b75 = -1.37487907794, b76 = 1.35858058491;

  const double c3  = 0.51611072831,  c4  = -0.14606356393, c5  = 0.23473048589;
  const double c6  = 0.27172234973,  c7  = kappa;
  #elif TIME_STEPPING == IERK45
  const double iota  = 0.25;
  const double a11 =  iota;
  const double a21 =  0.34114705729739, a22 =  iota;
  const double a31 =  0.80458720789763, a32 = -0.07095262154540, a33 = iota;
  const double a41 = -0.52932607329103, a42 =  1.15137638494253;
  const double a43 = -0.80248263237803, a44 =  iota;
  const double a51 =  0.11933093090075, a52 =  0.55125531344927;
  const double a53 = -0.1216872844994 , a54 =  0.20110104014943, a55 = iota;

  const double b21 =  0.39098372452428;
  const double b31 =  1.09436646160460, b32 =  0.33181504274704;
  const double b41 =  0.14631668003312, b42 =  0.69488738277516; 
  const double b43 =  0.46893381306619;
  const double b51 = -1.33389883143642, b52 =  2.90509214801204;
  const double b53 = -1.06511748457024, b54 =  0.27210900509137;

  const double c1 = a51, c2 = a52, c3 = a53, c4 = a54, c5 = a55;
  #endif

  RBox   box, box_tot;

  const  Ary4D &U0 = data->Uc0;
  static Ary4D Rc;
  static Ary4D S1, S2, S3, S4, S5, S6, S7;
  static Ary4D R1, R2, R3, R4, R5, R6, R7;
  static Ary4D Jf1, Jf2, Jf3, Jf4, Jf5, Jf6, Jf7;
  static Ary3D Vs0[6], Rs[6];
  static Ary3D Rs1[6],  Rs2[6], Rs3[6], Rs4[6], Rs5[6], Rs6[6], Rs7[6];
  // #if FALLBACK == YES
  // static Ary3D RsA1[3], RsA2[3], RsA3[3], RsA4[3], RsA5[3];
  // #endif

  // --------------------------------------------------------
  // 0a. Initialize boxes, allocate memory
  // --------------------------------------------------------

  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &box_tot);
  #pragma omp target enter data map(to: box, box_tot)
  #pragma acc enter data copyin (box, box_tot)

  if (Rc.isAllocated() == 0) {

    Rc.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R1.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R2.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R3.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R4.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R5.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    #if TIME_STEPPING == ARK4
    R6.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    R7.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    #endif

    S1.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S2.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S3.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S4.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S5.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    #if TIME_STEPPING == ARK4
    S6.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    S7.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    #endif 
    #pragma omp target enter data map(to: Rc, R1, R2, R3, R4, R5, R6, R7)
    #pragma omp target enter data map(to: S1, S2, S3, S4, S5, S6, S7)
    #pragma acc enter data copyin(Rc)
    #pragma acc enter data copyin(R1, R2, R3, R4, R5, R6, R7)
    #pragma acc enter data copyin(S1, S2, S3, S4, S5, S6, S7)

    #ifdef STAGGERED_MHD
    Jf1.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf2.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf3.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf4.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf5.create (nx1_tot,nx2_tot,nx3_tot,3);
    #if TIME_STEPPING == ARK4
    Jf6.create (nx1_tot,nx2_tot,nx3_tot,3);
    Jf7.create (nx1_tot,nx2_tot,nx3_tot,3);
    #endif

    int (*ns)[3];
    ns = grid->ntot_stag;

    CT_EXPAND(nv, Vs0[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                   Rs[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs1[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs2[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs3[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs4[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs5[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs6[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
                  Rs7[nv].create(ns[nv][IDIR], ns[nv][JDIR], ns[nv][KDIR]);
    )

    #pragma omp target enter data map(to: Vs0, Rs)
    #pragma omp target enter data map(to: Jf1, Jf2, Jf3, Jf4, Jf5, Jf6, Jf7)
    #pragma omp target enter data map(to: Rs1, Rs2, Rs3, Rs4, Rs5, Rs6, Rs7)
    
    #pragma acc enter data copyin(Vs0, Rs)
    #pragma acc enter data copyin(Jf1, Jf2, Jf3, Jf4, Jf5, Jf6, Jf7)
    #pragma acc enter data copyin(Rs1, Rs2, Rs3, Rs4, Rs5, Rs6, Rs7)
    #endif // STAGGERED_MHD 
  }

mynvtxstart_("RKSTEP",GREEN);

  // ------------------------------------------------------
  // 0b. Reset flag to zero and evaluate trouble cells
  //     of the new time step
  // ------------------------------------------------------

  g_intStage = 10;
  // #pragma acc update device(g_intStage)

  ConsToPrim3D (data, 0, grid, &box);
  Boundary (data, 0, grid);

  #if (SHOCK_FLATTENING == MULTID)
  FlagShock (data, grid, &box);
  #endif

  #if (HO_LAP_LIMITER != NO) || (HO_ORDER_REDUCTION != NO)
  DiscontinuityDetect(data, grid, &box);
  #endif

  // ------------------------------------------------------
  // 0c. Copy conservative variables Uc on U0 to start 
  // ------------------------------------------------------
  
  #pragma omp target teams distribute parallel for collapse(4) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(4) present(data, box)
  BOX_LOOP(box,i,j,k) NVAR_LOOP(nv) U0(i,j,k,nv) = data->Uc(i,j,k,nv);

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box_tot)
  BOX_LOOP(box_tot, i,j,k) {
    CT_EXPAND(nv, Vs0[nv](i,j,k) = data->Vs[nv](i,j,k);)
    // #if FALLBACK == YES      //  Reset vector potential           
    // DIM_EXPAND(                              ,
    //            data->A[KDIR](i,j,k) = 0.0;   ,
    //            data->A[IDIR](i,j,k) = 0.0;  
    //            data->A[JDIR](i,j,k) = 0.0;)
    // #endif
  }
  #endif // STAGGERED_MHD 

  // ------------------------------------------------------
  // 1. IMEX Stage #1 - Implicit contribution:
  //
  //    U(1) = U(n) + a11*dt*S1  [SSP(2,2,2), SSP(3,3,2), SSP(4,3,3), IERK45]
  //    U(1) = U(n)              [ARK4] (No implicit 1st step)  
  // ------------------------------------------------------

  #ifdef HIGH_ORDER
  // BoundaryConservative (data, grid);     // Set b.c. on data->Uc (1 ghost)  
  // PointValue(data, grid, &box);          // Obtain point values Upc, Vpc 
  Boundary (data, 1, grid);          
  #else

  #endif

  #if TIME_STEPPING != ARK4
  IMEXImplicitUpdate(data, dt*a11, S1, Jf1, grid, &box_tot); // data->Uc = U^1
  #elif TIME_STEPPING == ARK4
  IMEXExplicitSource (data, S1, grid);           // Compute source term      
  #endif

  #ifdef HIGH_ORDER
  IMEXSimpsonSource  (S1, &box, grid);    // High-order source term evaluation
  #endif

  // ------------------------------------------------------
  // 2a. IMEX stage #2 - Explicit contribution:
  //
  //     U(2*) = U(n) + dt*(b21*R1 + a21*S1)  [SSP(2,2,2), SSP(3,3,2), 
  //                                           ARK4,       IERK45]
  //     U(2*) = U(n) - alpha*dt*S1           [SSP(4,3,3)]
  // ------------------------------------------------------

  g_intStage = 1;
  // #pragma acc update device(g_intStage)

  #if TIME_STEPPING != RK4
  UpdateStage (data, R1, Rs1, grid, Dts, &box, dt);
  #elif TIME_STEPPING == RK4
  UpdateStage (data, R1, Rs1, grid, Dts, &box, 0.0);

  #pragma omp target teams distribute parallel for collapse(4) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(4) present(box_tot)
  NVAR_LOOP(nv) BOX_LOOP(box_tot, i,j,k) R1(i,j,k,nv) = 0.0;
    
  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(box_tot)
  BOX_LOOP(box_tot, i,j,k) {
    CT_EXPAND(nv, Rs1[nv](i,j,k) = 0.0;)
  }
  #endif // STAGGERED_MHD 
  #endif // TIME_STEPPING != RK4 

  mynvtxstart_("Stage#1",BLUE);

  // 2a1. Update zone-centered variables
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i,j,k){
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b21*R1(i,j,k,nv);
    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*a21*S1(i,j,k,nv);

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*a21*S1(i,j,k,ENTR);
    #endif
  }

  // 2a2. Update staggered variables

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b21*Rs1[nv](i,j,k);)

    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv) data->Vs[EX1s+nv](i,j,k) -= dt*a21*Jf1(i,j,k,nv);
    #endif 
  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

  mynvtxstop_();

  // ------------------------------------------------------
  // 2b. IMEX stage #2 (Implicit part):
  //
  //     U(2) = U(2*) + dt*a22*S2    [SSP(2,2,2), SSP(3,3,2), 
  //                                  SSP(4,3,3), ARK4, IERK45]
  // ------------------------------------------------------

  g_intStage = 20;
  // #pragma acc update device(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);            // Set b.c. on data->Vpc   
  IMEXImplicitUpdate(data, dt*a22, S2, Jf2, grid, &box_tot);  // data->Uc = U^2 
  IMEXSimpsonSource  (S2, &box, grid);  // High-order source term evaluation 
  #else
  IMEXImplicitUpdate(data, dt*a22, S2, Jf2, grid, &box);      // data->Uc = U^2 
  #endif

  // ------------------------------------------------------
  // 3a. IMEX stage #3 - Explicit contribution
  //
  //     U(n+1) = U(n) + dt*(R1  + R2)/2 + dt*(S1 + S2)/2     [SSP(2,2,2)]
  //     U(3*)  = U(n) + dt*(R2) + dt*(1-alpha)*S2            [SSP(4,3,3)]
  //     U(3*)  = U(n) + dt*(b31*R1 + b32*R2)                 [SSP(3,3,2), 
  //                   + dt*(a31*S1 + a32*S2)                  ARK4, IERK45]
  // ------------------------------------------------------

  g_intStage = 2;
  // #pragma acc update device(g_intStage)

  UpdateStage (data, R2, Rs2, grid, Dts, &box, dt);

  // 3a1. Update zone-centered variables 

  mynvtxstart_("Stage#1",BLUE);
    
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i,j,k){

    #if TIME_STEPPING == RK4 // -- Optimized version -- 
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b32*R2(i,j,k,nv);
    #else
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b31*R1(i,j,k,nv) 
                                                     + b32*R2(i,j,k,nv);
    #endif

    // Optimized summation since some coeffs are == 0.0 
    
    #if TIME_STEPPING == RK3

    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*a31*S1(i,j,k,nv);

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*a31*S1(i,j,k,ENTR);
    #endif

    #elif TIME_STEPPING == RK4
    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*a32*S2(i,j,k,nv);

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*a32*S2(i,j,k,ENTR);
    #endif
    #else
    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  a31*S1(i,j,k,nv) 
                                            + a32*S2(i,j,k,nv));

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  a31*S1(i,j,k,ENTR) 
                                + a32*S2(i,j,k,ENTR));
    #endif
    #endif
  }
  
  // 3a2. Update staggered variables

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present (data,box)
  BOX_STAG_LOOP(box, i,j,k){
    #if TIME_STEPPING == RK4   // -- Optimized version -- 
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b32*Rs2[nv](i,j,k);) 
    #else
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b31*Rs1[nv](i,j,k) 
                                                       + b32*Rs2[nv](i,j,k);)
    #endif

    // #if FALLBACK  == YES
    // DIM_EXPAND(                                                        ,
    //            data->A[KDIR](i,j,k) += b52*dt*data->emf->Ex3e(i,j,k);  ,
    //            data->A[IDIR](i,j,k) += b52*dt*data->emf->Ex1e(i,j,k);
    //            data->A[JDIR](i,j,k) += b52*dt*data->emf->Ex2e(i,j,k); )
    // #endif
 
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv) {
      #if TIME_STEPPING == RK3
      data->Vs[EX1s + nv](i,j,k) -= dt*a31*Jf1(i,j,k,nv);
      #elif TIME_STEPPING == RK4
      data->Vs[EX1s + nv](i,j,k) -= dt*a32*Jf2(i,j,k,nv);
      #else
      data->Vs[EX1s + nv](i,j,k) -= dt*(  a31*Jf1(i,j,k,nv) 
                                        + a32*Jf2(i,j,k,nv));
      #endif
    }
    #endif
  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

  mynvtxstop_();

#if (TIME_STEPPING != RK2)

  // ------------------------------------------------------
  // 3b. IMEX stage #3 - Implicit contribution
  //
  //     U(3) = U(3*) + dt*a33*S3    [SSP(3,3,2), SSP(4,3,3), ARK4, IERK45]
  // ------------------------------------------------------

  g_intStage = 30;
  // #pragma acc update device(g_intStage)

  #ifdef HIGH_ORDER  
  Boundary (data, 1, grid);                    // Set b.c. on data->Vpc    
  IMEXImplicitUpdate(data, dt*a33, S3, Jf3, grid, &box_tot);  // data->Uc = U^3
  IMEXSimpsonSource  (S3, &box, grid);  // High-order source term evaluation 
  #else
  IMEXImplicitUpdate(data, dt*a33, S3, Jf3, grid, &box);  // data->Uc = U^3
  #endif

  // ------------------------------------------------------
  // 4a. IMEX stage #4 - Explicit contribution
  //
  //     U(n+1) = Un + dt*(R1 + R2 + 4*R3)/6
  //                 + dt*(S1 + S2 + 4*S3)/6           [SSP(3,3,2)]
  //     U(4*)  = Un + dt*(R2 + R3)/4
  //                 + dt*(beta*S1 + eta*S2 + chi*S3)  [SSP(4,3,3)]
  //     U(4*)  = Un + dt*(b41*R1 + b42*R2 + b43*R3)
  //                 + dt*(a41*S1 + a42*S2 + a43*S3)   [ARK4]
  // ------------------------------------------------------

  g_intStage = 3;
  // #pragma acc update device(g_intStage)

  UpdateStage (data, R3, Rs3, grid, Dts, &box, dt);

  // 4a1. Update zone-centered variables 

  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present( data, box)
  BOX_LOOP(box, i,j,k){

    #if TIME_STEPPING == RK4
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) +
                                                     + b42*R2(i,j,k,nv) 
                                                     + b43*R3(i,j,k,nv);
    #else
    #pragma acc loop
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b41*R1(i,j,k,nv)
                                                     + b42*R2(i,j,k,nv) 
                                                     + b43*R3(i,j,k,nv);
    #endif

//ADD_ResRMHD_SOURCES(nv, data->Uc(i,j,k,nv) += dt*(...););
    #pragma acc loop
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  a41*S1(i,j,k,nv) 
                                            + a42*S2(i,j,k,nv)
                                            + a43*S3(i,j,k,nv));   

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  a41*S1(i,j,k,ENTR) 
                                + a42*S2(i,j,k,ENTR)
                                + a43*S3(i,j,k,ENTR));   
    #endif                                        
  }

  // 4a2. Update staggered variables 

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    #if TIME_STEPPING == RK4 // Optimized version 
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b42*Rs2[nv](i,j,k)
                                                       + b43*Rs3[nv](i,j,k);)
    #else
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b41*Rs1[nv](i,j,k)
                                                       + b42*Rs2[nv](i,j,k)
                                                       + b43*Rs3[nv](i,j,k);)
    #endif

    // #if FALLBACK  == YES
    // DIM_EXPAND(                                                    ,
    //            data->A[KDIR](i,j,k) += b53*dt*data->emf->Ex3e(i,j,k);  ,
    //            data->A[IDIR](i,j,k) += b53*dt*data->emf->Ex1e(i,j,k);
    //            data->A[JDIR](i,j,k) += b53*dt*data->emf->Ex2e(i,j,k); )
    // #endif

    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    #pragma acc loop seq
    DIM_LOOP(nv) data->Vs[EX1s+nv](i,j,k) -= dt*(  a41*Jf1(i,j,k,nv)
                                                 + a42*Jf2(i,j,k,nv)
                                                 + a43*Jf3(i,j,k,nv));
    #endif // DIVE_CONTROL == CONSTRAINED_TRANSPORT
  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

#endif  // TIME_STEPPING != RK2

#if (TIME_STEPPING == RK4) || (TIME_STEPPING == ARK4) || (TIME_STEPPING == IERK45)

  // ------------------------------------------------------
  // 4b. IMEX stage #4 - Implicit contribution
  //
  //     U(4) = U(4*) + dt*a44*S4       [SSP(4,3,3), ARK4, IERK45]
  // ------------------------------------------------------

  g_intStage = 40;
  // #pragma acc update device(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);            // Set BC on data->Vpc     
  IMEXImplicitUpdate(data, dt*a44, S4, Jf4, grid, &box_tot);  // data->Uc = U^4 
  IMEXSimpsonSource  (S4, &box, grid); // High-order source term evaluation 
  #else
  IMEXImplicitUpdate(data, dt*a44, S4, Jf4, grid, &box);  // data->Uc = U^4
  #endif

  // ------------------------------------------------------
  // 5a. IMEX stage #5 - Explicit contribution
  //
  //    U(n+1) = Un + dt*(R2 + R3 + 4*R4)/6
  //                + dt*(S2 + S3 + 4*S4)/6                   [SSP(4,3,3)]
  //    U(5*)  = Un + dt*(b51*R1 + b52*R2 + b53*R3 + b54*R4)
  //                + dt*(a51*S1 + a52*S2 + a53*S3 + a54*S4)  [ARK4, IERK45]
  // ------------------------------------------------------

  g_intStage = 4;
  // #pragma acc update device(g_intStage)

  UpdateStage (data, R4, Rs4, grid, Dts, &box, dt); // Add to current rhs 

  // 5a1. Update zone-centered variables 
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i,j,k){  

    #if TIME_STEPPING == RK4 
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) 
                                                     + b52*R2(i,j,k,nv) 
                                                     + b53*R3(i,j,k,nv) 
                                                     + b54*R4(i,j,k,nv);
    
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  
                                              a52*S2(i,j,k,nv) 
                                            + a53*S3(i,j,k,nv)
                                            + a54*S4(i,j,k,nv));

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  
                                  a52*S2(i,j,k,ENTR) 
                                + a53*S3(i,j,k,ENTR)
                                + a54*S4(i,j,k,ENTR));
    #endif
    #else
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b51*R1(i,j,k,nv) 
                                                     + b52*R2(i,j,k,nv) 
                                                     + b53*R3(i,j,k,nv) 
                                                     + b54*R4(i,j,k,nv);
    
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  a51*S1(i,j,k,nv)
                                            + a52*S2(i,j,k,nv) 
                                            + a53*S3(i,j,k,nv)
                                            + a54*S4(i,j,k,nv));  

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  a51*S1(i,j,k,ENTR)  
                                + a52*S2(i,j,k,ENTR) 
                                + a53*S3(i,j,k,ENTR)
                                + a54*S4(i,j,k,ENTR));
    #endif

    #endif  
  }

  // 5a2. Update staggered variables 

  #ifdef STAGGERED_MHD 
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    #if TIME_STEPPING == RK4
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b52*Rs2[nv](i,j,k)
                                                       + b53*Rs3[nv](i,j,k)
                                                       + b54*Rs4[nv](i,j,k);)

    #else
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b51*Rs1[nv](i,j,k) 
                                                       + b52*Rs2[nv](i,j,k)
                                                       + b53*Rs3[nv](i,j,k) 
                                                       + b54*Rs4[nv](i,j,k);)
    #endif // TIME_STEPPING  

    // #if FALLBACK  == YES
    // DIM_EXPAND(                                                    ,
    //            data->A[KDIR](i,j,k) += b54*dt*data->emf->Ex3e(i,j,k);  ,
    //            data->A[IDIR](i,j,k) += b54*dt*data->emf->Ex1e(i,j,k);
    //            data->A[JDIR](i,j,k) += b54*dt*data->emf->Ex2e(i,j,k); )
    // #endif

    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT 
      #error DIVE_CONTROL == CONSTRAINED_TRANSPORT not permitted with high order
    #endif
  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

#endif //  (TIME_STEPPING == RK4) || (TIME_STEPPING == ARK4)  || (TIME_STEPPING == IERK45)
  
#if (TIME_STEPPING == ARK4) || (TIME_STEPPING == IERK45)

  // ------------------------------------------------------
  // 5b. IMEX stage #5 - Implicit contribution
  //
  //     U(5) = U(5*) + dt*a55*S5       [ARK4, IERK45]
  // ------------------------------------------------------

  g_intStage = 50;
  // #pragma acc update device(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);            // Set b.c. on data->Vpc
  IMEXImplicitUpdate(data, dt*a55, S5, Jf5, grid, &box_tot);  // data->Uc = U^5
  IMEXSimpsonSource  (S5, &box, grid);  // High-order source term evaluation 
  #else
  IMEXImplicitUpdate(data, dt*a55, S5, Jf5, grid, &box);  // data->Uc = U^5
  #endif  

  // ------------------------------------------------------
  // 6a. IMEX stage #6 - Explicit contribution
  //
  //     U(n+1) = Un + dt*(c1*R1 + c2*R2 + c3*R3 + c4*R4 + c5*R5)       [IERK45]
  //                 + dt*(c1*S1 + c2*S2 + c3*S3 + c4*S4 + c5*S5)
  //     U(6*)  = Un + dt*(b61*R1 + b62*R2 + b63*R3 + b64*R4 + b65*R5)
  //                 + dt*(a61*S1 + a62*S2 + a63*S3 + a64*S4 + a65*S5)  [ARK4]
  // ------------------------------------------------------
  
  g_intStage = 5;
  // #pragma acc update device(g_intStage)

  UpdateStage (data, R5, Rs5, grid, Dts, &box, dt); 

  // 6a1. Update zone-centered variables
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i,j,k){  

    #if TIME_STEPPING == IERK45 
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + c1*R1(i,j,k,nv) 
                                                     + c2*R2(i,j,k,nv) 
                                                     + c3*R3(i,j,k,nv) 
                                                     + c4*R4(i,j,k,nv)
                                                     + c5*R5(i,j,k,nv);
    
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  c1*S1(i,j,k,nv) 
                                            + c2*S2(i,j,k,nv) 
                                            + c3*S3(i,j,k,nv)
                                            + c4*S4(i,j,k,nv)
                                            + c5*S5(i,j,k,nv));

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  c1*S1(i,j,k,ENTR) 
                                + c2*S2(i,j,k,ENTR) 
                                + c3*S3(i,j,k,ENTR)
                                + c4*S4(i,j,k,ENTR)
                                + c5*S5(i,j,k,ENTR));
    #endif

    #elif TIME_STEPPING == ARK4
    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b61*R1(i,j,k,nv) 
                                                     + b62*R2(i,j,k,nv) 
                                                     + b63*R3(i,j,k,nv) 
                                                     + b64*R4(i,j,k,nv)
                                                     + b65*R5(i,j,k,nv);
    
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(  a61*S1(i,j,k,nv) 
                                            + a62*S2(i,j,k,nv) 
                                            + a63*S3(i,j,k,nv)
                                            + a64*S4(i,j,k,nv)
                                            + a65*S5(i,j,k,nv));  

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  a61*S1(i,j,k,ENTR) 
                                + a62*S2(i,j,k,ENTR) 
                                + a63*S3(i,j,k,ENTR)
                                + a64*S4(i,j,k,ENTR)
                                + a65*S5(i,j,k,ENTR));
    #endif

    #endif  
  }

  // 6a2. Update staggered variables 

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    #if TIME_STEPPING == IERK45
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + c1*Rs1[nv](i,j,k)
                                                       + c2*Rs2[nv](i,j,k)
                                                       + c3*Rs3[nv](i,j,k)
                                                       + c4*Rs4[nv](i,j,k)
                                                       + c5*Rs5[nv](i,j,k);)
    #elif TIME_STEPPING == ARK4
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b61*Rs1[nv](i,j,k)
                                                       + b62*Rs2[nv](i,j,k)
                                                       + b63*Rs3[nv](i,j,k)
                                                       + b64*Rs4[nv](i,j,k)
                                                       + b65*Rs5[nv](i,j,k);)
    #endif

    // #if FALLBACK  == YES
    // DIM_EXPAND(                                                    ,
    //            data->A[KDIR](i,j,k) += b65*dt*data->emf->Ex3e(i,j,k);  ,
    //            data->A[IDIR](i,j,k) += b65*dt*data->emf->Ex1e(i,j,k);
    //            data->A[JDIR](i,j,k) += b65*dt*data->emf->Ex2e(i,j,k); )
    // #endif
  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

#endif  //  (TIME_STEPPING == ARK4)  || (TIME_STEPPING == IERK45) 

#if (TIME_STEPPING == ARK4)

  // ------------------------------------------------------
  // 6b. IMEX stage #6 - Explicit contribution
  //
  //     U(6) = U(6*) + dt*a66*S6       [ARK4]
  // ------------------------------------------------------
  
  g_intStage = 60;
  // #pragma acc update device(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);            // Set b.c. on data->Vpc     
  IMEXImplicitUpdate(data, dt*a66, S6, Jf6, grid, &box_tot);  // data->Uc = U^6
  IMEXSimpsonSource  (S6, &box, grid);  // High-order source term evaluation 
  #else
  IMEXImplicitUpdate(data, dt*a66, S6, Jf6, grid, &box);  // data->Uc = U^6
  #endif

  // ------------------------------------------------------
  // 7a. IMEX stage #7 - Explicit contribution [ARK4]
  //
  //    U(6*)  = Un + dt*(b71*R1 + b72*R2 + b73*R3 + b74*R4 + b75*R5 + b76*R6)
  //                + dt*(a73*S3 + a74*S4 + a75*S5 + a76*S6)  
  // ------------------------------------------------------

  g_intStage = 6;
  // #pragma acc update device(g_intStage)

  UpdateStage (data, R6, Rs6, grid, Dts, &box, dt); // Add to current rhs 

  // 7a1. Update zone-centered variables
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present( data, box)
  BOX_LOOP(box, i,j,k){  

    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) + b71*R1(i,j,k,nv) 
                                                     + b72*R2(i,j,k,nv) 
                                                     + b73*R3(i,j,k,nv) 
                                                     + b74*R4(i,j,k,nv)
                                                     + b75*R5(i,j,k,nv)
                                                     + b76*R6(i,j,k,nv);
    
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(
                                            + a73*S3(i,j,k,nv)
                                            + a74*S4(i,j,k,nv)
                                            + a75*S5(i,j,k,nv)
                                            + a76*S6(i,j,k,nv));  

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  a73*S3(i,j,k,ENTR)
                                + a74*S4(i,j,k,ENTR)
                                + a75*S5(i,j,k,ENTR)
                                + a76*S6(i,j,k,ENTR));
    #endif
  }

  // 7a2. Update staggered variables 

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data)
  BOX_STAG_LOOP(box, i,j,k){
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + b71*Rs1[nv](i,j,k)
                                                       + b72*Rs2[nv](i,j,k)
                                                       + b73*Rs3[nv](i,j,k)
                                                       + b74*Rs4[nv](i,j,k)
                                                       + b75*Rs5[nv](i,j,k)
                                                       + b76*Rs6[nv](i,j,k);)

    // #if FALLBACK  == YES
    // DIM_EXPAND(                                                    ,
    //            data->A[KDIR](i,j,k) += b76*dt*data->emf->Ex3e(i,j,k);  ,
    //            data->A[IDIR](i,j,k) += b76*dt*data->emf->Ex1e(i,j,k);
    //            data->A[JDIR](i,j,k) += b76*dt*data->emf->Ex2e(i,j,k); )
    // #endif

  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

  // ------------------------------------------------------
  // 7b. IMEX stage #7 - Implicit contribution [ARK4]
  //
  //     U(7) = U(7*) + dt*a77*S7       [ARK4]
  // ------------------------------------------------------
  
  g_intStage = 70;
  // #pragma acc update device(g_intStage)

  #ifdef HIGH_ORDER
  Boundary (data, 1, grid);            // Set b.c. on data->Vpc        
  IMEXImplicitUpdate(data, dt*a77, S7, Jf7, grid, &box_tot);  // data->Uc = U7 
  IMEXSimpsonSource  (S7, &box, grid);  // High-order source term evaluation 
  #else
  IMEXImplicitUpdate(data, dt*a77, S7, Jf7, grid, &box);  // data->Uc = U7
  #endif

  // ------------------------------------------------------
  // 8a. IMEX stage #8 - Explicit contribution [ARK4]
  //
  //
  //     U(6*)  = Un + dt*(b71*R1 + b72*R2 + b73*R3 + b74*R4 + b75*R5 + b76*R6)
  //                 + dt*(a73*S3 + a74*S4 + a75*S5 + a76*S6)  
  // ------------------------------------------------------

  g_intStage = 7;
  // #pragma acc update device(g_intStage)

  UpdateStage (data, R7, Rs7, grid, Dts, &box, dt); 

  // 8a1. Update zone-centered variables
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(box, i,j,k){  

    NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = U0(i,j,k,nv) 
                                                     + c3*R3(i,j,k,nv) 
                                                     + c4*R4(i,j,k,nv)
                                                     + c5*R5(i,j,k,nv)
                                                     + c6*R6(i,j,k,nv)
                                                     + c7*R7(i,j,k,nv);
    
    EVEC_LOOP(nv) data->Uc(i,j,k,nv) += dt*(
                                            + c3*S3(i,j,k,nv)
                                            + c4*S4(i,j,k,nv)
                                            + c5*S5(i,j,k,nv)
                                            + c6*S6(i,j,k,nv)  
                                            + c7*S7(i,j,k,nv));  

    #if ENTROPY_SWITCH
    data->Uc(i,j,k,ENTR) += dt*(  c3*S3(i,j,k,ENTR)
                                + c4*S4(i,j,k,ENTR)
                                + c5*S5(i,j,k,ENTR)
                                + c6*S6(i,j,k,ENTR)
                                + c7*S7(i,j,k,ENTR));
    #endif
  }

  // 8a2. Update staggered variables

  #ifdef STAGGERED_MHD
  #pragma omp target teams distribute parallel for collapse(3) map ( to: data[0:1] )
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_STAG_LOOP(box, i,j,k){
    CT_EXPAND(nv, data->Vs[nv](i,j,k) = Vs0[nv](i,j,k) + c3*Rs3[nv](i,j,k)
                                                       + c4*Rs4[nv](i,j,k)
                                                       + c5*Rs5[nv](i,j,k)
                                                       + c6*Rs6[nv](i,j,k)
                                                       + c7*Rs7[nv](i,j,k);)

    // #if FALLBACK  == YES
    // DIM_EXPAND(                                                       ,
    //            data->A[KDIR](i,j,k) += c7*dt*data->emf->Ex3e(i,j,k);  ,
    //            data->A[IDIR](i,j,k) += c7*dt*data->emf->Ex1e(i,j,k);
    //            data->A[JDIR](i,j,k) += c7*dt*data->emf->Ex2e(i,j,k); )
    // #endif
  }

  CT_AverageStaggeredFields (data, 1, &box, grid);
  #endif

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);   // Set b.c. on data->Uc (1 ghost)  
  PointValue(data, grid, &box);        // Obtain point values Upc, Vpc 
  #endif

#endif // TIME_STEPPING == ARK4 

  ConsToPrim3D (data, 0, grid, &box);

  #pragma acc exit data delete (box, box_tot)
  
  mynvtxstop_();

  return 0;
}

////////////////////////////////////////////////////////////////////////
void IMEXSimpsonSource(Ary4D &S, RBox *box, Grid *grid)
/*! 
  Compute the volume average of the source term. 
  On output S will be replace with its volume average.
*/
////////////////////////////////////////////////////////////////////////
{
  int i,j,k,nv;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  static Ary4D d2S;

  if (d2S.isAllocated() == 0) {
    d2S.create (nx1_tot,nx2_tot,nx3_tot,NVAR);
    #pragma omp target enter data map(to: d2S)
    #pragma acc enter data copyin(d2S)
  }

  // ------------------------------------------------------
  // 1. Compute Laplacian of source term 
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(4) 
  #pragma acc parallel loop collapse(4) present(S, d2S, box)
  NVAR_LOOP (nv) BOX_LOOP(*box, i,j,k){
    // d2S(i,j,k,nv) =   HO_LAPLXc(S,i,j,k,nv) 
    //                 + HO_LAPLYc(S,i,j,k,nv) 
    //                 + HO_LAPLZc(S,i,j,k,nv);
    d2S(i,j,k,nv) = 0.0;
    DIM_EXPAND(d2S(i,j,k,nv) += HO_Average4D<IDIR>(S, i,j,k,nv);  ,
               d2S(i,j,k,nv) += HO_Average4D<JDIR>(S, i,j,k,nv);  ,
               d2S(i,j,k,nv) += HO_Average4D<KDIR>(S, i,j,k,nv);)
  }

  // ------------------------------------------------------
  // 2. Modify input source term
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(4)
  #pragma acc parallel loop collapse(4) present(S, d2S, box)
  NVAR_LOOP (nv) BOX_LOOP(*box, i,j,k){
    S(i,j,k,nv) += d2S(i,j,k,nv)/24.0;
  }
}
