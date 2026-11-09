////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Loop on the computational cells to assign initial conditions.

  This function is called anyway, even if restart from file is enabled.
  This is useful to initialized a number of global variables and/or
  user-defined parameters by calling Init().

  \author A. Mignone  (andrea.mignone@unito.it)
          A. Suriano  (alession.suriano@unito.it)
          M. Rossazza (marco.rossazza@unito.it)
          V. Berta    (vittoria.berta@unito.it)

  \date   Jan 15, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifdef HIGH_ORDER
  #define Ng  HO_STARTUP_NGAUSS
#else 
  #define Ng 1
#endif

#if (ASSIGN_VECTOR_POTENTIAL == YES) && !(defined STAGGERED_MHD)
void FieldAverage(Ary3D *, Ary4D &, Grid *);
#endif

////////////////////////////////////////////////////////////////////////
/*!
  Assign initial conditions
*/
////////////////////////////////////////////////////////////////////////
void Startup (Data *data, Grid *grid)
{
  int i, j, k, nv;
  int nx1_tot=grid->ntot[IDIR];
  int nx2_tot=grid->ntot[JDIR];
  int nx3_tot=grid->ntot[KDIR];
  double x1,  x2,  x3;
  RBox box;
  RunConfig &run_config = *(data->run_config);

  // 0. Initialize function by allocating memory
  //    and setting labels.

  print ("> Assigning initial conditions (Startup) ...\n");

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);
  #pragma acc enter data copyin(box)
  #pragma omp target enter data map(to:box)

  #ifdef FARGO
  FARGO_Initialize(grid);
  double **wA = FARGO_Velocity();
  #endif

#if ACC_INIT
  #pragma acc enter data copyin(run_config)
  #pragma acc enter data copyin(grid)
  #pragma acc enter data copyin(data->Dts)
  #pragma acc enter data copyin(data)
  // g_useGPU = 1;
  // #pragma acc update device(g_useGPU)

  #pragma omp target enter data map(to:run_config)
  #pragma omp target enter data map(to:grid[:1])
  #pragma omp target enter data map(to:data[:1])
  
  // #pragma omp target update to(g_useGPU)

  #ifdef FARGO
  #pragma acc data present(wA)
  {
  #endif
#endif

  // ------------------------------------------------------
  // 1. Assign initial conditions zone by zone using
  //    primitive or conservative variables
  // ------------------------------------------------------

#if ACC_INIT
  #pragma acc parallel loop collapse(3) present(data, grid, run_config)
  #pragma omp target teams distribute parallel for collapse(3)
#endif
  for (k = 0; k < nx3_tot; k++){
  for (j = 0; j < nx2_tot; j++){
  for (i = 0; i < nx1_tot; i++){

    State state;
    double up[NVAR+5], vp[NVAR+5];

    // 1a. Define coordinates, zone lower and 
    //     upper bounds

    x1  = grid->xC[IDIR][i];
    x2  = grid->xC[JDIR][j];
    x3  = grid->xC[KDIR][k];
    #if (PHYSICS == GRMHD) && BARYCENTER
    x1 = grid->xgc[0](i,j,k);
    x2 = grid->xgc[1](i,j,k);
    x3 = grid->xgc[2](i,j,k);
    #endif

    // 1b. Reset arrays 

    STATE_SET(state, vp, up, NULL, x1, x2, x3);
    NVAR_LOOP(nv) data->Vc(i,j,k,nv) = up[nv] = vp[nv] = 0.0;

    #ifdef GLM_MHD
    vp[PSI_GLM] = up[PSI_GLM] = 0.0;
    #ifdef PHI_GLM
    vp[PHI_GLM] = up[PHI_GLM] = 0.0;
    #endif
    #endif

    // 1c. Simple cell-center assignment used in 
    //     traditional 2nd order schemes.

#ifndef HIGH_ORDER
    Init (vp, x1, x2, x3, run_config);

    #if PHYSICS != GRMHD
    PrimToCons (vp, up, run_config);
    #else
    Metric G = DEFAULT_METRIC;
    if (!SingularMetric(x1,x2,x3)){
      ComputeLocalMetric(x1, x2, x3, &G, run_config);
      ComputeSqrtDet(&G);
      // ComputeGcon(&G);     // useful only for double-check ConsToPrim 
    }
    PrimToCons(vp, up, run_config, &G);
    #endif

    NVAR_LOOP(nv) data->Vc(i,j,k,nv) = vp[nv];
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) = up[nv];
   
    #ifdef FARGO
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    wA[k][i]= vp[FARGO_W];
    #elif GEOMETRY == SPHERICAL
    wA[j][i]= vp[FARGO_W];
    #endif
    #endif

    #if ASSIGN_VECTOR_POTENTIAL == YES
    {
      double x1m = grid->xL[IDIR][i];
      double x2m = grid->xL[JDIR][j];
      double x3m = grid->xL[KDIR][k];
      
      // 1c1. Assign vector potential ? 
      //      Staggered MHD requires A[] to be staggered 
      //      at (Ax1[0++], Ax2[+0+], Ax3[++0]).
      //      All 3 components must be given here.

      Init (vp, x1, x2m, x3m, run_config);
      data->A[IDIR](i,j,k) = vp[AX1];

      Init (vp, x1m, x2, x3m, run_config);
      data->A[JDIR](i,j,k) = vp[AX2];

      Init (vp, x1m, x2m, x3, run_config);
      data->A[KDIR](i,j,k) = vp[AX3];

      // 1c2. Staggered electric field must be assigned directly at 
      //      zone faces.

      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_EXPAND(Init (vp, x1m, x2, x3, run_config);
                 data->Vs[EX1s](i,j,k) = vp[EX1];        ,

                 Init (vp, x1, x2m, x3, run_config);
                 data->Vs[EX2s](i,j,k) = vp[EX2];         ,

                 Init (vp, x1, x2, x3m, run_config);
                 data->Vs[EX3s](i,j,k) = vp[EX3];)
      #endif
    }
    #endif // ASSIGN_VECTOR_POTENTIAL == YES 

    #if (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == NO)
    {
      DIM_EXPAND(double x1m = grid->xL[IDIR][i];,
                 double x2m = grid->xL[JDIR][j];,
                 double x3m = grid->xL[KDIR][k];)
               
      // 1c3. Assign staggered fields components directly from
      //      Init()

      #if INCLUDE_IDIR
      Init (vp, x1m, x2, x3, run_config);
      CT_ASSIGN(data->Vs[BX1s](i,j,k) = vp[BX1];  ,
                data->Vs[EX1s](i,j,k) = vp[EX1];)
      #endif
      #if INCLUDE_JDIR
      Init (vp, x1, x2m, x3, run_config);
      CT_ASSIGN(data->Vs[BX2s](i,j,k) = vp[BX2];  ,
                data->Vs[EX2s](i,j,k) = vp[EX2];)
      #endif
      #if INCLUDE_KDIR
      Init (vp, x1, x2, x3m, run_config);
      CT_ASSIGN(data->Vs[BX3s](i,j,k) = vp[BX3];  ,
                data->Vs[EX3s](i,j,k) = vp[EX3];)
      #endif
    }
    #endif // (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == NO)  

#else // HIGH_ORDER

    double x1p = grid->xR[IDIR][i];
    double x2p = grid->xR[JDIR][j];
    double x3p = grid->xR[KDIR][k];
    double x1m = grid->xL[IDIR][i];
    double x2m = grid->xL[JDIR][j];
    double x3m = grid->xL[KDIR][k];
    
    double u_av[NVAR+5], v_av[NVAR+5];
    #ifdef STAGGERED_MHD
    intList nv_list;
    #endif

    // 1d. Perform Gaussian Integration on cell-centered 
    //     values using Ng Gaussian points.
    //     Convert zone-averaged value to primitive vars.

    #if PHYSICS != GRMHD
    InitGaussAverage (NULL, x1m, x1p, x2m, x2p, x3m, x3p, u_av, 
                      Ng, Ng, Ng, run_config);
    #else
    InitGaussAverage (NULL, x1m, x1p, x2m, x2p, x3m, x3p, u_av, 
                      Ng, Ng, Ng, run_config, grid);
    #endif
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) = u_av[nv];

    STATE_SET(state, v_av, u_av, NULL, x1, x2, x3);
    state.flag = &(data->flag(i,j,k));
    int local_err;
    #if PHYSICS != GRMHD
    local_err = ConsToPrim (state, run_config);
    #else
    local_err = ConsToPrim (state, run_config, &G);
    #endif
    NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v_av[nv];

    // 1e. Repeat Gaussian integration for face-centered values

    #if (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == NO)

    int nstag = 1;  // Number of staggered variable at a given face 
    CT_ASSIGN(nstag = 1;  ,
              nstag = 2);

    #if INCLUDE_IDIR
    nv_list.nvar    = nstag;
    CT_ASSIGN(nv_list.indx[0] = BX1;   ,
              nv_list.indx[1] = EX1;)
    InitGaussAverage (&nv_list, x1m, x1m, x2m, x2p, x3m, x3p, u_av, 
                      Ng, Ng, Ng, run_config);
    CT_ASSIGN(data->Vs[BX1s](i,j,k) = u_av[BX1];  ,
              data->Vs[EX1s](i,j,k) = u_av[EX1];)
    #endif

    #if INCLUDE_JDIR
    nv_list.nvar    = nstag;
    CT_ASSIGN(nv_list.indx[0] = BX2;  ,
              nv_list.indx[1] = EX2;)
    InitGaussAverage (&nv_list, x1m, x1p, x2m, x2m, x3m, x3p, u_av, 
                      Ng, Ng, Ng, run_config);
    CT_ASSIGN(data->Vs[BX2s](i,j,k) = u_av[BX2];,
              data->Vs[EX2s](i,j,k) = u_av[EX2];)
    #endif

    #if INCLUDE_KDIR
    nv_list.nvar    = nstag;
    CT_ASSIGN(nv_list.indx[0] = BX3;  ,
              nv_list.indx[1] = EX3;)
    InitGaussAverage (&nv_list, x1m, x1p, x2m, x2p, x3m, x3m, u_av, 
                      Ng, Ng, Ng, run_config);
    CT_ASSIGN(data->Vs[BX3s](i,j,k) = u_av[BX3];,
              data->Vs[EX3s](i,j,k) = u_av[EX3];)
    #endif // INCLUDE_KDIR 
    #endif // (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == NO)  

    // 1f. Repeat for line-averaged vector potential
    //     (ASSIGN_VECTOR_POTENTIAL == YES)

    #if (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == YES)

    #if INCLUDE_IDIR && INCLUDE_JDIR
    nv_list.nvar = 1;
    nv_list.indx[0] = AX3;
    InitGaussAverage (&nv_list, x1m, x1m, x2m, x2m, x3m, x3p, u_av, 
                      Ng, Ng, Ng, run_config);
    data->A[KDIR](i,j,k) = u_av[AX3];
    #endif

    #if INCLUDE_IDIR && INCLUDE_KDIR
    nv_list.nvar = 1;
    nv_list.indx[0] = AX2;
    InitGaussAverage (&nv_list, x1m, x1m, x2m, x2p, x3m, x3m, u_av, 
                      Ng, Ng, Ng, run_config);
    data->A[JDIR](i,j,k) = u_av[AX2];
    #endif

    #if INCLUDE_JDIR && INCLUDE_KDIR
    nv_list.nvar = 1;
    nv_list.indx[0] = AX1;
    InitGaussAverage (&nv_list, x1m, x1p, x2m, x2m, x3m, x3m, u_av, 
                      Ng, Ng, Ng, run_config);
    data->A[IDIR](i,j,k) = u_av[AX1];
    #endif
    #endif // (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == YES) 

#endif // HIGH_ORDER 

  }}}
#if !ACC_INIT
  #pragma acc enter data copyin(run_config)
  #pragma acc enter data copyin(grid)
  #pragma acc enter data copyin(data->Dts)
  #pragma acc enter data copyin(data)
  #pragma omp target enter data map(to:run_config[:1])
  #pragma omp target enter data map(to:grid[:1])
  #pragma omp target enter data map(to:data[:1])
#else
  #pragma acc update host(data->run_config)
  #pragma omp target update from(data->run_config)
  
  #ifdef FARGO
  #pragma acc update host(wA)
  }
  #endif
#endif

  // 2. Assign initial conditions by calling InitDomain()

  InitDomain(data, grid);

  // 3. Compute staggered fields from vector potential if required.
  //    The vector potential data->A[]() is always stored at zone edges, 
  //    regardless of the div-free method:
  // 
  //    - For CT, we directly compute the curl of A into data->Vs;
  //    - For cell-centered MHD, we store the face magnetic field
  //      into data->J[] and then average field components to the center.

  #if (ASSIGN_VECTOR_POTENTIAL == YES) && (defined STAGGERED_MHD)

  #if PHYSICS != GRMHD
  VectorCurl(&(data->A[0]), &(data->Vs[0]), grid, &box);  // B = Curl(A)
  #else
  VectorCurl(&(data->A[0]), &(data->Vs[0]), grid, &box, run_config);  // B = Curl(A)
  #endif

  #elif (ASSIGN_VECTOR_POTENTIAL == YES) && !(defined STAGGERED_MHD)

  data->J[IDIR].location_ = X1_FACE;
  data->J[JDIR].location_ = X2_FACE;
  data->J[KDIR].location_ = X3_FACE;
  #if PHYSICS != GRMHD
  VectorCurl(&(data->A[0]), &(data->J[0]), grid, &box);  // J = B = Curl(A)
  #else
  VectorCurl(&(data->A[0]), &(data->J[0]), grid, &box, run_config);  // J = B = Curl(A)
  #endif
  FieldAverage(data->J, data->Vc, grid);

  #endif //  (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == YES)

  // 4. Compute cell-centered magnetic field by averaging face values.
  //    Useful only for saving the first output, since the
  //    average will be re-computed anyway at the beginning
  //    of the computation.

  #ifndef HIGH_ORDER
  #ifdef STAGGERED_MHD
  RBox box_tot;
  RBoxDefine (0, nx1_tot-1,
              0, nx2_tot-1,
              0, nx3_tot-1, &box_tot);

  CT_AverageStaggeredFields (data, 0, &box_tot, grid);

  #pragma acc parallel loop collapse(3) present(data, grid) copyin(box_tot)
  #pragma omp target teams distribute parallel for collapse(3)
  BOX_LOOP(box_tot,i,j,k){
    double u[NVAR];

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet (metric, grid->rC[i], grid->sC[j]);
    #endif

    x1  = grid->xC[IDIR][i];
    x2  = grid->xC[JDIR][j];
    x3  = grid->xC[KDIR][k];

    NVAR_LOOP(nv) u[nv] = data->Uc(i,j,k,nv);
    #if GEOMETRY != CARTESIAN
    ConsToOrtho(u, metric);
    #endif

    #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    DIM_EXPAND(
      data->Vc(i,j,k,BX1) = u[BX1]; ,
      data->Vc(i,j,k,BX2) = u[BX2]; ,
      data->Vc(i,j,k,BX3) = u[BX3];
    )
    #endif
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(
      data->Vc(i,j,k,EX1) = u[EX1]; ,
      data->Vc(i,j,k,EX2) = u[EX2]; ,
      data->Vc(i,j,k,EX3) = u[EX3];
    )
    #endif

  }
  #endif // STAGGERED_MHD   
  #endif // #ifndef HIGH_ORDER 

  // 5. A few more operations:
  //    - compute charge for ResRMHD

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == DIV_CLEANING)
  {
    #pragma acc parallel loop collapse(3) present(data, grid, box)
    #pragma omp target teams distribute parallel for collapse(3)
    BOX_LOOP(box,i,j,k){
      DIM_EXPAND(double dx1 = grid->dx[IDIR][i];,
                 double dx2 = grid->dx[JDIR][j];,
                 double dx3 = grid->dx[KDIR][k]; )

      data->Vc(i,j,k,CRG) = DIM_EXPAND(  0.5*(data->Vc(i+1,j,k,EX1) - data->Vc(i-1,j,k,EX1))/dx1,
                                       + 0.5*(data->Vc(i,j+1,k,EX2) - data->Vc(i,j-1,k,EX2))/dx2,
                                       + 0.5*(data->Vc(i,j,k+1,EX3) - data->Vc(i,j,k-1,EX3))/dx3 );
    }
  }
  #endif  // (PHYSICS == ResRMHD) && (DIVE_CONTROL == DIV_CLEANING)

  // 6. Check physical admissibility

  uint32_t err_kernel = 0;  // Cumulative kernel error flag  

  #pragma acc parallel loop collapse(3) present(data, grid) copy(err_kernel) 
  #pragma omp target teams distribute parallel for collapse(3)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    double vp[NVAR];
    x1 = grid->xC[IDIR][i];
    x2 = grid->xC[JDIR][j];
    x3 = grid->xC[KDIR][k];

    NVAR_LOOP(nv) vp[nv] = data->Vc(i,j,k,nv);

    if (vp[RHO] <= 0.0) {
      ATOMIC_ERR_UPDATE(err_kernel,    ERR_NEGATIVE_RHO, ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_NEGATIVE_RHO, ERR_QUIT);
    }
    #if HAVE_ENERGY
    if (vp[PRS] <= 0.0) {
      ATOMIC_ERR_UPDATE(err_kernel,    ERR_NEGATIVE_PRS, ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_NEGATIVE_PRS, ERR_QUIT);
    }
    #endif

    #if (PHYSICS == RHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
    double scrh = vp[VX1]*vp[VX1] + vp[VX2]*vp[VX2] + vp[VX3]*vp[VX3];
    #if PHYSICS == GRMHD
    Metric G = DEFAULT_METRIC;
    ComputeLocalMetric(x1, x2, x3, &G, run_config);
    scrh = DotProduct(vp+VX1,vp+VX1,G.gCov);
    #endif
    if (scrh >= 1.0){
      ATOMIC_ERR_UPDATE(err_kernel,    ERR_SUPERLUMINAL_SPEED, ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_SUPERLUMINAL_SPEED, ERR_QUIT);
    }
    #endif

    #if RADIATION
    double f = sqrt(vp[FR1]*vp[FR1] + vp[FR2]*vp[FR2] + vp[FR3]*vp[FR3])/vp[ENR];
    if (f >= 1.0){
      ATOMIC_ERR_UPDATE(err_kernel,    ERR_SUPERLUMINAL_RADIATION_FLUX, ERR_QUIT);
      KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_SUPERLUMINAL_RADIATION_FLUX, ERR_QUIT);
    }
    #endif
  }}}

  Where ("Startup(): unphysical state", data, err_kernel, &box, grid);

  int quit_condition = (err_kernel & ERR_QUIT);
  if (quit_condition) QuitPLUTO("Startup()",err_kernel);

  // 7. Compute auxiliary units

  {
    const double rho0 = run_config.unit_density;
    const double vel0 = run_config.unit_velocity;
    const double len0 = run_config.unit_length;
    run_config.unit_time = len0/vel0;
    run_config.unit_Bmag = vel0*sqrt(4.0*CONST_PI*rho0);
    #pragma acc update device(run_config.unit_time, run_config.unit_Bmag)
    #pragma omp target update to(run_config.unit_time, run_config.unit_Bmag)
  }

  // 8. Convert to conservative

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost)
  #if (HO_LAP_LIMITER != NO) && (GEOMETRY == CARTESIAN)  
  DiscontinuityDetect(data, grid, &box);   
  #endif
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc
  #endif

  #ifndef HIGH_ORDER
  // RBoxDefine (grid->lbeg, grid->lend, &box);
  #pragma acc parallel loop collapse(3) present(data, grid, run_config, box)
  #pragma omp target teams distribute parallel for collapse(3)
  BOX_LOOP(box,i,j,k){
    double vc[NVAR], up[NVAR];
    x1 = grid->xC[IDIR][i]; 
    x2 = grid->xC[JDIR][j]; 
    x3 = grid->xC[KDIR][k]; 

    #if (GEOMETRY != CARTESIAN) && (PHYSICS != GRMHD) 
    Metric metric;
    MetricSet (metric, grid->rC[i], grid->sC[j]);    
    #elif PHYSICS == GRMHD
    Metric G = DEFAULT_METRIC;
    #endif

    NVAR_LOOP(nv) vc[nv] = data->Vc(i,j,k,nv);
    #if PHYSICS != GRMHD
    PrimToCons(vc, up, run_config);
    #if GEOMETRY != CARTESIAN
    OrthoToCons (up, metric);
    #endif
    #else
    if (!SingularMetric(x1,x2,x3)){
      ComputeLocalMetric(x1, x2, x3, &G, run_config);
      ComputeSqrtDet(&G);
      // ComputeGcon(&G);     // useful only for double-check ConsToPrim 
    }
    PrimToCons(vc, up, run_config, &G);
    #endif
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) = up[nv];

  }
  #endif  // !HIGH_ORDER

  #pragma omp target exit data map(delete:box)
  #pragma acc exit data delete(box) 
}

#if (ASSIGN_VECTOR_POTENTIAL == YES) && !(defined STAGGERED_MHD)
////////////////////////////////////////////////////////////////////////
/*!
  For central-type discretizaions (no CONSTRAINED_TRANSPORT),
  computes the volume-averaged magnetic fields from the face values.
  
  \return No return value.
*/
////////////////////////////////////////////////////////////////////////
void FieldAverage(Ary3D *Bf, Ary4D &Vc, Grid *grid)
{
  int i,j,k;

  int nx1_tot=grid->ntot[IDIR];
  int nx2_tot=grid->ntot[JDIR];
  int nx3_tot=grid->ntot[KDIR];

  DIM_EXPAND(Ary3D &Bx1 = Bf[IDIR];,
             Ary3D &Bx2 = Bf[JDIR];,
             Ary3D &Bx3 = Bf[KDIR]; )

  #pragma acc parallel loop collapse(3) present(Bf, Vc)
  #pragma omp target teams distribute parallel for collapse(3)
  for (k = 2*INCLUDE_KDIR; k < nx3_tot-2*INCLUDE_KDIR; k++){
  for (j = 2*INCLUDE_JDIR; j < nx2_tot-2*INCLUDE_JDIR; j++){
  for (i = 2*INCLUDE_IDIR; i < nx1_tot-2*INCLUDE_IDIR; i++){
    double bp[3], bm[3];
    double bavg[3];

    DIM_EXPAND(bm[IDIR] = Bx1(i,j,k); bp[IDIR] = Bx1(i+1,j,k);,
               bm[JDIR] = Bx2(i,j,k); bp[JDIR] = Bx2(i,j+1,k);,
               bm[KDIR] = Bx3(i,j,k); bp[KDIR] = Bx3(i,j,k+1); )

    // 1b. Form averages 

    #if   GEOMETRY == CARTESIAN
    DIM_EXPAND(bavg[IDIR] = 0.5*(bp[IDIR] + bm[IDIR]);,
               bavg[JDIR] = 0.5*(bp[JDIR] + bm[JDIR]);,
               bavg[KDIR] = 0.5*(bp[KDIR] + bm[KDIR]); )
    #elif GEOMETRY == CYLINDRICAL 
    double rp = grid->xR[IDIR][i];
    double rm = grid->xL[IDIR][i];
    DIM_EXPAND(bavg[IDIR] = (rp*bp[IDIR] + rm*bm[IDIR])/(rp + rm);,
               bavg[JDIR] = 0.5*(bp[JDIR] + bm[JDIR]);            ,
               bavg[KDIR] = 0.5*(bp[KDIR] + bm[KDIR]);             )

    #elif GEOMETRY == SPHERICAL
    double rp   = grid->xR[IDIR][i];
    double rm   = grid->xL[IDIR][i];
    double dx1  = grid->dx[IDIR][i];
    double dV1  = (rp*rp*rp - rm*rm*rm)/3.;
    double thp  = grid->xR[JDIR][j];
    double thm  = grid->xL[JDIR][j];
    double sp   = fabs(sin(thp));
    double sm   = fabs(sin(thm));
    double dx2  = grid->dx[JDIR][j];
    double dV2  = fabs(cos(thm) - cos(thp));

    DIM_EXPAND(bavg[IDIR] = 0.5*(rp*rp*bp[IDIR] + rm*rm*bm[IDIR])*dx1/dV1;  ,
               bavg[JDIR] = 0.5*(sp*bp[JDIR] + sm*bm[JDIR])*dx2/dV2;        ,
               bavg[KDIR] = 0.5*(bp[KDIR] + bm[KDIR]);                       )
    #endif

    DIM_EXPAND(Vc(i,j,k,BX1) = bavg[IDIR];,
               Vc(i,j,k,BX2) = bavg[JDIR];,
               Vc(i,j,k,BX3) = bavg[KDIR]; )
  }}}
}
#endif
