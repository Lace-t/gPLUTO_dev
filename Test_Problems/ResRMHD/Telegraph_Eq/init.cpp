////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Telegraph Equation

*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "rotate.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*!
  Assign initial conditions
*/
////////////////////////////////////////////////////////////////////////
{
  int jrot_x, jrot_y, krot_x, krot_z;
  // static int first_call = 1;
  double *inputParam = run_config.inputParam;
  double xvec[3] = {x,y,z};
  double Lx  = run_config.domEnd[IDIR] - run_config.domBeg[IDIR];
  double Ly  = run_config.domEnd[JDIR] - run_config.domBeg[JDIR];
  double Lz  = run_config.domEnd[KDIR] - run_config.domBeg[KDIR];
  DIM_EXPAND(double kx  = 2.0*CONST_PI/Lx*(DIMENSIONS >= 1);,
             double ky  = 2.0*CONST_PI/Ly*(DIMENSIONS >= 2);,
             double kz  = 2.0*CONST_PI/Lz*(DIMENSIONS >= 3); )
  double phi = DIM_EXPAND(kx*x, + ky*y, + kz*z);
  double theta_yz = inputParam[THETA_YZ]*CONST_PI/180.0; 

  // ----------------------------------------------
  // 0. Round cos(theta_yz) to 1 or 0 when THETA_YZ
  // (angle between magnetic field and y-axis is 
  // close to 0 or 90
  // ----------------------------------------------
   
  double ctheta = cos(theta_yz);
  double stheta = sin(theta_yz);
  
  if (fabs(theta_yz) < 1.e-5) {
    ctheta = 1.0; stheta = 0.0;
  }else if ( fabs(theta_yz - 0.5*CONST_PI) < 1.e-5) {
    ctheta = 0.0; stheta = 1.0;
  }else if ( fabs(theta_yz - CONST_PI) < 1.e-5) {
    ctheta = -1.0; stheta = 0.0;
  }

  // ----------------------------------------------
  // 1. Define rotation angles in terms of
  // integer shifts. Compute coordinates in
  // non-rotated frame (x1, y1, z1).
  // ----------------------------------------------

  jrot_y = 1;
  jrot_x = round(-jrot_y/Ly);

  krot_z = 1;
  krot_x = round(-krot_z/Lz);

  if (inputParam[APPLY_ROTATION] < 0.5) {
    jrot_x = krot_x = 0.0;
    DIM_EXPAND(,ky = 0.0;,kz = 0.0);
  }  

  RotateSet(jrot_x, jrot_y, krot_x, krot_z, run_config);
  RotateVector(xvec,1);

  // ----------------------------------------------
  // 2. Assign initial conditions
  // ----------------------------------------------

  v[RHO] = 1.e12;

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  double k     = sqrt(DIM_EXPAND(kx*kx, + ky*ky, + kz*kz));
  double eta   = inputParam[ETA];
  double sigma = 1.0/eta;
  double mu    = sqrt(k*k - 0.25*sigma*sigma);
  double Bz = cos(phi);

  v[BX1] = 0.0;
  v[BX2] = Bz*ctheta;
  v[BX3] = Bz*stheta;

  v[PRS] = 1.0;

  double Ey = 0.5*(2.0*mu*cos(phi) + sigma*sin(phi))/k;

  v[EX1] = 0.0;
  v[EX2] =  Ey*stheta;
  v[EX3] = -Ey*ctheta;
  
#if PHYSICS == ResRMHD
  #ifdef CRG
  v[CRG] = 0.0;
  #endif

  #ifdef PSI_GLM
  v[PSI_GLM] = 0.0;
  #endif
  #ifdef PHI_GLM
  v[PHI_GLM] = 0.0;
  #endif
#endif

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0; // 1.0/(k*k)*(mu*sin(phi) - 0.5*sigma*cos(phi));

  // ----------------------------------------------
  // 3. Rotate vectors
  // ----------------------------------------------

  RotateVector (v + AX1, -1);
  RotateVector (v + VX1, -1);
  RotateVector (v + BX1, -1);

  #if PHYSICS == ResRMHD
  RotateVector (v + EX1, -1);
  #endif
}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*! 
  [Not Used]
*/
////////////////////////////////////////////////////////////////////////
{ 

}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*! 
  
*/
////////////////////////////////////////////////////////////////////////
{
  int    i,j,k,nv;
  int    nx1_glob = grid->nint_glob[IDIR];
  int    nx2_glob = grid->nint_glob[JDIR];
  int    nx3_glob = grid->nint_glob[KDIR];
  int     ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int     jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int     kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  double *inputParam = data->run_config->inputParam;
  double theta_yz = inputParam[THETA_YZ]*CONST_PI/180.0;
  double sigma    = 1.0/inputParam[ETA];
  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];
  double damp, errEx, errEy, errBz;
  double v[256], vex[256];
  double Bz, Ey;
  double Bzex, Eyex;
  double time = data->Dts->time;
  RBox box;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box);

  // -------------------------------------------------------- 
  // 1. Compute error in the original 1D frame.
  // (We align the vectors so that B lies along the 
  // z-axis and E along the y-axis)
  // --------------------------------------------------------

  errEx = errEy = errBz = 0.0;
  damp = exp(-0.5*sigma*time);

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);
  PointValue(data, grid, &box);
  #endif

  BOX_LOOP(box,i,j,k){

    #ifdef HIGH_ORDER
    NVAR_LOOP(nv) v[nv] = data->Vpc(i,j,k,nv);
    #else
    NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k,nv);
    #endif

    // -------------------------------------------- 
    // 1a. Compute exact solution in 2D and rotate 
    // vectors clockwise around the z-axis 
    // --------------------------------------------
    
    Init (vex, x1[i], x2[j], x3[k], *(data->run_config));  // Exact solution 
    
    RotateVector (vex + BX1, 1);
    RotateVector (vex + EX1, 1);

    vex[BX1] *= damp;
    vex[BX2] *= damp;
    vex[BX3] *= damp;

    vex[EX1] *= damp;
    vex[EX2] *= damp;
    vex[EX3] *= damp;

    // -------------------------------------------- 
    // 1b. Further rotate exact solution vectors 
    // around the x-axis so that B lies along 
    // the z-direction and E along y.
    // --------------------------------------------

// theta_yz = M_PI*0.5;

    Eyex = vex[EX2]*sin(theta_yz) - vex[EX3]*cos(theta_yz);
    // Ezex = vex[EX2]*cos(theta_yz) + vex[EX3]*sin(theta_yz);

    // Byex = vex[BX2]*sin(theta_yz) - vex[BX3]*cos(theta_yz);
    Bzex = vex[BX2]*cos(theta_yz) + vex[BX3]*sin(theta_yz);

    // -------------------------------------------- 
    // 1c. Rotate vectors clockwise around z-axis
    // --------------------------------------------

    RotateVector (v + BX1, 1);
    RotateVector (v + EX1, 1);

    // -------------------------------------------- 
    // 1d. Forther rotate vectors around the 
    // x-axis in the yz plane
    // --------------------------------------------

    Ey = v[EX2]*sin(theta_yz) - v[EX3]*cos(theta_yz);
    // Ez = v[EX2]*cos(theta_yz) + v[EX3]*sin(theta_yz);

    // By = v[BX2]*sin(theta_yz) - v[BX3]*cos(theta_yz);
    Bz = v[BX2]*cos(theta_yz) + v[BX3]*sin(theta_yz);



// if (Ez > 1e-9 & By > 1e-9) {
//   printf ("[Ey, Ez]_num = %12.6e, %12.6e; [By, Bz]_num = %12.6e, %12.6e\n", 
//          Ey, Ez, By, Bz);
// }
// printf ("[Ey, Ez]_ref = %12.6e, %12.6e; [By, Bz]_ref = %12.6e, %12.6e\n",
//          Eyex, Ezex, Byex, Bzex);

    // -------------------------------------------- 
    // 1e. Compute error
    // --------------------------------------------

    errEx += fabs(v[EX1]);
    errEy += fabs(Eyex - Ey);
    errBz += fabs(Bzex - Bz);
  }

  // -------------------------------------------------------- 
  // 2. Reduce data among processors
  // --------------------------------------------------------

  #ifdef PARALLEL
  double err_glob;
  MPI_Allreduce (&errEx, &err_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  errEx = err_glob;
  MPI_Allreduce (&errEy, &err_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  errEy = err_glob;
  MPI_Allreduce (&errBz, &err_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  errBz = err_glob;
  #endif

  errEx = errEx/(nx1_glob*nx2_glob*nx3_glob);
  errEy = errEy/(nx1_glob*nx2_glob*nx3_glob);
  errBz = errBz/(nx1_glob*nx2_glob*nx3_glob);

  // -------------------------------------------------------- 
  // 3. Write to disk
  // --------------------------------------------------------
// printf("ANALYSIS\n");
  if (g_rank == 0){
    FILE *fp;
    if (g_stepNumber == 0){
      fp = fopen ("telegraph_eq.dat","w");
    }else{
      fp = fopen ("telegraph_eq.dat","a");
    }

    fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e\n", 
                 time, errEx, errEy, errBz);
    fclose(fp);
  }

}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
/*! 
  [Not Used]
*/
////////////////////////////////////////////////////////////////////////
{ 
}



