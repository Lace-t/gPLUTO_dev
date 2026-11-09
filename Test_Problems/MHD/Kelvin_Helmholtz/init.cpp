////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Kelvin-Helmholtz instability setup.

  We consider a simple 2D, shear layer configuration, with constant
  density and pressure.
  The velocity profile is simply given by \f[ v_y = M/2*\tanh(y/w), \f] where
  M is the relative Mach number, w is the shear width.
  Velocities are normalized to the sound speed, density to the fluid
  initial density and length to the total domain size.
  A passive tracer is used to distinguish the two fluids (+1 and -1 for the
  fluid above and below, respectively).
  
  For a magnetized plasma, two additional parameters, namely the Alfven
  velocity (vA) and the angle theta between the magnetic field and
  the x-direction are used to control the field strength and orientation.

  Input parameters are:
  
  - inputParam[AMPL] = perturbation amplitude
  - inputParam[MACH] = relative velocity (and Mach number)
                       between the two layers;
  - inputParam[WIDTH] = width of the shear;

  For MHD, 
 
  - inputParam[VALF]  = Alfven velocity in units of the sound speed (=1)
  - inputParam[THETA] = angle (in degrees) between the magnetic field 
                        and the x direction.
  
  Below the pre-defined configurations:

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM|   REC     |TIME_STEP  |DIVB| Solver
  -----|----|-----------|---|-----------|-----------|----|-------
   #01 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | roe
   #02 |MHD |CARTESIAN  | 2 | LimO3     | RK3       | CT | hlld
  </CENTER>

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Jul 31, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/* 
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
*/
////////////////////////////////////////////////////////////////////////
{
  int  i, j, k, nv;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  RunConfig &run_config = *(data->run_config);

  double *inputParam = run_config.inputParam;
  double B0, rnd, scrh;
  double x,y,z;
  double *xC   = grid->xC[IDIR], *yC = grid->xC[JDIR];
  double *xL   = grid->xL[IDIR], *yL = grid->xL[JDIR];
  double Ly    = grid->xend_glob[JDIR] - grid->xbeg_glob[JDIR];
  double M     = inputParam[MACH];
  double eps   = inputParam[AMPL]*M;  // Perturbation amplitude 
  double a     = inputParam[WIDTH];
  double vA    = inputParam[VALF];
  double theta = inputParam[THETA]/180.0*CONST_PI;
  double v[256];
  #ifdef HIGH_ORDER
  double u[256];
  #endif

  RBox   box, box_tot;

  // 0. Define boxes

  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine (0, grid->ntot[IDIR]-1, 
              0, grid->ntot[JDIR]-1, 
              0, grid->ntot[KDIR]-1,
              &box_tot);  
 
  // --------------------------------------------
  // 1. Assign zone-centered variables
  // --------------------------------------------

   BOX_LOOP(box_tot,i,j,k){

    // 1a. Global grid index (start at ibeg)

    int64_t ix = i - 2*grid->lbeg[IDIR] + grid->beg[IDIR];
    int64_t iy = j - 2*grid->lbeg[JDIR] + grid->beg[JDIR];
    int64_t iz = k - 2*grid->lbeg[KDIR] + grid->beg[KDIR];

    // Should be the same as
    // ix = floor((x - run_config.domBeg[IDIR])/dx);
    // iy = floor((y - run_config.domBeg[JDIR])/dy);
    // iz = floor((z - run_config.domBeg[KDIR])/dz);

    y = yC[j];

    // 1b. Avoid using the fixed-seed technique 
    //     to spawn always the same sequence of 
    //     pseudo-random numbers. Safer to use a 
    //     unique generator using triplets of 
    //     (x,y,z) and a hashing function. 

    rnd = GridRandomNumber(fabs(ix), fabs(iy), fabs(iz));  // Generate prn in [0,1]
    rnd = (2.0*rnd-1); // Map into [-1,1]

    scrh = y/(0.25*Ly);

    v[RHO] = 1.0;                      // Set constant density  
    v[VX1] = 0.5*M*tanh(y/a);          // Main flow             
    v[VX2] = eps*rnd*exp(-scrh*scrh);  // Vertical perturbation 
    v[VX3] = 0.0;

    #if EOS == IDEAL
    v[PRS] = 1.0/run_config.gamma;     // Set constant pressure (cs = 1) 
    #elif EOS == ISOTHERMAL
    g_isoSoundSpeed = 1.0;
    #endif
    v[TRC]  = (y < 0.0 ? 1.0:-1.0);    // Passive tracer (color)        

    #if PHYSICS == MHD
    B0     = vA*sqrt(v[RHO]);

    v[BX1] = B0*cos(theta);    // Magnetic field lies in the x-z plane.        
    v[BX2] = 0.0;              // Theta is the angle between B and x-direction
    v[BX3] = B0*sin(theta); 

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = y*v[BX1];
    #endif

    // 1c. Copy back to main array. 
    //     For high order, compute Upc

    NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];

    #ifdef HIGH_ORDER
    PrimToCons(v, u, run_config);
    NVAR_LOOP(nv) data->Upc(i,j,k,nv) = u[nv];
    #endif
  }

  // --------------------------------------------------------   
  // 2. Assign staggered fieds
  // --------------------------------------------------------

  #ifdef STAGGERED_MHD
  RBoxDefine (ibeg-INCLUDE_IDIR, iend+INCLUDE_IDIR, 
              jbeg-INCLUDE_JDIR, jend+INCLUDE_JDIR, 
              kbeg-INCLUDE_KDIR, kend+INCLUDE_KDIR, &box);

  BOX_LOOP(box, i,j,k){
    x = xL[i]; y = yC[j];                // Change coordinates to x1-face
    data->Vs[BX1s](i,j,k) = B0*cos(theta); 
  }

  BOX_LOOP(box, i,j,k){
    x = xC[i]; y = yL[j];                // Change coordinates to x2-face
    data->Vs[BX2s](i,j,k) = B0*sin(theta);  
  }

  #if ASSIGN_VECTOR_POTENTIAL
  BOX_LOOP(box,i,j,k){
    x = xL[i]; y = yL[j];                // Change coordinates to x3-edge 
    d->Ax1(i,j,k) = 0.0;
    d->Ax2(i,j,k) = 0.0;
    d->Ax3(i,j,k) = y*B0*cos(theta);
  }
  #endif

  #endif

  // 3a. Assign zone-centered volume-averages (high order only)

  #ifdef HIGH_ORDER
  RBoxDefine (ibeg, iend, 
              jbeg, jend, 
              kbeg, kend, &box);
  State state;
  BOX_LOOP(box, i,j,k){
    double d2ux = 0.0, d2uy = 0.0, d2uz = 0.0;
    double d2u;
    NVAR_LOOP(nv) {
      DIM_EXPAND(d2ux = HO_Average4D<IDIR>(data->Upc,i,j,k,nv);  ,
                 d2uy = HO_Average4D<JDIR>(data->Upc,i,j,k,nv);  ,
                 d2uz = HO_Average4D<KDIR>(data->Upc,i,j,k,nv);)      
      d2u = d2ux + d2uy + d2uz;
      data->Uc(i,j,k,nv) = u[nv] = data->Upc(i,j,k,nv) + d2u/24.0;
    }

    STATE_SET(state, v, u, NULL, x, y, z);
    state.flag = &(data->flag(i,j,k));
    int local_err = ConsToPrim (state, run_config);
    NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];
  }

  // 3b. Assign line-averaged volume-averages (high order only)

  // No need to assign line-averaged Ax3 since we're in 2D. 
  #endif
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*
  Runtime analysis.
*/
////////////////////////////////////////////////////////////////////////
{
  int    i,j,k,nv;  
  double max_vy  = -10.0;
  double min_vy  =  10.0;
  double max_By2 = 0.0;
  double pm = 0.0;
  double scrh, v[NVAR];
  double *dx = grid->dx[IDIR];
  double *dy = grid->dx[JDIR];
  #ifdef PARALLEL
  double vol =  (grid->xend_glob[IDIR] - grid->xbeg_glob[IDIR])
               *(grid->xend_glob[JDIR] - grid->xbeg_glob[JDIR]);
  #endif
  RBox box;

  RBoxDefine (grid->lbeg, grid->lend, &box);

  BOX_LOOP(box, i,j,k){
 
    NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k,nv);

    max_vy = MAX(max_vy, v[VX2]);
    min_vy = MIN(min_vy, v[VX2]);
    
    #if PHYSICS == MHD
    scrh    = v[BX2]*v[BX2];
    max_By2 = MAX(max_By2, scrh);

    pm += 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2])*dx[i]*dy[j];
    #endif
  }

// -- parallel reduce -- 

  #ifdef PARALLEL
  MPI_Allreduce (&max_vy, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  max_vy = scrh;

  MPI_Allreduce (&min_vy, &scrh, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  min_vy = scrh;

  MPI_Allreduce (&max_By2, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  max_By2 = scrh;

  MPI_Allreduce (&pm, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  pm = scrh/vol;

  MPI_Barrier (MPI_COMM_WORLD);
  #endif

// --  Proc #0 writes to disk -- 

  if (g_rank == 0){
    FILE *fp;
    const char *legend[4] = {"t","max-min","max(By^2)","pm"};
    if (g_stepNumber == 0) {
      fp = fopen("khi.dat","w");
      PrintColumnLegend (legend,4,fp);
    } else fp = fopen("khi.dat","a");

    fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e\n",data->Dts->time, 
                                             0.5*(max_vy - min_vy),
                                             max_By2, pm);
    fclose(fp);
  }

}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
/*
  
*/
////////////////////////////////////////////////////////////////////////
{ 
}

