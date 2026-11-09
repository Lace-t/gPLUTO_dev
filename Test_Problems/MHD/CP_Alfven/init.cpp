////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD Circularly polarized Alfven waves.

  This is a rotated version of a 1D setup where 

  \f[
     \rho = 1                      \,,\quad
     v_x  = V_0                    \,,\quad
     v_y  = |\epsilon|\sin(\phi)   \,,\quad
     v_z  = |\epsilon|\cos(\phi)   \,,\quad
     B_x  = 1                      \,,\quad 
     B_y  = \epsilon\sqrt{\rho}\sin(\phi)     \,,\quad
     B_z  = \epsilon\sqrt{\rho}\cos(\phi)     \,,\quad
     P    = P_0                    \,,
  \f]
 
  where \f$V_0\f$ is the translation velocity in the \f$x\f$ direction,
  \f$\phi\f$ is the phase (\f$k_xx\f$ in 1D and
  (\f$k_xx + k_yy\f$) in 2D), and \f$\epsilon\f$ is the wave amplitude
  (\f$\epsilon > 0\f$ implies right going waves;
  \f$\epsilon < 0\f$ implies left going waves).

  With this normalization, the Alfven velocity \f$V_A = B_x/\sqrt{\rho}\f$
  is always unity.

  The configuration is rotated by specifying \f$\tan\alpha = k_y/k_x\f$ and
  \f$\tan\beta = k_z/k_x\f$ which express the ratios between the \f$y\f$- and
  \f$z\f$- components of the wave vector with the \f$x\f$ component.
  In order to apply periodic boundary conditions everywhere,
  an integer number of wavelegnths must be contained in each direction,
  that is, \f$k_x = 2\pi/L_x\f$, \f$k_y = 2\pi/L_y\f$, \f$k_z = 2\pi/L_z\f$.
  
  We use the tools in *rotate.c* which requires to specify the four
  integer shifts such that
  \f[
      \tan\alpha = -\frac{j_x}{j_y}\frac{\Delta x}{\Delta y} = \frac{L_x}{L_y} ,
      \qquad
      \tan\beta  = -\frac{k_x}{k_z}\frac{\Delta x}{\Delta z} = \frac{L_x}{L_z}
  \f]

  The final time step is one period and is found from 
                                
  \f$V_A = \omega/|k|\f$ --> \f$1/T = \sqrt{1 + \tan^2\alpha + \tan^2\beta}\f$

  The runtime parameters that are read from \c pluto.ini are 
  - <tt>run_config.inputParam[EPS]</tt>:       sets the wave amplitude \f$\epsilon\f$;
  - <tt>run_config.inputParam[VEL0]</tt>:      sets \f$V_0\f$;
  - <tt>run_config.inputParam[PR0]</tt>:       sets the pressure of the 1D solution;
  - <tt>run_config.inputParam[ALPHA_GLM]</tt>: ; 

  Configurations:

  <CENTER>
  Conf.|PHYS| GEOMETRY  |DIM|   REC     |TIME_STEP  |DIVB| Solver 
  -----|----|-----------|---|-----------|-----------|----|--------
   #01 |MHD |CARTESIAN  | 2 | LINEAR    | RK3       | CT | roe  
   #02 |MHD |CARTESIAN  | 2 | LINEAR    | RK3       | DC | roe    
   #03 |MHD |CARTESIAN  | 2 | PARABOLIC | RK3       | CT | hllc   
   #04 |MHD |CARTESIAN  | 3 | WENOZ     | RK3       | DC | hlld    
   #05 |MHD |CARTESIAN  | 2 | WENOZ     | RK3       | DC | hll 
   #06 |MHD |CARTESIAN  | 2 | LINEAR    | RK2       | CT | hlld    
   #07 |MHD |CARTESIAN  | 2 | MP5       | RK4 [HO]  | CT | hlld   
   #08 |MHD |CARTESIAN  | 3 | WENOZ     | RK4 [HO]  | CT | hlld   
  </CENTER>

  \author A. Mignone (andrea.mignone@unito.it)\n
          V. Berta   (vittoria.berta@unito.it)

  \date   June 11, 2019

  \b References: 
     - "High-order conservative finite difference GLM-MHD 
        schemes for cell-centered MHD"
        Mignone, Tzeferacos & Bodo JCP (2010)
     - "A 4th-order accurate finite volume method
        for ideal classical and special relativistic MHD 
        based on pointwise reconstructions"
        Berta, Mignone, Bugli & Mattia JCP (2024)
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "rotate.hpp"

////////////////////////////////////////////////////////////////////////
/*!
   Assign initial conditions in terms of primitive variables
*/
////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
{
  int jshift_x, jshift_y, kshift_x, kshift_z;
  double eps;
  double kx, ky,kz;
  double phi;
  double Ly = run_config.domEnd[JDIR] - run_config.domBeg[JDIR];
  double Lz = run_config.domEnd[KDIR] - run_config.domBeg[KDIR];
  double kvec[3];
  double xvec[3];

  // ------------------------------------
  // 0. Define rotation shift vectors
  // from the domain size.
  // Assume jshift_y = kshift_z = 1
  // and square cells (dx=dy=dz).
  // ------------------------------------

  jshift_y = 1;
  jshift_x = round(-jshift_y/Ly);

  kshift_z = 1;
  kshift_x = round(-kshift_z/Lz);
  RotateSet(jshift_x, jshift_y, kshift_x, kshift_z, run_config);

  eps = run_config.inputParam[EPS];   // wave amplitude 
  kx  = 2.0*CONST_PI/1.0;    // Wave vector is more conveniently expressed 
  ky  = 2.0*CONST_PI/Ly;     // in the rotated frame   
  kz  = 2.0*CONST_PI/Lz;

  phi = DIM_EXPAND(kx*x, + ky*y, + kz*z);  // k.x is invariant under rotations 

  // ------------------------------------
  // 1. Define 1D solution
  // ------------------------------------

  v[RHO] = 1.0;
  v[PRS] = run_config.inputParam[PR0];
  v[TRC] = 0.0;
  v[VX1] = run_config.inputParam[VEL0];  // translational velocity 
  v[VX2] = fabs(eps)*sin(phi);
  v[VX3] = fabs(eps)*cos(phi);
  v[BX1] = 1.0;
  v[BX2] = eps*sqrt(v[RHO])*sin(phi);
  v[BX3] = eps*sqrt(v[RHO])*cos(phi);

  // ------------------------------------
  // 2. Rotate vectors
  // ------------------------------------

  RotateVector(v + VX1,-1);
  RotateVector(v + BX1,-1);
  
  kvec[0] = kx; kvec[1] = ky; kvec[2] = kz;
  xvec[0] =  x; xvec[1] =  y; xvec[2] =  z;
  RotateVector(kvec,1);
  RotateVector(xvec,1);

  v[AX1] = 0.0;
  v[AX2] = eps*sin(phi)/kvec[0];
  v[AX3] = eps*cos(phi)/kvec[0] + xvec[1];  // Assume rho = B = 1 

  RotateVector(v+AX1,-1);
}

////////////////////////////////////////////////////////////////////////
/*!
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
  
*/
////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
{
}

////////////////////////////////////////////////////////////////////////
/*!
  Perform runtime data analysis.
  
  \param [in] d the PLUTO Data structure
  \param [in] grid   pointer to array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
{
  int i, j, k, nv;

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  double By_max, By2;

  DIM_EXPAND(double Lx = grid->xend_glob[IDIR] - grid->xbeg_glob[IDIR];, 
             double Ly = grid->xend_glob[JDIR] - grid->xbeg_glob[JDIR];,
             double Lz = grid->xend_glob[KDIR] - grid->xbeg_glob[KDIR];)

  double Vol_tot = DIM_EXPAND(Lx, *Ly, *Lz);
  double v[NVAR];
  double time = data->Dts->time;

  FILE *fp;
  RBox box;
  
  // ----------------------------------------------
  // 0. Define Box, initialize variables
  // ----------------------------------------------

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box);
  By_max = By2 = 0.0;

  // ----------------------------------------------
  // 1. Reduce data on single processor
  // (here we focus on decay)
  // ----------------------------------------------

  BOX_LOOP(box,i,j,k) {
    double dV = grid->dV(i, j, k);

    #ifdef HIGH_ORDER
    NVAR_LOOP(nv) v[nv] = data->Vpc(i,j,k,nv);
    #else
    NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k,nv);
    #endif
    By_max  = MAX(fabs(v[BX3]), By_max);
    By2    += v[BX3]*v[BX3]*dV;
  }

  // ----------------------------------------------
  // 2. Parallel reduction
  // ----------------------------------------------

  #ifdef PARALLEL
  double gscrh;
  MPI_Allreduce (&By_max, &gscrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  By_max = gscrh;

  MPI_Allreduce (&By2, &gscrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  By2 = gscrh;
  #endif

  By2 /= Vol_tot;

  if (g_rank == 0){
    if (g_stepNumber == 0) {
      fp = fopen ("cp_alfven_decay.dat","w");
      fprintf (fp,"# ----------------------------------------------------");
      fprintf (fp,"------------------------------------------------------\n");
      fprintf (fp,"# CP Alfven\n");
      fprintf (fp,"# Resolution = %d X %d X %d\n", RuntimeGet()->npoint[IDIR],
                                                   RuntimeGet()->npoint[JDIR],
                                                   RuntimeGet()->npoint[KDIR]);
    } else {
      fp = fopen ("cp_alfven_decay.dat","a");
    }

    fprintf (fp, "%12.6e  %12.6e  %12.6e\n", time, 
             By_max/data->run_config->inputParam[EPS], By2);
    fclose(fp);
  }

}

////////////////////////////////////////////////////////////////////////
/*! 
  Assign user-defined boundary conditions.
  
  \param [in,out] data  pointer to the PLUTO data structure
  \param [in]     box   pointer to a RBox structure containing the 
                        boundary box of the corresponding to side
  \param [in]     side  side of computational domain (special value side == 0 
                        is used to control INTERNAL_BOUNDARY)
  \param [in]     grid  pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
{
}
#if BODY_FORCE != NO
////////////////////////////////////////////////////////////////////////
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
////////////////////////////////////////////////////////////////////////
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
////////////////////////////////////////////////////////////////////////
double BodyForcePotential(const State& state, const RunConfig& run_config)
////////////////////////////////////////////////////////////////////////
{
  return 0.0;
}
#endif
