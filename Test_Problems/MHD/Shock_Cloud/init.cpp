////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief MHD Shock Cloud interaction.

  Set the initial condition for the 2D or 3D shock-cloud problem.
  The shock propagates in the x-direction and it is initially located
  at \c x=0.6 with post- and -pre shock value equal to
  \f[
    \left(\begin{array}{l}
      \rho \\[0.5em]
       v_x  \\[0.5em]
       p   \\[0.5em]
       B_y  \\[0.5em]
       B_z
    \end{array}\right)
     =
    \left(\begin{array}{l}
        3.86859      \\[0.5em]
          0          \\[0.5em]
        167.345      \\[0.5em]
        B_{\rm post} \\[0.5em]
       -B_{\rm post} \\[0.5em]
    \end{array}\right)
     \quad\mathrm{for}\quad x < 0.6\,,\qquad\qquad
    \left(\begin{array}{l}
      \rho \\[0.5em]
       v_x  \\[0.5em]
       p   \\[0.5em]
       B_y  \\[0.5em]
       B_z
    \end{array}\right)
     =
    \left(\begin{array}{l}
          1         \\[0.5em]
        -11.2536    \\[0.5em]
         1          \\[0.5em]
        B_{\rm pre} \\[0.5em]
       -B_{\rm pre} \\[0.5em]
    \end{array}\right)
     \quad\mathrm{for}\quad x > 0.6\,,\qquad
  \f]
  while the remaining vector components are 0.
  The cloud has radius \f$R\f$ centered at
  \f$r_0 = (x_0,y_0,z_0) = (0.8, 0.5, 0.5)\f$ and larger density
  \f$\rho = 10\f$.

  The interaction may be divided into two phases: 1) the collapse stage
  where the front of the cloud is strongly compressed and two fast shocks
  are generated and 2) the reexpansion phase which begins when the
  transmitted fast shock overtake the back of the cloud.

  The runtime parameters that are read from \c pluto.ini are
  - <tt>run_config.inputParam[B_POST]</tt>: sets the post-shock magnetic field;
  - <tt>run_config.inputParam[B_PRE]</tt>:  sets the pre-shock magnetic field;
  - <tt>run_config.inputParam[RADIUS]</tt>: sets the radius of the cloud;

  A list of the available configurations is given in the following table:

  <CENTER>
  Conf.| GEOMETRY  |DIM| T. STEPPING|RECONSTRUCTION|divB|  AMR
  -----|-----------|---|------------|--------------| ---|-------
   #01 |CARTESIAN  | 2 |  RK2       |  LINEAR      | CT | NO
   #02 |CARTESIAN  | 2 |  ChTr      |  LINEAR      | CT | NO
   #03 |CARTESIAN  | 2 |  HANCOCK   |  LINEAR      | 8W | NO
   #04 |CARTESIAN  | 3 |  ChTr      |  LINEAR      | CT | NO
   #05 |CARTESIAN  | 3 |  RK3       |  LINEAR      | CT | NO
   #06 |CARTESIAN  | 3 | ChTr       |  LINEAR      | GLM| NO
   #07 |CARTESIAN  | 3 |  RK3       |  WENO3_FD    | GLM| NO
   #08 |CARTESIAN  | 3 | HANCOK     |  LINEAR      | GLM| YES
   #09 |CARTESIAN  | 2 |  RK2       |  LINEAR      | CT | NO
  </CENTER>

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 17, 2014

  \b References:
     - "Interactions between magnetohydrodynamical shocks and denser clouds"
        W. Dai an P.R. Woodward, ApJ (1994) 436, 776.
     - "The divB = 0 Constraint in Shock-Capturing MHD codes"
        G. Toth, JCP 161, 605-652 (2000)
     - "A central constrained transport scheme for ideal magnetohydrodynamics"
        U. Ziegler, JCP 196 (2004), 393
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double r;
  double x0, y0, z0;

  x0 = 0.8;
  y0 = 0.5;
  z0 = 0.5;

  if (x < 0.6) {
    v[RHO] = 3.86859;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
#if PHYSICS==MHD
    v[BX1] = 0.0;
    v[BX2] =  run_config.inputParam[B_POST];
    v[BX3] = -run_config.inputParam[B_POST];
#endif
    v[PRS] = 167.345;
  }else{
    v[RHO] = 1.0;
    v[VX1] = -11.2536;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
#if PHYSICS==MHD
    v[BX1] = 0.0;
    v[BX2] = run_config.inputParam[B_PRE];
    v[BX3] = run_config.inputParam[B_PRE];
#endif
    v[PRS] = 1.0;

  }

  //  ----  CLOUD  ----  

  r = DIM_EXPAND((x-x0)*(x-x0), + (y-y0)*(y-y0) , + (z-z0)*(z-z0));

  if (sqrt(r) < run_config.inputParam[RADIUS]) v[RHO] = 10.0;

// no need for potential vector,
   // since CT_VEC_POT_INIT is set to NO
// 
  // v[AX1] = x3*v[BX2] - x2*v[BX3];
  // v[AX2] = 0.0;
  // v[AX3] = 0.0;
// 

  #if !INCLUDE_KDIR
  #endif
}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*!
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
  
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int  i, j, k, nv;
  RunConfig &run_config = *(data->run_config);
  #ifdef HIGH_ORDER 
  const Ary4D& Vc      = data->Vpc;  
  const Ary3D (&Vs)[6] = data->Vps;
  #else 
  const Ary4D& Vc      = data->Vc;  
  #ifdef STAGGERED_MHD
  const Ary3D (&Vs)[6] = data->Vs;
  #endif
  #endif

  if (side == X1_END){ // -- select the boundary side -- 

    #pragma acc parallel loop collapse(3) present(data, run_config,box) 
    BOX_LOOP(*box,i,j,k) {
      double v[NVAR];
      
      v[RHO] =  1.0;
      v[VX1] = -11.2536;
      v[VX2] =  0.0;
      v[VX3] =  0.0;
      v[PRS] =  1.0;
      #if PHYSICS == MHD
      v[BX1] =  0.0;
      v[BX2] =  run_config.inputParam[B_PRE];
      v[BX3] =  run_config.inputParam[B_PRE];
      #ifdef GLM_MHD
      v[PSI_GLM] = 0.0;
      #endif
      #endif
    
      NVAR_LOOP(nv) Vc(i,j,k,nv) = v[nv];
      
      // ----------------------------------------------------
      // Define conservative variables (at 2nd order level). 
      // Useful only during a Boundary_Conservative() call.
      // ----------------------------------------------------

      #ifdef HIGH_ORDER
      double u[NVAR];
      PrimToCons(v, u, run_config);
      NVAR_LOOP(nv) data->Uc(i,j,k,nv) = u[nv];
      #endif
    }

    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(data, box)
    BOXJ_LOOP(*box, i,j,k) {
      Vs[BX2s](i,j,k) =  run_config.inputParam[B_PRE];
      #ifdef HIGH_ORDER 
      data->Vs[BX2s](i,j,k) = data->Vps[BX2s](i,j,k);
      #endif
    }
    #endif // STAGGERED_MHD

    #if (DIMENSIONS == 3) && (defined STAGGERED_MHD)
    int iend = grid->lend[IDIR];
    int nx1_tot = grid->ntot[IDIR];
    int nx2_tot = grid->ntot[JDIR];
    int nx3_tot = grid->ntot[KDIR];

    #pragma acc parallel loop collapse(3) present(data)
    for (k = 0;      k < nx3_tot+1; k++){
    for (j = 0;      j < nx2_tot;   j++){
    for (i = iend+1; i < nx1_tot;   i++){
      Vs[BX3s](i,j,k) = run_config.inputParam[B_PRE];
      #ifdef HIGH_ORDER 
      data->Vs[BX3s](i,j,k) = data->Vps[BX3s](i,j,k);
      #endif
    }}}
    #else
    #endif

  }
}
