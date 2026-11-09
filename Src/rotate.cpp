////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Rotate a 1D problem

  Given the initial condition in the 1D, non-rotated frame (\f$\Sigma'\f$),
  we rotate quantities by an angle \f$\gamma \f$ around the y-axis and an
  angle \f$\alpha\f$  around the z-axis.
  The resulting transformation leaves scalar quantities invariant and produce
  vector rotations \f$\vec{q} = \tens{R}_{\gamma\alpha}\vec{q}'\f$ where
  \f[
    \tens{R}_{\gamma\alpha} = \left(\begin{array}{ccc}
       \cos\alpha\cos\gamma  &   -\sin\alpha  &   -\cos\alpha\sin\gamma
       \\[0.5em]
       \sin\alpha\cos\gamma  &    \cos\alpha  &   -\sin\alpha\sin\gamma
       \\[0.5em]
       \sin\gamma            &       0        &   \cos\gamma
    \end{array}\right)
    \,,\qquad
    \tens{R}^{-1}_{\gamma\alpha} = \left(\begin{array}{ccc}
       \cos\alpha\cos\gamma  &  \sin\alpha\cos\gamma  &  \sin\gamma
       \\[0.5em]
        -\sin\alpha          &   \cos\alpha           &    0
       \\[0.5em]
       -\cos\alpha\sin\gamma &  -\sin\alpha\sin\gamma  &  \cos\gamma
      \end{array}\right)
  \f]
  where \f$\tan\gamma = \cos\alpha\tan\beta\f$.
  In 2D, simply set \f$\gamma = 0\f$.

  In a discrete domain, the angles \f$\alpha\f$ and \f$\beta\f$ must be chosen
  in such a way that an integer shift of cells satisfies, for any flow
  quantity q, the translational invariance expressed by
  \f$ q(\vec{x}+\vec{s}) = q (\vec{x})\f$ where
  \f$\vec{s} = (n_x\Delta x, n_y \Delta y, n_z \Delta z)\f$ is orthogonal
   to \f$\vec{n}_p\f$:
  \f[
    \vec{n}_p\cdot\vec{s} = n_x\Delta x + n_y\Delta y\tan\alpha
                                        + n_z\Delta z\tan\beta = 0
  \f]
  The mesh spacing does need to be the same in all directions.
  For a given pair of rotation angles, there can be more than one set of
  integers that satisfies this equation.
  For instance, taking \f$\tan\alpha = 1,\;\tan\beta=2 \f$ one may choose
  (if mesh spacing is the same)
  \f$\vec{s} = \pm(1,1,-1) \f$, \f$\vec{s} = \pm (-1,1,0)\f$ (useful at a
  y-boundary) or  \f$\vec{s} = \pm (2,0,-1) \f$ (useful at a z-boundary).
  In practice, we find convenient to provide two pairs of integer numbers,
  \f$(j_x,\, j_y)\f$ and \f$ (k_x,\, k_z)\f$ giving the integer shifts
  along two directions at a time.
  These numbers are used to find the rotation angles:
  \f[
     \left\{\begin{array}{lcl}
      \DS j_x\Delta x + j_y\Delta y \tan\alpha &=& 0  \\[0.5em]
      \DS k_x\Delta x + k_z\Delta z \tan\beta  &=& 0  \,,\qquad
    \end{array}\right.
    \qquad\rightarrow\qquad
      \tan\alpha = -\frac{j_x}{j_y}\frac{\Delta x}{\Delta y} ,\quad
      \tan\beta  = -\frac{k_x}{k_z}\frac{\Delta x}{\Delta z}
  \f]
  This choice allows us to assign boundary conditions using only 2 shifts
  (rather than 3): one along the x direction and the other in the direction
  normal to the boundary plane.

  <b>Plane of Discontinuity.</b>
  A plane of discontinuity has equation \f$x' = 0\f$ in the non-rotated frame.
  Using the inverse transformation we find the equation of the plane in the
  rotated \f$\Sigma\f$-plane:
  \f[
    x' = \left(\tens{R}^{-1}_{\gamma} \vec{x}\right)_x
       =   (\cos\alpha\cos\gamma)\, x
         + (\sin\alpha\sin\gamma)\, y
         + (\sin\gamma) z = 0
    \quad\rightarrow\quad
    \boxed{x + y\tan\alpha + z\tan\beta = 0}
  \f]
  The normal to this plane is the vector \f$\vec{n}_p = (1, \tan\alpha, \tan\beta)\f$.

  <b>Note on periodic domains.</b>
  In a periodic domain, an integer number of wavelengths must be contained
  in each direction.
  Hence the rotation is usuallt applied by specifying
  \f$\tan\alpha = k_y/k_x\f$ and \f$\tan\beta = k_z/k_x\f$ which express
  the ratios between the \f$y\f$- and  \f$z\f$- components of the wavevector
  with the \f$x\f$ component.
  For one wavelength in each direction,
  \f$k_x = 2\pi/L_x\f$, \f$k_y = 2\pi/L_y\f$, \f$k_z = 2\pi/L_z\f$.
  The integer shifts can then be specified using
  \f[
      \tan\alpha = -\frac{j_x}{j_y}\frac{\Delta x}{\Delta y} = \frac{L_x}{L_y} ,
      \qquad
      \tan\beta  = -\frac{k_x}{k_z}\frac{\Delta x}{\Delta z} = \frac{L_x}{L_z}
  \f]


  <b>Note on the final time.</b>
  Usually, the rotated domain size differs from the typical 1D configuration
  since solutions are compared by looking at the profiles in the x direction.
  The final time is then adjusted in such a way that the intersection
  of the rotated plane with the x axis travels the same distance covered
  in the 1D case.
  Imagine a point in 1D that travels a distance \f$\delta L_1 = v\delta t_1\f$:
  in 2- or 3-D the same point will be traveling along the direction given by
  \f$\DS \vec{x}_p(t) = \hat{\bf n}_p vt\f$ where
  \f$ \hat{\bf n}_p = \vec{n}_p/|\vec{n}_p| \f$.
  In the rotated problem the solution is constant on the plane of equation
  \f$ \hat{\bf n}_p\cdot(\vec{x}-\vec{x}_p(t)) = 0\f$ and its intersection
  with the x axis is (\f$ y = z = 0 \f$) covers the distance \f$\delta L_1\f$
  \f[
      \hat{n}_x\delta L_1 = v\delta t
      \qquad\rightarrow\qquad
      \delta t = n_x\delta t_1 =
      \frac{\delta t_1}{\sqrt{1 + \tan^2\alpha + \tan^2\beta}}
  \f]

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Feb 27, 2025

  \b Reference:
     - [MT10] "A second-order unsplit Godunov scheme for cell-centered MHD: The CTU-GLM scheme"
        Mignone \& Tzeferacos, JCP (2010) 229, 2217
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "rotate.hpp"

static int    rot_jx, rot_jy, rot_kx, rot_kz;
static double rot_ta, rot_tb;
#if ACC_INIT
  #pragma acc declare create(rot_ta, rot_tb)
  #pragma acc declare create(rot_jx, rot_jy, rot_kx, rot_kz)
#endif

////////////////////////////////////////////////////////////////////////
double RotateGet_ta(void)
/*!
  Return the rotation angle tan(al.hppa)
*/
////////////////////////////////////////////////////////////////////////
{
  return rot_ta;
}

////////////////////////////////////////////////////////////////////////
double RotateGet_tb(void)
/*!
  
  Return the rotation angle tan(beta)
*/
////////////////////////////////////////////////////////////////////////
{
  return rot_tb;
}

////////////////////////////////////////////////////////////////////////
void RotateSet(int jx, int jy, int kx, int kz, RunConfig &run_config)
/*!
  For a given choice of the input shifts, initialize the rotation
  angles rot_ta \f$ = \tan\alpha\f$ and rot_tb \f$ = \tan\beta \f$.
*/
////////////////////////////////////////////////////////////////////////
{
  // static int first_call = 1;
  int    nx = run_config.Nx1;
  int    ny = run_config.Nx2;
  int    nz = run_config.Nx3;
  double ta;
  double tb;
  // double tg, cg, sg;

  double dx = (run_config.domEnd[IDIR] - run_config.domBeg[IDIR])/nx;
  double dy = (run_config.domEnd[JDIR] - run_config.domBeg[JDIR])/ny;
  double dz = (run_config.domEnd[KDIR] - run_config.domBeg[KDIR])/nz;

  // --------------------------------------------------------
  // 1. Define rotation angles in terms of integer shifts. 
  // Compute coordinates in non-rotated frame (x1, y1, z1).
  // --------------------------------------------------------

  if (jy == 0 || kz == 0){
    // printLog ("! RotateSet(): specify a non-zero value for jy and kz\n");
    // QuitPLUTO(" ",1);
  }
  #if DIMENSIONS == 1
  jx = kx = 0;
  #elif DIMENSIONS == 2
  kx = 0;
  #endif

  ta = -(double)jx/(double)jy*dx/dy;
  // ca = 1.0/sqrt(1.0 + ta*ta);
  // sa = ta*ca;

  tb = -(double)kx/(double)kz*dx/dz;
  // cb = 1.0/sqrt(1.0 + tb*tb);
  // sb = tb*cb;

  // tg = tb*ca;
  // cg = 1.0/sqrt(1.0 + tg*tg);
  // sg = cg*tg;

  rot_ta = ta;
  rot_tb = tb;

  rot_jx = jx; rot_jy = jy;
  rot_kx = kx; rot_kz = kz;

  // if (first_call){
  //   print ("> RotateSet(): ta = %f; tb = %f\n", ta, tb);
  //   print ("               tg = %f; sg = %f, cg = %f\n",tg,sg,cg);
  //   print ("               jx, jy = %d, %d\n",jx,jy);
  //   print ("               kx, kz = %d, %d\n",kx,kz);
  //   first_call = 0;
  // }
}
////////////////////////////////////////////////////////////////////////
void RotateVector(double *v, int s)
/*!
  Rotate a vector <v[0], v[1], v[2]>;
  - s =  1 --> clockwise rotation
  [transform vector to the 1D unrotated frame]
  - s = -1 --> counter-clockwise rotation
  [transform vector to the 2/3D frame]
*/
////////////////////////////////////////////////////////////////////////
{
  double vx = v[0];
  double vy = v[1];
  double vz = v[2];
  double ta = rot_ta;
  double tb = rot_tb;
  double sa, ca, tg, cg, sg;

// -- Compute tangent, sin() and cos() of al.hppa and beta -- 

  ca = 1.0/sqrt(ta*ta + 1.0);
  // cb = 1.0/sqrt(tb*tb + 1.0);
  sa = ta*ca;
  // sb = tb*cb;

  tg = tb*ca;
  cg = 1.0/sqrt(tg*tg + 1.0);
  sg = cg*tg;

// -- Rotate -- 

  if (s == 1){    // Clockwise rotation 
    v[0] =  vx*ca*cg + vy*sa*cg + vz*sg;
    v[1] = -vx*sa    + vy*ca;
    v[2] = -vx*ca*sg - vy*sa*sg + vz*cg;
  }else if (s == -1){ // Counter-clockwise rotation 
    v[0] = vx*ca*cg - vy*sa - vz*ca*sg;
    v[1] = vx*sa*cg + vy*ca - vz*sa*sg;
    v[2] = vx*sg            + vz*cg;
  }else{
    // print ("! RotateVector(): invalid value of s = %d\n", s);
    // QuitPLUTO(" ",1);
  }
}

////////////////////////////////////////////////////////////////////////
void RotateBoundary (const Ary4D &Vc,const Ary3D (&Vs)[6], RBox *box,
                     int side, Grid *grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int nv, i, j, k;
  int i1,j1,k1;
  int s;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  // --------------------------------------------------------
  // 2a. Set BC on X2_BEG
  // Note: Staggered B are assigned using the same box as 
  // cell centered variables, meaning that BX1 & BX3 will
  // never have the rightmost staggered value set correctly.
  // This is done for conveniency and also because using 
  // rotated BC requires more ghost zones than usual. 
  // --------------------------------------------------------

  if (side == X2_BEG){  // -- X2_BEG boundary -- 

    #pragma acc parallel loop collapse(3) seq present(Vc, Vs, grid) async(1)
    for (k = kbeg; k <= kend; k++){
      for (j = jend; j >= jbeg; j--){
        for (i = ibeg; i <= iend; i++){
          if (rot_jy == 0){
            i1 = i; j1 = grid->lbeg[JDIR];
          }else{
            s  = (rot_jy > 0 ? 1:-1);
            i1 = MAX(i + s*rot_jx, 0); i1 = MIN(i1, nx1_tot-1);
            j1 = MAX(j + s*rot_jy, 0); j1 = MIN(j1, nx2_tot-1);
          }

          NVAR_LOOP(nv) Vc(i,j,k,nv) = Vc(i1,j1,k,nv);

          #ifdef STAGGERED_MHD
          #if INCLUDE_IDIR
          CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i1,j1,k);  ,
                    Vs[EX1s](i,j,k) = Vs[EX1s](i1,j1,k);)
          #endif
          #if INCLUDE_JDIR
          CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i1,j1,k);  ,
                    Vs[EX2s](i,j,k) = Vs[EX2s](i1,j1,k);)
          #endif
          #if INCLUDE_KDIR
          CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i1,j1,k); ,
                    Vs[EX3s](i,j,k) = Vs[EX3s](i1,j1,k);)
          #endif
          #endif // STAGGERED_MHD 
        }
      }
    }
    #pragma acc wait
  }  

  // --------------------------------------------------------
  // 2b. Set BC on X2_END. See note @2a.
  // --------------------------------------------------------

  if (side == X2_END){  // -- X2_END boundary --   

    #pragma acc parallel loop collapse(3) seq present(Vc, Vs, grid) async(1)
    for (k = kbeg; k <= kend; k++){
      for (j = jbeg; j <= jend; j++){
        for (i = ibeg; i <= iend; i++){
          if (rot_jy == 0){
            i1 = i; j1 = grid->lend[JDIR];
          }else{
            s  = (rot_jy > 0 ? 1:-1);
            i1 = MAX(i - s*rot_jx, 0); i1 = MIN(i1, nx1_tot-1);
            j1 = MAX(j - s*rot_jy, 0); j1 = MIN(j1, nx2_tot-1);
          }

          NVAR_LOOP(nv) Vc(i,j,k,nv) = Vc(i1,j1,k,nv);

          #ifdef STAGGERED_MHD
          #if INCLUDE_IDIR
          CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i1,j1,k);  ,
                    Vs[EX1s](i,j,k) = Vs[EX1s](i1,j1,k);)
          #endif
          #if INCLUDE_JDIR
          CT_ASSIGN(Vs[BX2s](i,j+1,k) = Vs[BX2s](i1,j1+1,k);  ,
                    Vs[EX2s](i,j+1,k) = Vs[EX2s](i1,j1+1,k);)
          #endif
          #if INCLUDE_KDIR
          CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i1,j1,k); ,
                    Vs[EX3s](i,j,k) = Vs[EX3s](i1,j1,k);)
          #endif
         #endif
        }
      }
    }
    #pragma acc wait
  }

  // --------------------------------------------------------
  // 3a. Set BC on X3_BEG. See note @2a.
  // --------------------------------------------------------

  if (side == X3_BEG){  // -- X3_BEG boundary -- 

    #pragma acc parallel loop collapse(3) seq present(Vc, Vs, grid) async(1)
    for (k = kend; k >= kbeg; k--){
      for (j = jbeg; j <= jend; j++){
        for (i = ibeg; i <= iend; i++){
          if (rot_kz == 0){
            i1 = i; k1 = grid->lbeg[KDIR];
          }else{
            s  = (rot_kz > 0 ? 1:-1);
            i1 = MAX(i + s*rot_kx, 0); i1 = MIN(i1, nx1_tot-1);
            k1 = MAX(k + s*rot_kz, 0); k1 = MIN(k1, nx3_tot-1);
          }
          
          NVAR_LOOP(nv) Vc(i,j,k,nv) = Vc(i1,j,k1,nv);

          #ifdef STAGGERED_MHD
          #if INCLUDE_IDIR
          CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i1,j,k1);  ,
                    Vs[EX1s](i,j,k) = Vs[EX1s](i1,j,k1);)
          #endif
          #if INCLUDE_JDIR
          CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i1,j,k1);  ,
                    Vs[EX2s](i,j,k) = Vs[EX2s](i1,j,k1);)
          #endif
          #if INCLUDE_KDIR
          CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i1,j,k1); ,
                    Vs[EX3s](i,j,k) = Vs[EX3s](i1,j,k1);)
          #endif
          #endif // STAGGERED_MHD  
	
        }
      }
    } 
    #pragma acc wait
  }

  // --------------------------------------------------------
  // 3b. Set BC on X3_END. See note @2a.
  // --------------------------------------------------------

  if (side == X3_END){  // -- X3_END boundary -- 

    #pragma acc parallel loop collapse(3) seq present(Vc, Vs, grid) async(1)
    for (k = kbeg; k <= kend; k++){
      for (j = jbeg; j <= jend; j++){
        for (i = ibeg; i <= iend; i++){
          if (rot_kz == 0){
            i1 = i; k1 = grid->lend[KDIR];
          }else{
            s  = (rot_kz > 0 ? 1:-1);
            i1 = MAX(i - s*rot_kx, 0); i1 = MIN(i1, nx1_tot-1);
            k1 = MAX(k - s*rot_kz, 0); k1 = MIN(k1, nx3_tot-1);
          }

          NVAR_LOOP(nv) Vc(i,j,k,nv) = Vc(i1,j,k1,nv);

          #ifdef STAGGERED_MHD
          #if INCLUDE_IDIR
          CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i1,j,k1);  ,
                    Vs[EX1s](i,j,k) = Vs[EX1s](i1,j,k1);)
          #endif
          #if INCLUDE_JDIR
          CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i1,j,k1);  ,
                    Vs[EX2s](i,j,k) = Vs[EX2s](i1,j,k1);)
          #endif
          #if INCLUDE_KDIR
          CT_ASSIGN(Vs[BX3s](i,j,k+1) = Vs[BX3s](i1,j,k1+1); ,
                    Vs[EX3s](i,j,k+1) = Vs[EX3s](i1,j,k1+1);)
          #endif
          #endif  // STAGGERED_MHD 
        }
      }
    }
    #pragma acc wait
  }

  #ifndef STAGGERED_MHD
  #endif
}
