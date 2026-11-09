////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Compute the curl of a vector.

  Computes the curl (J) of a vector (V) using a cross-like stencil: 

  - when V is edge-centered, J will be face-centered
  - when V is face-centered, J will be edge-centered
  - when V is cell-centered, J will be also face-centered

  - In Cartesian geometry:
    \f[
    B_x = \pd{A_z}{y} - \pd{A_y}{z} \,,\quad
    B_y = \pd{A_x}{z} - \pd{A_z}{x} \,,\quad
    B_z = \pd{A_y}{x} - \pd{A_x}{y}
    \f]

  - In Cylindrical geometry:
    \f[
    B_R    = \frac{1}{R}\pd{A_z}{\phi} - \pd{A_\phi}{z}  \,,\quad
    B_\phi = \pd{A_R}{z} - \pd{A_z}{R}             \,,\quad
    B_z    = \frac{1}{R}\pd{(RA_\phi)}{R} - \frac{1}{R}\pd{A_R}{\phi}
    \f]

  - In spherical geometry:
    \f[
    B_r      =   \frac{1}{r\sin\theta}\pd{(\sin\theta A_\phi)}{\theta}
               - \frac{1}{r\sin\theta}\pd{A_\theta}{\phi}           \,,\quad
    B_\theta =   \frac{1}{r\sin\theta}\pd{A_r}{\phi}
               - \frac{1}{r}\pd{(rA_\phi)}{r}                         \,,\quad
    B_\phi   =   \frac{1}{r}\pd{(rA_\theta)}{r}
               - \frac{1}{r}\pd{A_r}{\theta}
    \f]

  For instance, when V == A (vector potential) in CT-MHD, the convention 
  is the following:

  - <tt> Ax1[i,j,k] --> Ax1(0++) </tt> \f$= A_{x_1, i, j+\HALF, k+\HALF}\f$
  - <tt> Ax2[i,j,k] --> Ax2(+0+) </tt> \f$= A_{x_2, i+\HALF, j, k+\HALF}\f$
  - <tt> Ax3[i,j,k] --> Ax3(++0) </tt> \f$= A_{x_3, i+\HALF, j+\HALF, k}\f$

  \verbatim
            ______________________
           /                     /|
          /                     / |
         /       z face        /  |
        /                     Ax  |
       /                     /    |
      /                     /     Az
     ----------Ay-----------      |
     |                     |      |
     |                     |   y  |
     |                     | face |
     |                     |     /
     Az       x face      Az    /
     |                     |   Ax
     |                     |  /
     |                     | /
     |                     |/
     ----------Ay-----------
  \endverbatim

  \author  A. Mignone (andrea.mignone@unito.it)
  \date    Nov 08, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"


#define SET_COORDINATES(x1c, x2c, x3c, x1, x2, x3)  x1c = x1; \
                                                    x2c = x2; \
                                                    x3c = x3; 

/**
 * \brief Compute the curl of a staggered or cell-centered vector field.
 *
 * This routine computes the curl of a generic vector field **V = (Vx1, Vx2, Vx3)** 
 * and stores the result in **J = (Jx1, Jx2, Jx3)**. 
 * It supports both cell-centered and face/edge-centered fields using PLUTO/gPLUTO's
 * staggered indexing conventions.  
 *
 * The function is used for:
 *   - **B = curl(A)** when the vector potential **A** is edge-centered
 *   - **J = curl(B)** when the magnetic field **B** is face-centered
 *   - generic applications where V may be cell-centered (location = CELL_CENTER)
 *
 * Depending on the output array location_, backward-difference indices 
 * '(im, jm, km)' refer to:
 *   - **face/edge-centered V** -> backward index equals forward index (im = i)
 *   - **cell-centered V**      ->  backward index shifted by 1 (im = i-1)
 *
 * The routine is geometry-aware and supports:
 *   - Cartesian
 *   - Cylindrical (r, phi, z)
 *   - Spherical (r, theta, phi)
 *
  \param[in]  V     Array of three Ary3D objects containing the components
 *                   of the input vector field V.
  \param[out] J     Array of three Ary3D objects where the components of
 *                   curl(V) will be stored.
  \param[in]  grid  Pointer to grid structure containing metric, spacing, and 
 *                   precomputed cell-interface coordinates (x1p, x1m, ...).
  \param[in]  box   Pointer to RBox defining the loop bounds (ibeg,iend,...)
 *
 * \return None (results written into *J).
 */
#if (RESISTIVITY == EXPLICIT) || (ASSIGN_VECTOR_POTENTIAL == YES)
void VectorCurl(Ary3D *V, Ary3D *J, Grid *grid, RBox *box)
{
  int  i,  j, k;
  DIM_EXPAND(int ip;, int jp;, int kp;)
  DIM_EXPAND(int im;, int jm;, int km;)

  int  inc_p, inc_m;  // Increments for central differencing

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  int location;
  
  DIM_EXPAND(double *dx1 = grid->dx[IDIR];,
             double *dx2 = grid->dx[JDIR];,
             double *dx3 = grid->dx[KDIR];) 

  #if GEOMETRY != CARTESIAN
  double **xC = grid->xC, **xL = grid->xL, **xR = grid->xR;
  double *x1p, *x1m, *x1c;
  double *x2p, *x2m, *x2c;
  double *x3p, *x3m, *x3c;
  #endif

  double Vx1_p = 0.0, Vx2_p = 0.0, Vx3_p = 0.0;
  double Vx1_m = 0.0, Vx2_m = 0.0, Vx3_m = 0.0;
  double curlV, dt = 1.0;


  Ary3D &Vx1 = V[IDIR], &Jx1 = J[IDIR];
  Ary3D &Vx2 = V[JDIR], &Jx2 = J[JDIR];
  Ary3D &Vx3 = V[KDIR], &Jx3 = J[KDIR];

  // ------------------------------------------------------
  // 0. Check Compatibility
  // ------------------------------------------------------

  // ------------------------------------------------------
  // 1. Compute curl(V), x1 component: 
  //
  //    curl(V)_x = dVz/dy - dVy/dz
  //
  //              Vy_p                 
  //               o                   z^ 
  //               |                    |
  //    Vz_m       |         Vz_p       |
  //     o---------c---------o          +----> y
  //               |            
  //               |            
  //               o
  //              Vy_m          
  //      
  //    V   |      V_plus          |      V_minus         |  Jx1 = curlV 
  //  ------+----------------------+----------------------+-------------
  //  (A = edge-centered)                                  (B = face-centered)
  //    Ay  | (i-1/2, j,     k+1/2)| (i-1/2, j,     k-1/2)| (i-1/2, j, k)
  //    Az  | (i-1/2, j+1/2, k)    | (i-1/2, j-1/2, k)    | (i-1/2, j, k)
  //  ------+----------------------+----------------------+-------------
  //  (B = face-centered)                                  (J = edge-centered)
  //    By  | (i, j-1/2, k)        | (i, j-1/2, k-1)      | (i, j-1/2, k-1/2)
  //    Bz  | (i, j,     k-1/2)    | (i, j-1,   k-1/2)    | (i, j-1/2, k-1/2)
  //  ------+----------------------+----------------------+-------------
  //  (B = cell-centered)                                  (J = cell-centered)
  //    By  | (i, j,   k+1)        | (i, j,   k-1)        | (i, j, k)
  //    Bz  | (i, j+1, k)          | (i, j-1, k)          | (i, j, k)
  //  ------+----------------------+----------------------+-------------

  // 1a. Determine coordinates location for Jx1 component.
  //     Here {x1c, x2c, x3c} are the coordinates of Jx1, while 
  //     {x1p, x2p, x3p} and {x1m, x2m, x3m} are the coordinates used in 
  //     differentiating V.

  inc_p = inc_m = 0;

  if (Jx1.isAllocated()){

    location = Jx1.location_;
    dt       = 1.0;
    ibeg = box->ibeg, iend = box->iend;
    jbeg = box->jbeg, jend = box->jend;
    kbeg = box->kbeg, kend = box->kend;

    if (location == X1_FACE) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xL[IDIR], xC[JDIR], xC[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xL[IDIR], xL[JDIR], xL[KDIR]);
      SET_COORDINATES (x1p, x2p, x3p, xL[IDIR], xR[JDIR], xR[KDIR]);
      #endif

      iend += INCLUDE_IDIR;  // Increment by one the staggered index
      inc_p = 1;
    }else if (location == X1_EDGE) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xC[IDIR], xL[JDIR],   xL[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xC[IDIR], xC[JDIR]-1, xC[KDIR]-1);
      SET_COORDINATES (x1p, x2p, x3p, xC[IDIR], xC[JDIR],   xC[KDIR]);
      #endif

      jend += INCLUDE_JDIR; // Increment by one the staggered indices
      kend += INCLUDE_KDIR; // Increment by one the staggered indices
    
      inc_m = 1;
    }else if (location == CELL_CENTER) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xC[IDIR],   xC[JDIR],   xC[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xC[IDIR]-1, xC[JDIR]-1, xC[KDIR]-1);
      SET_COORDINATES (x1p, x2p, x3p, xC[IDIR]+1, xC[JDIR]+1, xC[KDIR]+1);
      #endif

      dt    = 0.5;
      inc_p = inc_m = 1;
    }else{
      QuitPLUTO("VectorCurl(): Invalid location_ for Jx1",1);
    }

    // 1b. Compute Jx1

    #pragma acc parallel loop collapse(3) present(V, J, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){

      // 1b1. Define forward and backward indices to compute central derivatives. 

      DIM_EXPAND(                                 ,
                 jp = j + inc_p; jm = j - inc_m;  ,
                 kp = k + inc_p; km = k - inc_m;
      )

      DIM_EXPAND(                                           ,
                 Vx3_p = Vx3(i,jp,k); Vx3_m = Vx3(i,jm,k);  ,
                 Vx2_p = Vx2(i,j,kp); Vx2_m = Vx2(i,j,km);
      )

      // 1b2. Compute x1-component of curl(V) in different geometries 

      #if ((GEOMETRY == CARTESIAN) && (PHYSICS != GRMHD)) || \
          ((PHYSICS == GRMHD) && (GR_METRIC == CARTESIAN)) 

      DIM_EXPAND(                         ,
                 double dtdx2 = dt/dx2[j];, 
                 double dtdx3 = dt/dx3[k];)

      curlV = DIM_EXPAND(0.0, 
                         + dtdx2*(Vx3_p - Vx3_m), 
                         - dtdx3*(Vx2_p - Vx2_m));

      #elif GEOMETRY == CYLINDRICAL
      
      DIM_EXPAND(                                        , 
                 double dtdx2 = dt/(fabs(x1c[i])*dx2[j]);, 
                 double dtdx3 = dt/dx3[k];)

      // HOT FIX: avoid singularities on axis.
      #if INCLUDE_JDIR
      if (x1c[i] == 0.0) dtdx2 = 0.0;
      #endif
      
      curlV = DIM_EXPAND(   0.0, + dtdx2*(Vx3_p - Vx3_m),
                                 - dtdx3*(Vx2_p - Vx2_m));

      #elif GEOMETRY == SPHERICAL

      double dV2   = fabs(cos(x2m[j]) - cos(x2p[j]));
      DIM_EXPAND(                               ,
                 double dtdx2 = dt/(x1c[i]*dV2);, 
                 double dtdx3 = dt*dx2[j]/(x1c[i]*dV2*dx3[k]);)

      double Ax2p = fabs(sin(x2p[j]));
      double Ax2m = fabs(sin(x2m[j]));

      curlV = DIM_EXPAND(0.0, + dtdx2*(Ax2p*Vx3_p - Ax2m*Vx3_m),
                              - dtdx3*(Vx2_p - Vx2_m));
      #endif

      Jx1(i,j,k) = curlV;
    }}}
  }

  // -------------------------------------------------------------------
  // 2. Compute curl(V), x2 component
  //
  //    curl(V)_y = dVx/dz - dVz/dx
  //
  //              Vz_p                 
  //               o                   x^ 
  //               |                    |
  //    Vx_m       |         Vx_p       |
  //     o---------+---------o          +----> z
  //               |            
  //               |            
  //               o
  //              Vz_m          
  //      
  //    V   |      V_plus          |      V_minus         |  Jx2 = curlV 
  //  ------+----------------------+----------------------+-----------------
  //  (A = edge-centered)                                  (B = face-centered)
  //    Az  | (i+1/2, j-1/2, k)    | (i-1/2, j-1/2, k)    | (i, j-1/2, k)
  //    Ax  | (i,     j-1/2, k+1/2)| (i,     j-1/2, k-1/2)| (i, j-1/2, k)
  //  ------+----------------------+----------------------+-----------------
  //  (B = face-centered)                                  (J = edge-centered)
  //    Bx  | (i-1/2, j,     k)    | (i-1/2, j,     k-1)  | (i-1/2, j, k-1/2)
  //    Bz  | (i+1/2, j,     k-1/2)| (i-1/2, j,     k-1/2)| (i-1/2, j, k-1/2)
  //  ------+----------------------+----------------------+-----------------
  //  (B = cell-centered)                                  (J = cell-centered)
  //    Bx  | (i,   j, k+1)        | (i,   j, k-1)        | (i, j, k)
  //    Bz  | (i+1, j, k)          | (i-1, j, k)          | (i, j, k)
  //  ------+----------------------+----------------------+-----------------

  // 2a. Determine coordinates location Jx2 component.
  //     Here {x1c, x2c, x3c} are the coordinates of Jx2, while 
  //     {x1p, x2p, x3p} and {x1m, x2m, x3m} are the coordinates used in 
  //     differentiating V.

  if (Jx2.isAllocated()){
    location = Jx2.location_;
    dt       = 1.0;
    ibeg = box->ibeg, iend = box->iend;
    jbeg = box->jbeg, jend = box->jend;
    kbeg = box->kbeg, kend = box->kend;

    if (location == X2_FACE) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xC[IDIR], xL[JDIR], xC[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xL[IDIR], xL[JDIR], xL[KDIR]);
      SET_COORDINATES (x1p, x2p, x3p, xR[IDIR], xL[JDIR], xR[KDIR]);
      #endif

      jend += INCLUDE_JDIR;  // Increment by one the staggered index
      inc_p = 1;
    }else if (location == X2_EDGE) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xL[IDIR],   xC[JDIR], xL[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xC[IDIR]-1, xC[JDIR], xC[KDIR]-1);
      SET_COORDINATES (x1p, x2p, x3p, xC[IDIR],   xC[JDIR], xC[KDIR]);
      #endif

      iend += INCLUDE_IDIR; // Increment by one the staggered indices
      kend += INCLUDE_KDIR; // Increment by one the staggered indices
      inc_m = 1;

    }else if (location == CELL_CENTER) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xC[IDIR],   xC[JDIR],   xC[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xC[IDIR]-1, xC[JDIR]-1, xC[KDIR]-1);
      SET_COORDINATES (x1p, x2p, x3p, xC[IDIR]+1, xC[JDIR]+1, xC[KDIR]+1);
      #endif

      dt  = 0.5;
      inc_p = inc_m = 1;

    }else{
      QuitPLUTO("VectorCurl(): Invalid location_ for Jx2",1);
    }

    // 2b. Compute Jx2

    #pragma acc parallel loop collapse(3) present(V, J, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){

      // 2b1. Define forward and backward indices to compute central derivatives. 

      DIM_EXPAND(ip = i + inc_p; im = i - inc_m;  ,
                                                  , 
                 kp = k + inc_p; km = k - inc_m;
      )

      DIM_EXPAND(Vx3_p = Vx3(ip,j,k); Vx3_m = Vx3(im,j,k);  ,
                                                            ,
                 Vx1_p = Vx1(i,j,kp); Vx1_m = Vx1(i,j,km);
      )

      // 2b2. Compute x2-component of curl(V) in different geometries 

      #if GEOMETRY == CARTESIAN
      DIM_EXPAND(double dtdx1 = dt/dx1[i];, 
                                          ,
                 double dtdx3 = dt/dx3[k];)
      
      curlV = DIM_EXPAND(- dtdx1*(Vx3_p - Vx3_m),  
                                                , 
                         + dtdx3*(Vx1_p - Vx1_m));

      #elif GEOMETRY == CYLINDRICAL

      DIM_EXPAND(double dtdx1 = dt/dx1[i];, 
                                          ,
                 double dtdx3 = dt/dx3[k];)

      curlV =  DIM_EXPAND(- dtdx1*(Vx3_p - Vx3_m),  
                                                 ,
                          + dtdx3*(Vx1_p - Vx1_m));

      #elif GEOMETRY == SPHERICAL
      
      DIM_EXPAND(double dtdx1 = dt/(x1c[i]*dx1[i]);    , 
                                                       ,
                 double Ax2 = fabs(sin(x2c[j]));
                 double dtdx3 = dt/(x1c[i]*Ax2*dx3[k]);)

      curlV = DIM_EXPAND(- dtdx1*(x1p[i]*Vx3_p - x1m[i]*Vx3_m),
                                                              ,
                         + dtdx3*(Vx1_p - Vx1_m));
      #endif

      Jx2(i,j,k) = curlV;
    }}}
  }

  // ------------------------------------------------------
  // 3. Compute curl(V), x3 component
  //
  //    curl(V)_z = dVy/dx - dVx/dy
  //
  //              Vx_p                 
  //               o                   y^ 
  //               |                    |
  //    Vy_m       |         Vy_p       |
  //     o---------+---------o          +----> x
  //               |            
  //               |            
  //               o
  //              Vx_m          
  //      
  //    V   |      V_plus          |      V_minus         |  Jx3 = curlV 
  //  ------+----------------------+----------------------+-----------------
  //  (A = edge-centered)                                  (B = face-centered)                                
  //    Ax  | (i,     j+1/2, k-1/2)| (i,     j-1/2, k-1/2)| (i,   j,   k-1/2)
  //    Ay  | (i+1/2, j,     k-1/2)| (i-1/2, j,     k-1/2)| (i,   j,   k-1/2)
  //  ------+----------------------+----------------------+-----------------
  //  (B = face-centered)                                  (J = edge-centered)
  //    By  | (i,     j-1/2, k)    | (i-1,   j-1/2, k)    | (i-1/2, j-1/2, k)
  //    Bx  | (i-1/2, j,     k)    | (i-1/2, j-1,   k)    | (i-1/2, j-1/2, k)
  //  ------+----------------------+----------------------+-----------------
  //  (B = cell-centered)                                  (J = cell-centered)
  //    By  | (i+1, j,   k)        | (i-1, j,   k)        | (i, j, k)   
  //    Bx  | (i,   j+1, k)        | (i,   j-1, k)        | (i, j, k)
  //  ------+----------------------+----------------------+-----------------
  
  // 3a. Determine coordinates location for Jx3 component.
  //     Here {x1c, x2c, x3c} are the coordinates of Jx3, while 
  //     {x1p, x2p, x3p} and {x1m, x2m, x3m} are the coordinates used in 
  //     differentiating V.

  if (Jx3.isAllocated()){
    location = Jx3.location_;
    dt       = 1.0;

    ibeg = box->ibeg, iend = box->iend;
    jbeg = box->jbeg, jend = box->jend;
    kbeg = box->kbeg, kend = box->kend;

    if (location == X3_FACE) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xC[IDIR], xC[JDIR], xL[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xL[IDIR], xL[JDIR], xL[KDIR]);
      SET_COORDINATES (x1p, x2p, x3p, xR[IDIR], xR[JDIR], xL[KDIR]);
      #endif

      kend += INCLUDE_KDIR;  // Increment by one the staggered index
      inc_p = 1;
    }else if (location == X3_EDGE) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xL[IDIR],   xL[JDIR],   xC[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xC[IDIR]-1, xC[JDIR]-1, xC[KDIR]);
      SET_COORDINATES (x1p, x2p, x3p, xC[IDIR],   xC[JDIR],   xC[KDIR]);
      #endif

      iend += INCLUDE_IDIR; // Increment by one the staggered indices
      jend += INCLUDE_JDIR; // Increment by one the staggered indices
      inc_m = 1;

    }else if (location == CELL_CENTER) {
      #if GEOMETRY != CARTESIAN
      SET_COORDINATES (x1c, x2c, x3c, xC[IDIR],   xC[JDIR],   xC[KDIR]);
      SET_COORDINATES (x1m, x2m, x3m, xC[IDIR]-1, xC[JDIR]-1, xC[KDIR]-1);
      SET_COORDINATES (x1p, x2p, x3p, xC[IDIR]+1, xC[JDIR]+1, xC[KDIR]+1);
      #endif

      inc_p = inc_m = 1;
      dt  = 0.5;
    }else{
      QuitPLUTO("VectorCurl(): Invalid location_ for Jx3",1);
    }
  
    // 3b. Compute Jx3

    #pragma acc parallel loop collapse(3) present(V, J, grid)
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){

      // 3b1. Define forward and backward indices to compute central derivatives. 

      DIM_EXPAND(ip = i + inc_p; im = i - inc_m;  ,
                 jp = j + inc_p; jm = j - inc_m;  ,
      )

      DIM_EXPAND(Vx2_p = Vx2(ip,j,k); Vx2_m = Vx2(im,j,k);  ,
                 Vx1_p = Vx1(i,jp,k); Vx1_m = Vx1(i,jm,k);  ,
      )

      // 3b2. Compute x3-component of curl(V) in different geometries 

      #if GEOMETRY == CARTESIAN

      DIM_EXPAND(double dtdx1 = dt/dx1[i];, 
                 double dtdx2 = dt/dx2[j];, );

      curlV =  DIM_EXPAND(  dtdx1*(Vx2_p - Vx2_m), 
                          - dtdx2*(Vx1_p - Vx1_m), );

      #elif GEOMETRY == CYLINDRICAL

      double dV1 = fabs(x1c[i])*dx1[i];
      DIM_EXPAND(double dtdx1 = dt/dV1;,
                 double dtdx2 = dt/(x1c[i]*dx2[j]);,
                 );

      curlV =  DIM_EXPAND( dtdx1*(fabs(x1p[i])*Vx2_p - fabs(x1m[i])*Vx2_m) ,
                         - dtdx2*(Vx1_p - Vx1_m), );
                         
      #elif GEOMETRY == SPHERICAL

      DIM_EXPAND(double dtdx1 = dt/(x1c[i]*dx1[i]);, 
                 double dtdx2 = dt/(x1c[i]*dx2[j]);, );

      curlV = DIM_EXPAND(  dtdx1*(x1p[i]*Vx2_p - x1m[i]*Vx2_m),
                         - dtdx2*(Vx1_p - Vx1_m) ,
                                                 );
      #endif

      Jx3(i,j,k) = curlV;
    }}}
  }
}
#endif
