/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shock finding algorithm.

  Search and flag computational zones lying in a shock wave.
  The flagging strategy is based on two switches designed to detect
  the presence of compressive motion or shock waves in the fluid:
  \f[
     \nabla\cdot\vec{v} < 0 \qquad{\rm and}\qquad
     \Delta x\frac{|\nabla p|}{p} > \epsilon_p
  \f]
  where \f$\epsilon_p\f$ sets the shock strength.
  At the discrete level we replace the two conditions by
  \f[
    \sum_d \frac{ A_{\vec{i}+\HALF\hvec{e}_d}v_{d,\vec{i}+\HALF\hvec{e}_d}
                 -A_{\vec{i}-\HALF\hvec{e}_d}v_{d,\vec{i}-\HALF\hvec{e}_d} }
                {\Delta{\cal V}_{d,\vec{i}}}  < 0
                \qquad{\rm and}\qquad
    \sum_{d} \left|p_{\vec{i}+\hvec{e}_d} - p_{\vec{i}-\hvec{e}_d}\right|
             <
     \epsilon_p \min_d\left(p_{\vec{i}+\hvec{e}_d},
                            p_{\vec{i}-\hvec{e}_d},p_{\vec{i}}\right)
  \f]
  where \f$\hvec{i} = (i,j,k)\f$ is a vector of integer numbers
  giving the position of a computational zone, while \f$\hvec{e}_d =
  (\delta_{1d},\delta_{2d},\delta_{3d})\f$ is a unit vector in the direction
  given by \c d.
  Once a zone has been tagged as lying in a shock, different flags may be
  switched on or off to control the update strategy in these critical regions.

  This function can be called called when:
  - \c SHOCK_FLATTENING has been set to \c MULTID: in this case shocked zones
    are tagged with \c FLAG_MINMOD and \c FLAG_HLL that will later
    be used to force the reconstruction with the minmod limiter
    and the Riemann solver to HLL.
  - \c ENTROPY_SWITCH has been turned on: this flag will be checked later in
    the ConsToPrim() functions in order to recover pressure from the
    entropy density rather than from the total energy density.
    The update process is:

    - start with a vector of primitive variables  <tt> {V,s} </tt>
      where \c s is the entropy;
    - set boundary condition on \c {V};
    - compute \c {s} from \c {V};
    - flag zones where entropy may be used (flag = 1);
    - evolve equations for one time step;
    - convert <tt> {U,S} </tt> to primitive:
      \code
        if (flag == 0) {  // Use energy
          p = p(E)
          s = s(p)
        }else{            // use entropy
          p = p(S)
          E = E(p)
        }
        \endcode

  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    March 10, 2026
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"gpluto.hpp"

#ifndef EPS_PSHOCK_FLATTEN
 #define EPS_PSHOCK_FLATTEN  5.0
#endif

#ifndef EPS_PSHOCK_ENTROPY
 #define EPS_PSHOCK_ENTROPY  0.05
#endif

#if (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
/* *************************************************************** */
/*!
 * \param [in,out] d     pointer to data structure
 * \param [in]     grid  pointer to grid structure
 *
 ***************************************************************** */
void FlagShock (Data *data, Grid *grid, RBox *box)
{
  int  i, j, k;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  const Ary3D &Ax1 = grid->A[IDIR];
  const Ary3D &Ax2 = grid->A[JDIR];
  const Ary3D &Ax3 = grid->A[KDIR];

  static Ary3D pt;

  RBox boxt;
  RunConfig &run_config = *(data->run_config);
  const Ary4D &Vc = data->Vc;

  // 0. Define boxes, allocate memory

  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &boxt);  

  if (pt.isAllocated() == 0) {
    pt.create(nx1_tot,nx2_tot,nx3_tot);
    #pragma acc enter data copyin(pt)
  }

  // ------------------------------------------------------
  // 1. Compute pressure and flag, initially all
  //    zones with ENTROPY_SWITCH.
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(data, run_config, pt) \
                                        copyin(boxt)
  BOX_LOOP(boxt,i,j,k){
    #if HAVE_ENERGY
    pt(i,j,k) = Vc(i,j,k,PRS);
    #else
    double a  = run_config.isoSoundSpeed;
    pt(i,j,k) = Vc(i,j,k,RHO)*a*a;
    #endif
    #if (ENTROPY_SWITCH == SELECTIVE) || (ENTROPY_SWITCH == ALWAYS)
    data->flag(i,j,k) |= FLAG_ENTROPY;
    #endif
  }

  // ------------------------------------------------------
  // 2. For user-defined boundary, enlarge the active box 
  //    by one zone into the interior (this requires at 
  //    least 2 ghost zones)
  // ------------------------------------------------------

  if (1){
    int set_bc = 0, dir;
    if (grid->lbound[IDIR] == USERDEF) ibeg -= INCLUDE_IDIR;
    if (grid->lbound[JDIR] == USERDEF) jbeg -= INCLUDE_JDIR;
    if (grid->lbound[KDIR] == USERDEF) kbeg -= INCLUDE_KDIR;

    if (grid->rbound[IDIR] == USERDEF) iend += INCLUDE_IDIR;
    if (grid->rbound[JDIR] == USERDEF) jend += INCLUDE_JDIR;
    if (grid->rbound[KDIR] == USERDEF) kend += INCLUDE_KDIR;
  }

  // ------------------------------------------------------
  // 3. Track zones lying in a shock and activate or  
  //    deactivate flag bits.
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(data, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    double divv, gradp, pt_min;
    double dpx1, pt_min1, dvx1;
    double dpx2, pt_min2, dvx2;
    double dpx3, pt_min3, dvx3;
    
    // 3A. Compute divergence of velocity 

    #if GEOMETRY == CARTESIAN
    DIM_EXPAND(dvx1 = (Vc(i+1,j,k,VX1) - Vc(i-1,j,k,VX1))/grid->dx[IDIR][i];   ,
               dvx2 = (Vc(i,j+1,k,VX2) - Vc(i,j-1,k,VX2))/grid->dx[JDIR][j];   ,
               dvx3 = (Vc(i,j,k+1,VX3) - Vc(i,j,k-1,VX3))/grid->dx[KDIR][k];)
    divv = DIM_EXPAND(dvx1, + dvx2, + dvx3);
    #else
    DIM_EXPAND(dvx1 =   Ax1(i+1,j,k)*(Vc(i+1,j,k,VX1) + Vc(i,j,k,VX1))
                      - Ax1(i,j,k)  *(Vc(i-1,j,k,VX1) + Vc(i,j,k,VX1));,

               dvx2 =   Ax2(i,j+1,k)*(Vc(i,j+1,k,VX2) + Vc(i,j,k,VX2))
                      - Ax2(i,j,k)  *(Vc(i,j-1,k,VX2) + Vc(i,j,k,VX2));,

               dvx3 =   Ax3(i,j,k+1)*(Vc(i,j,k+1,VX3) + Vc(i,j,k,VX3))
                      - Ax3(i,j,k)  *(Vc(i,j,k-1,VX3) + Vc(i,j,k,VX3));)
    divv = (DIM_EXPAND(dvx1, + dvx2, + dvx3))/grid->dV(i,j,k);
    #endif

    // // 3B. Flatten on axis 

    // #if SHOCK_FLATTENING_AXIS == YES
    // #if GEOMETRY == CYLINDRICAL
    // if (grid->lbound[IDIR] != 0) {
    //   if (i == ibeg || i == ibeg-1) data->flag(i,j,k) |= FLAG_FLAT;
    // }
    // #elif GEOMETRY == SPHERICAL
    // if (grid->lbound[JDIR] != 0) {
    //   if (j == jbeg || j == jbeg-1) data->flag(i,j,k) |= FLAG_FLAT;
    // }
    // #endif
    // #endif

    if (divv < 0.0){

      // 3C. Compute undivided difference of the total
      //     pressure and minimum value in neighbour zones

      pt_min = pt(i,j,k);
      DIM_EXPAND(pt_min1 = MIN(pt(i+1,j,k), pt(i-1,j,k));  ,
                 pt_min2 = MIN(pt(i,j+1,k), pt(i,j-1,k));  ,
                 pt_min3 = MIN(pt(i,j,k+1), pt(i,j,k-1)); )

      DIM_EXPAND(pt_min = MIN(pt_min, pt_min1);  ,
                 pt_min = MIN(pt_min, pt_min2);  ,
                 pt_min = MIN(pt_min, pt_min3);)

      DIM_EXPAND(dpx1 = fabs(pt(i+1,j,k) - pt(i-1,j,k));  ,
                 dpx2 = fabs(pt(i,j+1,k) - pt(i,j-1,k));  ,
                 dpx3 = fabs(pt(i,j,k+1) - pt(i,j,k-1));)

      gradp = DIM_EXPAND(dpx1, + dpx2, + dpx3);

      // 3D. Flag shocked zones and expand the buffer by
      //     one zone to the left and to the right.

      #if SHOCK_FLATTENING == MULTID
      if (gradp > EPS_PSHOCK_FLATTEN*pt_min) {
        uint32_t *f   = &data->flag(i,j,k);
        #pragma acc atomic update
        *f   |= FLAG_HLL;

        int i1,j1,k1;
        int imin = MAX(i-INCLUDE_IDIR,0), imax = MIN(i+INCLUDE_IDIR,nx1_tot-1);
        int jmin = MAX(j-INCLUDE_JDIR,0), jmax = MIN(j+INCLUDE_JDIR,nx2_tot-1);
        int kmin = MAX(k-INCLUDE_KDIR,0), kmax = MIN(k+INCLUDE_KDIR,nx3_tot-1);

        #if INCLUDE_KDIR
        #pragma acc loop seq
        for (k1 = kmin; k1 <= kmax; k1++){
          uint32_t *f = &data->flag(i,j,k1);
          #pragma acc atomic update
          *f |= FLAG_MINMOD;
        }
        #endif
        #if INCLUDE_JDIR
        #pragma acc loop seq
        for (j1 = jmin; j1 <= jmax; j1++){
          uint32_t *f = &data->flag(i,j1,k);
          #pragma acc atomic update
          *f |= FLAG_MINMOD;
        }
        #endif
        #if INCLUDE_IDIR
        #pragma acc loop seq
        for (i1 = imin; i1 <= imax; i1++){
          uint32_t *f = &data->flag(i1,j,k);
          #pragma acc atomic update
          *f |= FLAG_MINMOD;
        }
        #endif

      }
      #endif /* SHOCK_FLATTENING == MULTID */

      // 3E. When using entropy, unflag those zones lying
      //     in a shock as well as one neighbour cells to the
      //     left and to the right for each dimension.

      #if ENTROPY_SWITCH == SELECTIVE
      if (gradp > EPS_PSHOCK_ENTROPY*pt_min) { /* -- unflag zone -- */
        DIM_EXPAND(uint32_t *f   = &data->flag(i,j,k);
                   uint32_t *fim = &data->flag(i-1,j,k);
                   uint32_t *fip = &data->flag(i+1,j,k); ,
                   uint32_t *fjm = &data->flag(i,j-1,k);
                   uint32_t *fjp = &data->flag(i,j+1,k); ,
                   uint32_t *fkm = &data->flag(i,j,k-1);
                   uint32_t *fkp = &data->flag(i,j,k+1);)

        #pragma acc atomic update
        *f   |= FLAG_ENERGY;

        #pragma acc atomic update
        *fip |= FLAG_ENERGY;
        #pragma acc atomic update
        *fim |= FLAG_ENERGY;

        #if INCLUDE_JDIR
        #pragma acc atomic update
        *fjp |= FLAG_ENERGY;
        #pragma acc atomic update
        *fjm |= FLAG_ENERGY;
        #endif

        #if INCLUDE_KDIR
        #pragma acc atomic update
        *fkp |= FLAG_ENERGY;
        #pragma acc atomic update
        *fkm |= FLAG_ENERGY;
        #endif

      }
      #endif

    } // End if (divV < 0)

  }}} // End for(i,j,k)

  // ------------------------------------------------------
  // 4. Apply b.c. on flag. Particular useful to preserve 
  //    correct axial symmetry across r = 0 (cylindrical) 
  //    or theta = 0 (spherical coordinates).
  // ------------------------------------------------------

  // if (1){
  //   int set_bc = 0, dir;

  //   DIM_LOOP(dir) {
  //     int lbound = grid->lbound[dir];
  //     int rbound = grid->rbound[dir];
  //     if (   (lbound == PERIODIC) 
  //         || (lbound == REFLECTIVE) 
  //         || (lbound == AXISYMMETRIC)
  //         || (lbound == EQTSYMMETRIC)) set_bc = 1;
  //     if (   (rbound == PERIODIC) 
  //         || (rbound == REFLECTIVE) 
  //         || (rbound == AXISYMMETRIC)
  //         || (rbound == EQTSYMMETRIC)) set_bc = 1;
  //   }
  //   if (set_bc) FlagBoundary(data->flag, grid);
  // }

  // // 5. Exchange flag across processors

  // #ifdef PARALLEL
  // MPI_Barrier (MPI_COMM_WORLD);
  // ParallelExchangeFlagNew<FLAG_OR>(data->flag, grid);
  // #endif

  FlagBoundary(data->flag, FLAG_OR, grid);

  // 6. Unflag entropy where FLAG_ENERGY is enabled

  #if ENTROPY_SWITCH == SELECTIVE
  BOX_LOOP(boxt,i,j,k){
    if (data->flag(i,j,k) & FLAG_ENERGY) {
      data->flag(i,j,k) &= ~(FLAG_ENTROPY);
    } 
  }
  #endif

}
#undef  EPS_PSHOCK_FLATTEN
#undef  EPS_PSHOCK_ENTROPY
#endif // (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
