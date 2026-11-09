// ///////////////////////////////////////////////////////////////////// 
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
    \sum_d \frac{ A_{\vec{i}+\HALF\hat{\bf e}_d}v_{d,\vec{i}+\HALF\hat{\bf e}_d}
                 -A_{\vec{i}-\HALF\hat{\bf e}_d}v_{d,\vec{i}-\HALF\hat{\bf e}_d} }
                {\Delta{\cal V}_{d,\vec{i}}}  < 0
                \qquad{\rm and}\qquad
    \sum_{d} \left|p_{\vec{i}+\hat{\bf e}_d} - p_{\vec{i}-\hat{\bf e}_d}\right|
             <
     \epsilon_p \min_d\left(p_{\vec{i}+\hat{\bf e}_d},
                            p_{\vec{i}-\hat{\bf e}_d},p_{\vec{i}}\right)
  \f]
  where \f$\hat{\bf i} = (i,j,k)\f$ is a vector of integer numbers
  giving the position of a computational zone, while \f$\hat{\bf e}_d =
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
    -  compute \c {s} from \c {V};
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
  \date    Aug 28, 2025
*/
// ///////////////////////////////////////////////////////////////////// 
#include"gpluto.hpp"

#ifndef EPS_PSHOCK_FLATTEN
 #define EPS_PSHOCK_FLATTEN  5.0
#endif

#ifndef EPS_PSHOCK_ENTROPY
 #define EPS_PSHOCK_ENTROPY  0.05
#endif

#if (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
// *************************************************************** 
/*!
 * \param [in,out] d     pointer to data structure
 * \param [in]     grid  pointer to grid structure
 *
 ***************************************************************** */
void FlagShock (Data *d, Grid *grid, RBox *box)
{
  int  i, j, k;
  #if SHOCK_FLATTENING_AXIS == YES
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg;
  #endif
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];

  #if GEOMETRY != CARTESIAN
  const Ary3D &Ax1 = grid->A[IDIR];
  const Ary3D &Ax2 = grid->A[JDIR];
  const Ary3D &Ax3 = grid->A[KDIR];
  #endif

  static Ary3D pt;

  RBox boxt, boxtm;
  RunConfig &run_config = *(d->run_config);

  // 0. Define boxes, allocate memory

  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &boxt);  
  RBoxDefine (INCLUDE_IDIR, nx1_tot-1-INCLUDE_IDIR, 
              INCLUDE_JDIR, nx2_tot-1-INCLUDE_JDIR, 
              INCLUDE_KDIR, nx3_tot-1-INCLUDE_KDIR, &boxtm);

  if (pt.isAllocated() == 0) {
    pt.create(nx1_tot,nx2_tot,nx3_tot);
    #pragma acc enter data copyin(pt)
  }

  // ------------------------------------------------------
  // 1. Compute pressure and flag, initially all
  //    zones with ENTROPY_SWITCH.
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(d, run_config, pt) \
                                        copyin(boxt)
  BOX_LOOP(boxt,i,j,k){
    #if HAVE_ENERGY
    pt(i,j,k) = d->Vc(i,j,k,PRS);
    #else
    double a = run_config.isoSoundSpeed;
    pt(i,j,k) = d->Vc(i,j,k,RHO)*a*a;
    #endif

    #if (ENTROPY_SWITCH == SELECTIVE) || (ENTROPY_SWITCH == ALWAYS)
    d->flag(i,j,k) |= FLAG_ENTROPY;
    #endif
  }

  // ------------------------------------------------------
  // 1a. Track zones lying in a shock
  // ------------------------------------------------------

  const double rj = run_config.inputParam[JET_RAD];
  const double rj2 = rj*rj;

  #pragma acc parallel loop collapse(3) present(d, grid) copyin(boxtm)
  BOX_LOOP(boxtm,i,j,k){
    double divv, gradp, pt_min;
    double dpx1, pt_min1, dvx1;
    double dpx2, pt_min2, dvx2;
    double dpx3, pt_min3, dvx3;
    
    // 2a. Compute divergence of velocity 

    #if GEOMETRY == CARTESIAN
    DIM_EXPAND(dvx1 = (d->Vc(i+1,j,k,VX1) - d->Vc(i-1,j,k,VX1))/grid->dx[IDIR][i];   ,
               dvx2 = (d->Vc(i,j+1,k,VX2) - d->Vc(i,j-1,k,VX2))/grid->dx[JDIR][j];   ,
               dvx3 = (d->Vc(i,j,k+1,VX3) - d->Vc(i,j,k-1,VX3))/grid->dx[KDIR][k];)
    divv = DIM_EXPAND(dvx1, + dvx2, + dvx3);
    #else
    DIM_EXPAND(dvx1 =   Ax1(i+1,j,k)*(d->Vc(i+1,j,k,VX1) + d->Vc(i,j,k,VX1))
                      - Ax1(i,j,k)  *(d->Vc(i-1,j,k,VX1) + d->Vc(i,j,k,VX1));,

               dvx2 =   Ax2(i,j+1,k)*(d->Vc(i,j+1,k,VX2) + d->Vc(i,j,k,VX2))
                      - Ax2(i,j,k)  *(d->Vc(i,j-1,k,VX2) + d->Vc(i,j,k,VX2));,

               dvx3 =   Ax3(i,j,k+1)*(d->Vc(i,j,k+1,VX3) + d->Vc(i,j,k,VX3))
                      - Ax3(i,j,k)  *(d->Vc(i,j,k-1,VX3) + d->Vc(i,j,k,VX3));)
    divv = (DIM_EXPAND(dvx1, + dvx2, + dvx3))/grid->dV(i,j,k);
    #endif

    // 1a1. Flatten on axis 

    #if SHOCK_FLATTENING_AXIS == YES
    #if GEOMETRY == CYLINDRICAL
    if (grid->lbound[IDIR] != 0) {
      if (i == ibeg || i == ibeg-1) d->flag(i,j,k) |= FLAG_FLAT;
    }
    #elif GEOMETRY == SPHERICAL
    if (grid->lbound[JDIR] != 0) {
      if (j == jbeg || j == jbeg-1) d->flag(i,j,k) |= FLAG_FLAT;
    }
    #endif
    #endif

    double x1 = grid->xC[IDIR][i], x2 = grid->xC[JDIR][j], x3 = grid->xC[KDIR][k];


    double r2 = x1*x1 + x2*x2 + x3*x3;
    // bool customFlag =   (-2.0*rj < x1) && (x1 < 2.0*rj)
    //                  && (-2.0*rj < x2) && (x2 < 2.0*rj) 
    //                                    && (x3 < 2.0*rj);

    bool customFlag = r2 < 4.0*rj2;
    if (divv < 0.0 || customFlag){

      // 1a2. Compute undivided difference of the total
      //      pressure and minimum value in neighbour zones

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

      // 1a3. Flag shocked zones and expand the buffer by
      //      one zone to the left and to the right.
      //      IMPORTANT: do not try to expand by more than
      //      one zone or otherwise serial and parallel
      //      computations will give different results.
      //      If you really need to do so, thif function should
      //      be restructured as:
      //
      //      DOM_LOOP(k,j,i) flag[][][] =..
      //      MPI_Exchange(flag)
      //      BOX_LOOP(i,j,k) expand buffer

      #if SHOCK_FLATTENING == MULTID
      if (gradp > EPS_PSHOCK_FLATTEN*pt_min || customFlag) {
        {
        uint32_t *f   = &d->flag(i,j,k);
        #pragma acc atomic update
        *f   |= FLAG_HLL;
        }

        #if INCLUDE_KDIR
        int kmin = MAX(k-INCLUDE_KDIR,0), kmax = MIN(k+INCLUDE_KDIR,nx3_tot-1);
        #pragma acc loop seq
        for (int k1 = kmin; k1 <= kmax; k1++){
          uint32_t *f = &d->flag(i,j,k1);
          #pragma acc atomic update
          *f |= FLAG_MINMOD;
        }
        #endif
        #if INCLUDE_JDIR
        int jmin = MAX(j-INCLUDE_JDIR,0), jmax = MIN(j+INCLUDE_JDIR,nx2_tot-1);
        #pragma acc loop seq
        for (int j1 = jmin; j1 <= jmax; j1++){
          uint32_t *f = &d->flag(i,j1,k);
          #pragma acc atomic update
          *f |= FLAG_MINMOD;
        }
        #endif
        #if INCLUDE_IDIR
        int imin = MAX(i-INCLUDE_IDIR,0), imax = MIN(i+INCLUDE_IDIR,nx1_tot-1);
        #pragma acc loop seq
        for (int i1 = imin; i1 <= imax; i1++){
          uint32_t *f = &d->flag(i1,j,k);
          #pragma acc atomic update
          *f |= FLAG_MINMOD;
        }
        #endif
      }
      #endif // SHOCK_FLATTENING == MULTID 

      // 1b. When using entropy, we unflag those zones lying
      //     in a shock as well as one neighbour cells to the
      //     left and to the right for each dimension.

      #if ENTROPY_SWITCH == SELECTIVE
      if (gradp > EPS_PSHOCK_ENTROPY*pt_min) { // -- unflag zone -- 
        DIM_EXPAND(uint32_t *f   = &d->flag(i,j,k);
                   uint32_t *fim = &d->flag(i-1,j,k);
                   uint32_t *fip = &d->flag(i+1,j,k);,
                   uint32_t *fjm = &d->flag(i,j-1,k);
                   uint32_t *fjp = &d->flag(i,j+1,k);,
                   uint32_t *fkm = &d->flag(i,j,k-1);
                   uint32_t *fkp = &d->flag(i,j,k+1);)

        #pragma acc atomic update
        *f   &= ~(FLAG_ENTROPY);

        #pragma acc atomic update
        *fip &= ~(FLAG_ENTROPY);
        #pragma acc atomic update
        *fim &= ~(FLAG_ENTROPY);

        #if INCLUDE_JDIR
        #pragma acc atomic update
        *fjp &= ~(FLAG_ENTROPY);
        #pragma acc atomic update
        *fjm &= ~(FLAG_ENTROPY);
        #endif

        #if INCLUDE_KDIR
        #pragma acc atomic update
        *fkp &= ~(FLAG_ENTROPY);
        #pragma acc atomic update
        *fkm &= ~(FLAG_ENTROPY);
        #endif

      }
      #endif

    } // End if (divV < 0)
  }

  #ifdef PARALLEL
  ParallelExchangeFlag<FLAG_COPY>(d->flag, grid);
  #endif

}
#undef  EPS_PSHOCK_FLATTEN
#undef  EPS_PSHOCK_ENTROPY
#endif // (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
