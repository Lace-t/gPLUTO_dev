////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Update staggered magnetic field.

  Update face-centered magnetic field in the constrained transport
  formulation using a discrete version of Stoke's theorem.
  The update consists of a single Euler step:
  \f[
    \mathtt{Bs} \Longleftarrow \mathtt{Bs} + \Delta t R
  \f]
  where \c data->Vs is the main staggered array used by PLUTO,
  \c Bs is the magnetic field to be updated and \c R is the
  right hand side alreadx2[j] computed during the unsplit integrator.
  \c data->Vs and \c Bs may be the same array or may be different.

  \b References
   - "A staggered mesh algorithm using high-order Godunov fluxes to
      ensure solenoidal magnetic field in MHD simulations"\n
      Balsara \& Spicer, JCP (1999) 149, 270

  \author  A. Mignone (andrea.mignone@unito.it)
  \date    Feb 18, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
   Update staggered magnetic field using discrete version of
   Stoke's theorem.
   Only \c data->Vs is updated, while \c Bs is the original array:\n
   \c data->Vs = \c Bs + \c dt * \c R, where R = curl(E) is the electric field.
 
   \param [in,out] d     pointer to PLUTO Data structure.
   \param [in]     data->dVs    array of staggered variables to be updated.
   \param [in]     dt    step size
   \param [in]     grid  pointer to Grid structure
 */ 
////////////////////////////////////////////////////////////////////////
void CT_Update(const Data *data, Ary3D (&Rs)[6], double dt, Grid *grid, RBox* box)
{
  int  i, j, k;

  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  double *dx1 = grid->dx[IDIR];
  double *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];

  #if GEOMETRY != CARTESIAN  
  double *x1m = grid->xL[IDIR];
  double *x1p = grid->xR[IDIR];
  double *x1 = grid->xC[IDIR];
  #endif
  #if GEOMETRY == SPHERICAL
  double Ax2, dV2;
  double *x2p = grid->xR[JDIR], *x2m = grid->xL[JDIR];
  #endif

  EMF *emf = data->emf;
  auto &Ex1e = emf->Ee[IDIR];
  auto &Ex2e = emf->Ee[JDIR];
  auto &Ex3e = emf->Ee[KDIR];
  #if PHYSICS == ResRMHD
  auto &Bx1e = emf->Be[IDIR];
  auto &Bx2e = emf->Be[JDIR];
  auto &Bx3e = emf->Be[KDIR];
  #endif

  RBox tot_box;

  RBoxDefine (0, grid->ntot[IDIR]-1,
              0, grid->ntot[JDIR]-1,
              0, grid->ntot[KDIR]-1, &tot_box);

  mynvtxstart_("CT_Update",YELLOW);
  
  // ------------------------------------------------------
  // 0a. Set boundary (if enabled) and handle singularities
  // ------------------------------------------------------

  #if CT_EMF_BOUNDARY == YES
  CT_EMF_Boundary(data, grid);
  #endif

  #if (GEOMETRY == CYLINDRICAL) && (INCLUDE_JDIR == YES)
  {
    double Ax1 = grid->A[IDIR](ibeg,jbeg,kbeg);
    static const bool singularity =    (fabs(Ax1) < 1.e-12) 
                                    && (grid->lbound[IDIR] == SINGULAR_AXIS);

    if (singularity){ 
      #pragma acc parallel loop present(data, grid)
      for (k = kbeg; k <= kend; k++){

        int ib = grid->lbeg[IDIR];
        int nphi = grid->nint[JDIR];

        // Average Ez at the pole 

        double Ez = 0.0; 
        #pragma acc loop vector reduction(+:Ez)
        for (j = grid->lbeg[JDIR]; j <= grid->lend[JDIR]+1; j++) {
          Ez += Ex3e(ib,j,k);
        }
        Ez /= nphi+1;
        #pragma acc loop vector
        for (j = grid->lbeg[JDIR]-1; j <= grid->lend[JDIR]+1; j++) {
          DIM_EXPAND(                    ,
                      Ex3e(ib,j,k) = Ez;  ,
                      Ex2e(ib,j,k) = 0.0;)
        }
      }
    }
  }
  #elif (GEOMETRY == SPHERICAL) && (INCLUDE_KDIR == YES)
  {
    Ax2 = grid->A[JDIR](ibeg,jbeg,kbeg);
    static const bool north_pole =    (fabs(Ax2) < 1.e-12)
                                   && (grid->lbound[JDIR] == SINGULAR_AXIS);

    if (north_pole){ 
      #pragma acc parallel loop present(data, grid)
      for (i = ibeg; i <= iend; i++){

        int jb = grid->lbeg[JDIR];
        int nphi = grid->nint[JDIR];

        // Average Ez at the pole 

        double Er = 0.0; 
        #pragma acc loop vector reduction(+:Er)
        for (k = grid->lbeg[KDIR]; k <= grid->lend[KDIR]+1; k++) { 
          Er += Ex1e(i,jb,k);
        }
        Er /= nphi+1; 

        #pragma acc loop vector
        for (k = grid->lbeg[KDIR]-1; k <= grid->lend[KDIR]+1; k++) {
          Ex1e(i,jb,k) = Er;  
          Ex3e(i,jb,k) = 0.0;
        }
      }
    }

    Ax2 = grid->A[JDIR](ibeg,jend+1,kbeg);
    static const bool south_pole =    ( fabs(Ax2) < 1.e-12)
                                   && (grid->rbound[JDIR] == SINGULAR_AXIS);

    if (south_pole){ 
      #pragma acc parallel loop present(data, grid)
      for (i = ibeg; i <= iend; i++){

        int je = grid->lend[JDIR]+1;
        int nphi = grid->nint[JDIR];

        // Average Ez at the pole 

        double Er = 0.0; 
        #pragma acc loop vector reduction(+:Er)
        for (k = grid->lbeg[KDIR]; k <= grid->lend[KDIR]+1; k++) {
          Er += Ex1e(i,je,k);
        }
        Er /= nphi+1;
        #pragma acc loop vector
        for (k = grid->lbeg[KDIR]-1; k <= grid->lend[KDIR]+1; k++) {
          Ex1e(i,je,k) = Er;  
          Ex3e(i,je,k) = 0.0;
        }
      }
    }
  }
  #endif

  // ------------------------------------------------------
  // 1. Update Bx1, Ex1 at (i+1/2, j, k) faces
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) map(to:  Rs)
  #pragma acc parallel loop collapse(3) present(data, emf, grid, Rs)
  for (k = kbeg; k <= kend;              k++){
  for (j = jbeg; j <= jend;              j++){
  for (i = ibeg; i <= iend+INCLUDE_IDIR; i++){
    CT_ASSIGN(double curlE;,
              double curlB;)
    
    #if ((GEOMETRY == CARTESIAN) && (PHYSICS != GRMHD)) || ((PHYSICS == GRMHD) && (GR_METRIC == CARTESIAN)) 
    DIM_EXPAND(                         ,
               double dtdx2 = dt/dx2[j];,
               double dtdx3 = dt/dx3[k];)
    curlE = DIM_EXPAND(  0.0                                    ,
                       + dtdx2*(Ex3e(i,j+1,k) - Ex3e(i,j,k))    ,
                       - dtdx3*(Ex2e(i,j,k+1) - Ex2e(i,j,k)) );

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB = DIM_EXPAND(  0.0                                            ,
                       + dt/dx2[j]*(Bx3e(i,j+1,k) - Bx3e(i,j,k))    ,
                       - dt/dx3[k]*(Bx2e(i,j,k+1) - Bx2e(i,j,k)) );
    #endif

    #elif GEOMETRY == CYLINDRICAL
    double Ax1 = fabs(x1m[i]);

    curlE = DIM_EXPAND(   0.0                                           ,
                       + dt/(Ax1*dx2[j])*(Ex3e(i,j+1,k) - Ex3e(i,j,k))  ,
                       - dt/dx3[k]      *(Ex2e(i,j,k+1) - Ex2e(i,j,k)));

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB = DIM_EXPAND(   0.0                                           ,
                       + dt/(Ax1*dx2[j])*(Bx3e(i,j+1,k) - Bx3e(i,j,k))  ,
                       - dt/dx3[k]      *(Bx2e(i,j,k+1) - Bx2e(i,j,k)));
    #endif

#if RING_AVERAGE != NO || 1
if (fabs(x1m[i]) < 1.e-12 && i == ibeg) curlE = 0.0; // !!!! HOT FIX
#endif

    #elif GEOMETRY == SPHERICAL

    double Ax2p = fabs(sin(x2p[j]));
    double Ax2m = fabs(sin(x2m[j]));

    dV2 = fabs(cos(x2m[j]) - cos(x2p[j]));
    curlE = DIM_EXPAND(0.0                                           ,
                       + dt/(x1m[i]*dV2)*(  Ax2p*Ex3e(i,j+1,k) 
                                          - Ax2m*Ex3e(i,j,k))   ,
                       - dt*dx2[j]/(x1m[i]*dV2*dx3[k])*(  Ex2e(i,j,k+1) 
                                                        - Ex2e(i,j,k)));

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB = DIM_EXPAND(0.0                                           ,
                       + dt/(x1m[i]*dV2)*(  Ax2p*Bx3e(i,j+1,k) 
                                          - Ax2m*Bx3e(i,j,k))   ,
                       - dt*dx2[j]/(x1m[i]*dV2*dx3[k])*(  Bx2e(i,j,k+1) 
                                                        - Bx2e(i,j,k)));
    #endif

    #elif ((PHYSICS == GRMHD) && (GR_METRIC != CARTESIAN)) 
    double dtdA1 = dt/grid->A[IDIR](i,j,k);
    double dL2m = grid->L[JDIR](i,j,k);
    double dL2p = grid->L[JDIR](i,j,k+1);
    double dL3m = grid->L[KDIR](i,j,k);
    double dL3p = grid->L[KDIR](i,j+1,k);

    curlE = DIM_EXPAND(  0.0                                                                      ,
                       + dtdA1*(Ex3e(i,j+1,k)*dL3p - Ex3e(i,j,k)*dL3m)  ,
                       - dtdA1*(Ex2e(i,j,k+1)*dL2p - Ex2e(i,j,k)*dL2m));
    // curlE = DIM_EXPAND(  0.0                                       ,
    //                    + dtdA2[j]*(Ex3e(i,j+1,k) - Ex3e(i,j,k))    ,
    //                    - dtdA3[k]*(Ex2e(i,j,k+1) - Ex2e(i,j,k)) );

    #endif

    CT_ASSIGN(Rs[BX1s](i,j,k) = -curlE;  ,
              Rs[EX1s](i,j,k) = +curlB;)
  }}}

  // ------------------------------------------------------
  // 2. Update Bx2, Ex2 at (i, j+1/2, k) faces
  // ------------------------------------------------------

#if INCLUDE_JDIR

  #pragma omp target teams distribute parallel for collapse(3) map (to:  Rs)
  #pragma acc parallel loop collapse(3) present(data, emf, grid, Rs)
  for (k = kbeg; k <= kend;              k++){
  for (j = jbeg; j <= jend+INCLUDE_JDIR; j++){
  for (i = ibeg; i <= iend;              i++){

    CT_ASSIGN(double curlE;,
              double curlB;)

    #if (GEOMETRY == CARTESIAN) || ((PHYSICS == GRMHD) && (GR_METRIC == CARTESIAN)) 
    DIM_EXPAND(double dtdx1 = dt/dx1[i];,
                                        ,
               double dtdx3 = dt/dx3[k];)

    curlE = DIM_EXPAND(- dtdx1*(Ex3e(i+1,j,k) - Ex3e(i,j,k))   ,
                                                               ,
                       + dtdx3*(Ex1e(i,j,k+1) - Ex1e(i,j,k)) );

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB = DIM_EXPAND(- dt/dx1[i]*(Bx3e(i+1,j,k) - Bx3e(i,j,k))   ,
                                                                             ,
                       + dt/dx3[k]*(Bx1e(i,j,k+1) - Bx1e(i,j,k)) );
    #endif

    #elif GEOMETRY == CYLINDRICAL

    curlE =  DIM_EXPAND(- dt/dx1[i]*(Ex3e(i+1,j,k) - Ex3e(i,j,k))    ,
                                                                     ,
                        + dt/dx3[k]*(Ex1e(i,j,k+1) - Ex1e(i,j,k)));

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB =  DIM_EXPAND(- dt/dx1[i]*(Bx3e(i+1,j,k) - Bx3e(i,j,k))    ,
                                                                     ,
                        + dt/dx3[k]*(Bx1e(i,j,k+1) - Bx1e(i,j,k)));
    #endif

    #elif GEOMETRY == SPHERICAL

    Ax2 = fabs(sin(x2m[j]));
    curlE = DIM_EXPAND(- dt/(x1[i]*dx1[i])*(  x1p[i]*Ex3e(i+1,j,k) 
                                            - x1m[i]*Ex3e(i,j,k))           ,
                                                                            ,
                       + dt/(x1[i]*Ax2*dx3[k])*(Ex1e(i,j,k+1) - Ex1e(i,j,k)));

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB = DIM_EXPAND(- dt/(x1[i]*dx1[i])*(  x1p[i]*Bx3e(i+1,j,k) 
                                            - x1m[i]*Bx3e(i,j,k))           ,
                                                                            ,
                       + dt/(x1[i]*Ax2*dx3[k])*(Bx1e(i,j,k+1) - Bx1e(i,j,k)));    
    #endif

#if RING_AVERAGE != NO || 1
if (Ax2 < 1.e-12 && j == jend+1) curlE = 0.0; // !!!! HOT FIX
if (Ax2 < 1.e-12 && j == jbeg)   curlE = 0.0; 
#endif

    #elif ((PHYSICS == GRMHD) && (GR_METRIC != CARTESIAN)) 
    double dtdA2 = dt/grid->A[JDIR](i,j,k);
    // double dL1 = grid->L[IDIR];
    // double dL3 = grid->L[KDIR];
    double dL1m = grid->L[IDIR](i,j,k);
    double dL1p = grid->L[IDIR](i,j,k+1);
    double dL3m = grid->L[KDIR](i,j,k);
    double dL3p = grid->L[KDIR](i+1,j,k);

    curlE = DIM_EXPAND(- dtdA2*(Ex3e(i+1,j,k)*dL3(i+1,j,k) - Ex3e(i,j,k)*dL3(i,j,k))  , 
                         0.0                                                            ,
                       + dtdA2*(Ex1e(i,j,k+1)*dL1(i,j,k+1) - Ex1e(i,j,k)*dL1(i,j,k)));

    #endif

    CT_ASSIGN(Rs[BX2s](i,j,k) = -curlE;,
              Rs[EX2s](i,j,k) = +curlB;)
  }}}
#endif 

  // ------------------------------------------------------
  //  3. Update Bx3, Ex3 at (i, j, k+1/2) faces
  // ------------------------------------------------------

#if INCLUDE_KDIR

  #pragma omp target teams distribute parallel for collapse(3) map (to:Rs)
  #pragma acc parallel loop collapse(3) present(data, emf, grid, Rs)
  for (k = kbeg; k <= kend+INCLUDE_KDIR; k++){
  for (j = jbeg; j <= jend;              j++){
  for (i = ibeg; i <= iend;              i++){

    CT_ASSIGN(double curlE;,
              double curlB;)

    #if (GEOMETRY == CARTESIAN) || ((PHYSICS == GRMHD) && (GR_METRIC == CARTESIAN)) 
    DIM_EXPAND(double dtdx1 = dt/dx1[i];,
               double dtdx2 = dt/dx2[j];,
                                        )

    curlE =  DIM_EXPAND( dtdx1*(Ex2e(i+1,j,k) - Ex2e(i,j,k))  ,
                       - dtdx2*(Ex1e(i,j+1,k) - Ex1e(i,j,k))  ,
                       );

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)

    curlB =   dt/dx1[i]*(Bx2e(i+1,j,k) - Bx2e(i,j,k))
            - dt/dx2[j]*(Bx1e(i,j+1,k) - Bx1e(i,j,k));
    #endif

    #elif GEOMETRY == CYLINDRICAL
    double dV1 = fabs(x1[i])*dx1[i];
    curlE =  DIM_EXPAND(dt/dV1*(  fabs(x1p[i])*Ex2e(i+1,j,k) 
                                - fabs(x1m[i])*Ex2e(i,j,k))     ,
                       - dt/(x1[i]*dx2[j])*(  Ex1e(i,j+1,k) 
                                            - Ex1e(i,j,k));     ,
                       );

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB =  DIM_EXPAND(dt/dV1*(  fabs(x1p[i])*Bx2e(i+1,j,k) 
                                - fabs(x1m[i])*Bx2e(i,j,k))     ,
                       - dt/(x1[i]*dx2[j])*(  Bx1e(i,j+1,k) 
                                            - Bx1e(i,j,k));     ,
                       );
    #endif

    #elif GEOMETRY == SPHERICAL
    curlE = DIM_EXPAND(  dt/(x1[i]*dx1[i])*(  x1p[i]*Ex2e(i+1,j,k) 
                                            - x1m[i]*Ex2e(i,j,k)),
                       - dt/(x1[i]*dx2[j])*(  Ex1e(i,j+1,k) 
                                            - Ex1e(i,j,k)) ,
                      );

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    curlB = DIM_EXPAND(  dt/(x1[i]*dx1[i])*(  x1p[i]*Bx2e(i+1,j,k) 
                                            - x1m[i]*Bx2e(i,j,k)),
                       - dt/(x1[i]*dx2[j])*(  Bx1e(i,j+1,k) 
                                            - Bx1e(i,j,k)) ,
                      );
    #endif

    #elif ((PHYSICS == GRMHD) && (GR_METRIC != CARTESIAN)) 
    double dtdA3 = dt/grid->A[KDIR](i,j,k);
    // double dL1 = grid->L[IDIR];
    // double dL2 = grid->L[JDIR];
    double dL1m = grid->L[IDIR](i,j,k);
    double dL1p = grid->L[IDIR](i,j+1,k);
    double dL2m = grid->L[JDIR](i,j,k);
    double dL2p = grid->L[JDIR](i+1,j,k);

    // curlE = DIM_EXPAND(+ dtdA3*(emf->Ex2e(i+1,j,k)*dL2(i+1,j,k) - emf->Ex2e(i,j,k)*dL2(i,j,k)) ,
    //                    - dtdA3*(emf->Ex1e(i,j+1,k)*dL1(i,j+1,k) - emf->Ex1e(i,j,k)*dL1(i,j,k)) ,
    curlE = DIM_EXPAND(+ dtdA3*(Ex2e(i+1,j,k)*dL2p - Ex2e(i,j,k)*dL2m) ,
                       - dtdA3*(Ex1e(i,j+1,k)*dL1p - Ex1e(i,j,k)*dL1m) ,
                         0.0                                                                  );                                                                      ,

    #endif

    CT_ASSIGN(Rs[BX3s](i,j,k) = -curlE;,
              Rs[EX3s](i,j,k) = +curlB;)
  }}}
#endif

  // ------------------------------------------------------
  // 4. Add q*v source term to Ampere's law (only with ResRMHD)
  // ------------------------------------------------------

#if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  int ni, nj;
  int rec = RECONSTRUCTION;
  double qv;
  static Ary3D q_D;
  RBox box_rec;
  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);
  
  if (q_D.isAllocated() == 0) {
    q_D.create(nx1_tot,nx2_tot,nx3_tot);
    #pragma acc enter data copyin(q_D)
    #pragma omp target enter data map(q_D[:nx1_tot][0:nx2_tot][0:nx3_tot])
  }

  #pragma omp target teams distribute parallel for collapse(3)  
  #pragma acc parallel loop collapse(3) present(data, q_D)
  BOX_LOOP(tot_box, i,j,k){
    double vx = data->Vc(i,j,k,VX1);
    double vy = data->Vc(i,j,k,VX2);
    double vz = data->Vc(i,j,k,VX3);
    double gamma = 1.0/sqrt(1.0 - vx*vx - vy*vy - vz*vz);
    double D  = data->Vc(i,j,k,RHO)*gamma;
    q_D(i,j,k) = data->q(i,j,k)/D;
  }

  // 4a. Update x1-face staggered Electric field

  #if DIMENSIONS >= 1

  ARRAY_STRIDES (grid, IDIR, ni, nj);

  ArrayReconstruct<IDIR> (data, q_D, IDIR, rec, grid, &box_rec, 0); 

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) \
          present(data, q_D, Rs) 
  for (k = kbeg; k <= kend;   k++){
  for (j = jbeg; j <= jend;   j++){
  for (i = ibeg; i <= iend+1; i++){
  
    long int offset  = ni*(j + nj*k);
    Ary1D q_DL(&data->sweep->vL[offset],ni);
    Ary1D q_DR(&data->sweep->vR[offset],ni);

    double Frho = data->emf->Frho[IDIR](i,j,k);
    if (Frho > 0.0)  qv = q_DL(i)*Frho;
    else             qv = q_DR(i)*Frho;

    Rs[EX1s](i,j,k) -= dt*qv;
  }}}
  #endif

  // 4b. Update x2-face staggered Electric field

  #if DIMENSIONS >= 2 && INCLUDE_JDIR == TRUE
  ARRAY_STRIDES (grid, JDIR, ni, nj);

  ArrayReconstruct<JDIR> (data, q_D, JDIR, rec, grid, &box_rec, 0);
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) \
          present(data, q_D, Rs) 
  for (k = kbeg; k <= kend; k++)  {
  for (i = ibeg; i <= iend; i++)  {
  for (j = jbeg; j <= jend+1; j++){
  
    long int offset = ni*(i + nj*k);
    Ary1D q_DL(&data->sweep->vL[offset],ni);
    Ary1D q_DR(&data->sweep->vR[offset],ni);

    double Frho = data->emf->Frho[JDIR](i,j,k);
    if (Frho > 0.0)  qv = q_DL(j)*Frho;
    else             qv = q_DR(j)*Frho;
    Rs[EX2s](i,j,k) -= dt*qv;
  }}}
  #endif

  // 4c. Update x3-face staggered Electric field

  #if DIMENSIONS == 3 && INCLUDE_KDIR == TRUE

  ARRAY_STRIDES (grid, KDIR, ni, nj);

  ArrayReconstruct<KDIR> (data, q_D, KDIR, rec, grid, &box_rec, 0);

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) \
          present(data, q_D, Rs) 
  for (j = jbeg; j <= jend;   j++){
  for (i = ibeg; i <= iend;   i++){
  for (k = kbeg; k <= kend+1; k++){
    
    long int offset = ni*(i + nj*j);
    Ary1D q_DL(&data->sweep->vL[offset],ni);
    Ary1D q_DR(&data->sweep->vR[offset],ni);

    double Frho = data->emf->Frho[KDIR](i,j,k);
    if (Frho > 0.0)  qv = q_DL(k)*Frho;
    else             qv = q_DR(k)*Frho;
    Rs[EX3s](i,j,k) -= dt*qv;
  }}}
  #endif

#endif


mynvtxstop_();
}

////////////////////////////////////////////////////////////////////////  
/*!
  Check the divergence-free condition of magnetic field in the
  constrained transport formalism.
  The solenoidal condition is discretized in a finite-volume sense:
  \f[
       \int \nabla\cdot\vec{B}\, d^3x
     = \int \vec{B}\cdot d\hat{\bf S}
     = \sum_d \Big(\bar{B}_{d,+}S_{d,+} - \bar{B}_{d,-}S_{d,-} \Big)
  \f]
  where \f$S_{d,\pm}\f$ denotes the right (+) or left (-) interface surface
  areas orthogonal to the \c d direction.
  Thus in Cartesian coordinates one has
  \f[
       \int \nabla\cdot\vec{B}\, d^3x
     =  \Delta y\Delta z \left(B_{x,+} - B_{x,-}\right)
      + \Delta x\Delta z \left(B_{y,+} - B_{y,-}\right)
      + \Delta x\Delta y \left(B_{z,+} - B_{z,-}\right)
  \f]
  while in spherical coordinates (\f$ dS_1 = r^2\,\sin\theta\,d\theta\,d\phi,
  \, dS_2 = r\sin\theta\, dr\,d\phi,\, dS_3 = r\,dr\,d\theta\f$) we have
  \f[
       \int \nabla\cdot\vec{B}\, d^3x
     =  \Delta \mu\Delta\phi
          \Big(r^2_+B_{r,+} - r^2_-B_{r,-}\Big)
      + \Delta\left(\frac{r^2}{2}\right) \Delta\phi
          \Big(  \sin\theta_+B_{\theta,+} - \sin\theta_-B_{\theta,-}\Big)
      + \Delta\left(\frac{r^2}{2}\right)\Delta\theta
          \Big(B_{\phi,+} - B_{\phi,-}\Big)
  \f]
  where \f$\mu = 1-\cos\theta\f$.
  Notice also that \f$\Delta (r^2/2) = r\Delta r\f$ where \c r is the cell
  center.
 
  \param [in]  Bs    an array of staggered magnetic field components
  \param [in]  grid  a pointer to Grid structure
  \param [in]  box   the box over which div(B)=0 has to be checked
*/
//////////////////////////////////////////////////////////////////////// 
void CT_CheckDivB (const Ary3D (&Bs)[6], Grid *grid, RBox *box)
{
  int i,j,k;
  #if GEOMETRY == CYLINDRICAL
  int cyl_axis = (grid->lbound[IDIR] == SINGULAR_AXIS)*INCLUDE_IDIR;
  #endif
  double divB;

  double dBx1, dBx2, dBx3;
  double dbmax=0.0;
  DIM_EXPAND(auto &Ax1 = grid->A[IDIR];,
             auto &Ax2 = grid->A[JDIR];,
             auto &Ax3 = grid->A[KDIR];)
  DIM_EXPAND(auto &Bx1 = Bs[IDIR];,
             auto &Bx2 = Bs[JDIR];,
             auto &Bx3 = Bs[KDIR];)

  // ------------------------------------------------------
  // Loop over computational domain
  // ------------------------------------------------------

  BOX_LOOP(*box, i,j,k){
    DIM_EXPAND(dBx1 = Ax1(i+1,j,k)*Bx1(i+1,j,k) - Ax1(i,j,k)*Bx1(i,j,k);  ,
               dBx2 = Ax2(i,j+1,k)*Bx2(i,j+1,k) - Ax2(i,j,k)*Bx2(i,j,k);  ,
               dBx3 = Ax3(i,j,k+1)*Bx3(i,j,k+1) - Ax3(i,j,k)*Bx3(i,j,k);)

    // NOTE (cylindrical r-phi origin): 
    // Ghost zones across r=0 represent a \pi-rotated continuation of the solution,
    // not physical control volumes. To allow a single scalar div(B) diagnostic
    // to detect inconsistencies/NaNs in ghost staggered fields, the azimuthal
    // contribution is flipped so that the result matches the FV divergence of
    // the \pi-rotated physical cell.  This is a diagnostic convention only and
    // does NOT represent a geometric area or normal reversal.

    #if GEOMETRY == CYLINDRICAL
    if (cyl_axis && i < grid->lbeg[IDIR]) dBx2 *= -1.0;
    #endif

    divB = DIM_EXPAND(dBx1, + dBx2, + dBx3);

    dbmax = MAX(dbmax,fabs(divB));

    if (fabs(divB) > 1.e-6) {
      printLog ("! CT_CheckDivB(): div(B) = %12.6e, ijk = %d %d %d, stage = %d\n",
              divB, i,j,k, g_intStage);
      printLog ("                  dBx = %12.6e, dBy = %12.6e\n", dBx1, dBx2);
      
      printf ("-------------- Bs ---------------\n");
      DIM_EXPAND(printLog ("  Bx1: %12.6e  %12.6e\n",Bx1(i,j,k),Bx1(i+1,j,k)); ,
                 printLog ("  Bx2: %12.6e  %12.6e\n",Bx2(i,j,k),Bx2(i,j+1,k)); ,
                 printLog ("  Bx3: %12.6e  %12.6e\n",Bx3(i,j,k),Bx3(i,j,k+1)); )


      QuitPLUTO("from CT_CheckDivB()",1);
    }
  }
}

