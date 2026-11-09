////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief   Add HO source terms to the conservative right hand side

  \authors A. Mignone (andrea.mignone@unito.it)\n
           V. Berta   (vittoria.berta@unito.it)
           
  \date    Nov 4, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if PHYSICS == MHD
  #if BACKGROUND_FIELD == YES
    #define TotBBi(v, b0, a, b) (v(i,a)*(v(i,b) + b0[b-BX1]) + b0[a-BX1]*v(i,b))
  #else
    #define TotBBi(v, b0, a, b) (v(i,a)*v(i,b))
  #endif
#else
  #define TotBBi(v, b0, a, b)  0.0
#endif

#if PHYSICS == MHD
  #if BACKGROUND_FIELD == YES
    #define TotBBj(v, b0, a, b) (v(j,a)*(v(j,b) + b0[b-BX1]) + b0[a-BX1]*v(j,b))
  #else
    #define TotBBj(v, b0, a, b) (v(j,a)*v(j,b))
  #endif
#else
  #define TotBBj(v, b0, a, b)  0.0
#endif

#if PHYSICS == MHD
  #if BACKGROUND_FIELD == YES
    #define TotBBk(v, b0, a, b) (v(k,a)*(v(k,b) + b0[b-BX1]) + b0[a-BX1]*v(k,b)])
  #else
    #define TotBBk(v, b0, a, b) (v(k,a)*v(k,b))
  #endif
#else
  #define TotBBk(v, b0, a, b)  0.0
#endif

////////////////////////////////////////////////////////////////////////
void Sources (Data *data, Ary4D &dUc, double dt, Grid *grid, RBox *box)
/*!
  Compute high-order quadrature of source terms.
  
  - Gravity / Body forces
  - Geometrical source terms
  - Powell / EGLM source terms
  - Shearing Box 
  - FARGO
  - DUST fluid source terms
*/
////////////////////////////////////////////////////////////////////////
 {
  int i, j, k, nv;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  #if GEOMETRY == SPHERICAL
  double *rt   = grid->rt;
  double *sp   = grid->sp;
  double *sm   = grid->sp-1;
  double *s    = grid->s;
  double *dmu  = grid->dmu;
  #endif

  // --------------------------------------------------------
  // 0. Defining integration box & allocate arrays
  // --------------------------------------------------------
  
#if (GEOMETRY != CARTESIAN) || (BODY_FORCE & VECTOR) \
    || (PHYSICS == ResRMHD) && (DIVE_CONTROL == NO)

  #pragma acc data present(data)
  #ifdef HIGH_ORDER
  const Ary4D& Vc = data->Vpc; 
  #else
  const Ary4D& Vc = data->Vc;
  #endif
  
  int ns = 1;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  int bibeg, biend;
  int bjbeg, bjend;
  int bkbeg, bkend;
  RBox   boxs;
  RBoxDefine(ibeg - ns*INCLUDE_IDIR, iend + ns*INCLUDE_IDIR,
             jbeg - ns*INCLUDE_JDIR, jend + ns*INCLUDE_JDIR,
             kbeg - ns*INCLUDE_KDIR, kend + ns*INCLUDE_KDIR, &boxs);

  bibeg = boxs.ibeg; biend = boxs.iend;
  bjbeg = boxs.jbeg; bjend = boxs.jend;
  bkbeg = boxs.kbeg; bkend = boxs.kend;
#endif

  static Ary4D Sr;

  if (Sr.isAllocated() == 0) {
    Sr.create(nx1_tot,nx2_tot,nx3_tot,NVAR);
    #pragma acc enter data copyin(Sr)
  }

  mynvtxstart_("Source",BLUE);

  // -------------------------------------------------------- 
  // 1. Compute pointwise source terms
  // --------------------------------------------------------

  #if (GEOMETRY != CARTESIAN) || (BODY_FORCE & VECTOR)
  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];
  #pragma acc parallel loop collapse(2) present(data, Sr, x1)
  for (k = bkbeg; k <= bkend; k++){
  for (j = bjbeg; j <= bjend; j++){
  #pragma acc loop vector
  for (i = bibeg; i <= biend; i++){

    double g[3] = {0.0, 0.0, 0.0};
    double v[NVAR];

    #pragma acc loop
    NVAR_LOOP(nv) {
      v[nv] = Vc(i, j, k, nv); // Point-wise value 
      Sr(i, j, k, nv) = 0.0;
    }

    // --------------------------------------------
    // 1a. Add body forces
    // --------------------------------------------
      
    #if (BODY_FORCE & VECTOR) 
    #pragma acc routine (BodyForceVector) seq
    BodyForceVector (v, g, x1[i], x2[j], x3[k]);
    #endif

    // --------------------------------------------
    // 1b. Add geometrical sources containing 1/R 
    // factors
    // --------------------------------------------

    #if EOS == IDEAL
    double pt = v[PRS];
    #elif EOS == ISOTHERMAL
    double cs = data->run_config->isoSoundSpeed;
    double pt = cs*cs*v[RHO];
    #endif

    #if (GEOMETRY == CYLINDRICAL)
    #if PHYSICS == MHD
    pt += 0.5*DOT_PRODUCT(v + BX1, v + BX1);
    #endif
    Sr(i,j,k,MX1) = v[RHO]*g[IDIR]*x1[i] + pt + v[RHO]*v[iVPHI]*v[iVPHI];
    Sr(i,j,k,MX2) = v[RHO]*g[JDIR]*x1[i]      - v[RHO]*v[iVR]*v[iVPHI];
    Sr(i,j,k,MX3) = v[RHO]*g[KDIR];
    #if PHYSICS == MHD
    Sr(i,j,k,MX1) -= v[iBPHI]*v[iBPHI];
    Sr(i,j,k,MX2) += v[iBPHI]*v[iBR];
    double Ez_R = -(v[iVR]*v[iBPHI] - v[iVPHI]*v[iBR]);
    Sr(i,j,k,BX2) -= Ez_R;
    #endif
    #endif
  
    #if (GEOMETRY == SPHERICAL)
    #if PHYSICS == MHD
    pt += 0.5*DOT_PRODUCT(v + BX1, v + BX1);
    #endif
    Sr(i,j,k,MX1) =   v[RHO]*g[IDIR]*x1[i] + 2.0*pt 
                      + v[RHO]*(v[VX2]*v[VX2] + v[VX3]*v[VX3]);
    Sr(i,j,k,MX2) = v[RHO]*g[JDIR]*x1[i] - v[RHO]*v[iVR]*v[iVPHI];
    Sr(i,j,k,MX3) = v[RHO]*g[KDIR];
    #if PHYSICS == MHD
    Sr(i,j,k,MX1) -= v[iBPHI]*v[iBPHI];
    Sr(i,j,k,MX2) += v[iBPHI]*v[iBR];
    double Ez_R = -(v[iVR]*v[iBPHI] - v[iVPHI]*v[iBR]);
    Sr(i,j,k,BX2) -= Ez_R;
    #endif
    #endif

    #if (BODY_FORCE & VECTOR) 
    #if HAVE_ENERGY
    Sr(i,j,k,ENG) = v[RHO]*DOT_PRODUCT(v+VX1, g)*x1[i];
    #endif
    #endif

  }}}
  #endif // (GEOMETRY != CARTESIAN) || (BODY_FORCE & VECTOR)

  // --------------------------------------------------------- 
  // 2a. Compute point-wise source term for ResRMHD: 
  // add charge term (q*v) to zone-centered Ampere's law 
  // ---------------------------------------------------------

#if (PHYSICS == ResRMHD) && (DIVE_CONTROL == NO)
  DIM_EXPAND(double *dx1 = grid->dx[IDIR];,
             double *dx2 = grid->dx[JDIR];,
             double *dx3 = grid->dx[KDIR];)
  #pragma acc parallel loop collapse(2) present(data, grid, Sr)
  for (k = bkbeg; k <= bkend; k++){
  for (j = bjbeg; j <= bjend; j++){
  #pragma acc loop vector
  for (i = bibeg; i <= biend; i++){
    double qx = 0.0, qy = 0.0, qz = 0.0;

    // -- 1st approach: simpler but requires 4 ghosts -- 
    
    DIM_EXPAND(qx = 8.0*(Vc(i+1,j,k,EX1)- Vc(i-1,j,k,EX1))  
                      - (Vc(i+2,j,k,EX1)- Vc(i-2,j,k,EX1));  
               qx /= 12.0*dx1[i];                                  ,

               qy  = 8.0*(Vc(i,j+1,k,EX2) - Vc(i,j-1,k,EX2))  
                       - (Vc(i,j+2,k,EX2) - Vc(i,j-2,k,EX2));
               qy /= 12.0*dx2[j];                                  ,

               qz  = 8.0*(Vc(i,j,k+1,EX3) - Vc(i,j,k-1,EX3))  
                       - (Vc(i,j,k+2,EX3) - Vc(i,j,k-2,EX3));
               qz /= 12.0*dx3[k]; )  

    Sr(i,j,k,EX1) = -(qx + qy + qz)*Vc(i,j,k,VX1);
    Sr(i,j,k,EX2) = -(qx + qy + qz)*Vc(i,j,k,VX2);
    Sr(i,j,k,EX3) = -(qx + qy + qz)*Vc(i,j,k,VX3);
  }}}
#endif // (PHYSICS == ResRMHD) && (DIVE_CONTROL == NO) 

  // --------------------------------------------------------- 
  // 2. Use Simpson rule to compute source terms
  // ---------------------------------------------------------

  #pragma acc parallel loop collapse(4) present(data, Sr, box)
  for (nv = 0; nv < NVAR; nv++){
  BOX_LOOP(*box,i,j,k){
      double S  =  Sr(i,j,k,nv);

      DIM_EXPAND(S += HO_Average4D<IDIR>(Sr, i,j,k,nv)/24.0;  ,
                 S += HO_Average4D<JDIR>(Sr, i,j,k,nv)/24.0;  ,
                 S += HO_Average4D<KDIR>(Sr, i,j,k,nv)/24.0;)

      // #if 1 // HO_WENO_LAPLACIAN == NO
      // DIM_EXPAND(
      //    S += HO_LAPLXc(Sr,i,j,k,nv)/24.0;   ,
      //    S += HO_LAPLYc(Sr,i,j,k,nv)/24.0;   ,
      //    S += HO_LAPLZc(Sr,i,j,k,nv)/24.0; 
      // )
      // #else
      // DIM_EXPAND(
      //    S += WENOLaplacianXc(Sr,i,j,k,nv)/24.0;   ,
      //    S += WENOLaplacianYc(Sr,i,j,k,nv)/24.0;   ,
      //    S += WENOLaplacianZc(Sr,i,j,k,nv)/24.0; 
      // )
      // #endif
      dUc(i,j,k,nv) += dt*S;
  }}

  mynvtxstop_();
}
