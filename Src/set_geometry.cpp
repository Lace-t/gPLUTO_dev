////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Initialize geometry-dependent grid quantities.

  Compute grid quantities (such as interface areas, volumes,
  centroid of volume, etc..) that depend on the geometry.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Aug 20, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  \param [in,out] grid  Pointer to an array of Grid structures;
*/
////////////////////////////////////////////////////////////////////////
void SetGeometry (Grid *grid)
{
  int i, j, k, idim;
  int iend = grid->lend[IDIR] + grid->nghost[IDIR];
  int jend = grid->lend[JDIR] + grid->nghost[JDIR];
  int kend = grid->lend[KDIR] + grid->nghost[KDIR];
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  DIM_EXPAND(double *x1  = grid->xC[IDIR]; double *dx1 = grid->dx[IDIR];,             
             double *dx2 = grid->dx[JDIR];,
             double *dx3 = grid->dx[KDIR];)
  #if GEOMETRY != CARTESIAN
  double *x1p = grid->xR[IDIR], *x1m = grid->xL[IDIR];             
  #if GEOMETRY == SPHERICAL
  double *x2p = grid->xR[JDIR], *x2m = grid->xL[JDIR];
  #endif
  double dVr;
  #endif

  // 0. Allocate memory for the Grid structure.
  //    Values are defined either at cell centers or
  //    interfaces.

  grid->dV.create(nx1_tot,nx2_tot,nx3_tot);

  grid->A[IDIR].create((nx1_tot+1),nx2_tot,nx3_tot);
  grid->A[JDIR].create(nx1_tot,(nx2_tot+1),nx3_tot);
  grid->A[KDIR].create(nx1_tot,nx2_tot,(nx3_tot+1));

  grid->rt     = ARRAY_1D(nx1_tot, double);  // USELESS ?
  grid->rC     = grid->xC[IDIR];
  grid->rL     = grid->xL[IDIR];

  grid->dmu    = ARRAY_1D(nx2_tot, double);
  grid->dsn    = ARRAY_1D(nx2_tot, double);
  grid->cot    = ARRAY_1D(nx2_tot, double);
  grid->sC     = ARRAY_1D(nx2_tot, double);
  grid->sL     = ARRAY_1D(nx2_tot, double);

  for (idim = 0; idim < 3; idim++) {
    grid->inv_dx[idim]  = ARRAY_1D(grid->ntot[idim], double);
    grid->inv_dxi[idim] = ARRAY_1D(grid->ntot[idim], double);
  }

  // 1a. Compute positions arrays in the X1 (IDIR) direction

  for (i = 0; i <= iend; i++){
   
    #if GEOMETRY == CARTESIAN
    grid->rt[i]  = x1[i];
    #elif GEOMETRY == CYLINDRICAL
    grid->rt[i]  = x1[i];
    #elif GEOMETRY == SPHERICAL
    double xL = x1m[i];
    double xR = x1p[i];
    grid->rt[i]  = (xR*xR*xR - xL*xL*xL)/(xR*xR - xL*xL)/1.5;
    #endif
  }

  // 1b. Compute non-primary position arrays in the X2 
  //     (JDIR) direction. At present this is only done for 
  //     SPHERICAL geometry.
  //     Note that inverse (s_inv) is set to zero insetad of 
  //     infinity (ok as long as the quantity multiplying 
  //     s_inv goes to zero more rapidly.)

  #if GEOMETRY == SPHERICAL
  double *x2  = grid->xC[JDIR];
  for (j = 0; j <= jend; j++){
    double xL = x2m[j];
    double xR = x2p[j];
    double xC = x2[j];
    grid->dmu[j] = fabs(cos(xL) - cos(xR));
    grid->dsn[j] = sin(xR) - sin(xL);
    grid->cot[j] = grid->dsn[j]/grid->dmu[j];
    grid->sL[j]  = fabs(sin(xL));
    grid->sC[j]  = fabs(sin(xC));
  }
  #endif

  // 2. Compute volumes

  for (k = 0; k <= kend; k++){
  for (j = 0; j <= jend; j++){
  for (i = 0; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    grid->dV(i,j,k)  = DIM_EXPAND(dx1[i], *dx2[j], *dx3[k]);  // = dx*dy*dz 
    #elif GEOMETRY == CYLINDRICAL
    dVr = fabs(x1[i])*dx1[i];
    grid->dV(i,j,k)  = DIM_EXPAND(dVr, *dx2[j], *dx3[k]);     // = |r|*dr*dphi*dz 
    #elif GEOMETRY == SPHERICAL
    dVr = fabs(x1p[i]*x1p[i]*x1p[i] - x1m[i]*x1m[i]*x1m[i])/3.0;
    double dmu = fabs(cos(x2m[j]) - cos(x2p[j]));
    grid->dV(i,j,k)  = DIM_EXPAND(dVr, *dmu, *dx3[k]);        // = dVr*dmu*dphi
    #endif
  }}}

  // 3a. Compute area in the x1-direction

  for (k =  0; k <= kend;   k++){
  for (j =  0; j <= jend;   j++){
  for (i =  0; i <= iend+1; i++){
    #if GEOMETRY == CARTESIAN
    grid->A[IDIR](i,j,k) = DIM_EXPAND(1.0, *dx2[j], *dx3[k]);         // = dy*dz 
    #elif GEOMETRY == CYLINDRICAL
    if (i == iend+1) {
      grid->A[IDIR](i,j,k) = DIM_EXPAND(fabs(x1p[iend]), *dx2[j], *dx3[k]); // = rp*dphi*dz 
    }else{
      grid->A[IDIR](i,j,k) = DIM_EXPAND(fabs(x1m[i])   , *dx2[j], *dx3[k]); // = rp*dphi*dz 
    }
    #elif GEOMETRY == SPHERICAL
    double dmu = fabs(cos(x2m[j]) - cos(x2p[j]));
    if (i == iend+1) {
      grid->A[IDIR](i,j,k) = DIM_EXPAND(x1p[iend]*x1p[iend], *dmu, *dx3[k]); // = rp^2*dmu*dphi 
    }else{
      grid->A[IDIR](i,j,k) = DIM_EXPAND(x1m[i]*x1m[i]      , *dmu, *dx3[k]); // = rp^2*dmu*dphi 
    }
    #endif
  }}}

  // 3b. Compute area in the x2-direction

  for (k =  0; k <= kend;   k++){
  for (j =  0; j <= jend+1; j++){
  for (i =  0; i <= iend;   i++){
    #if GEOMETRY == CARTESIAN
    grid->A[JDIR](i,j,k) = DIM_EXPAND(dx1[i], *1.0, *dx3[k]);        // = dx*dz 
    #elif GEOMETRY == CYLINDRICAL
    grid->A[JDIR](i,j,k) = DIM_EXPAND(dx1[i], *1.0, *dx3[k]);        // = dr*dz 
    #elif GEOMETRY == SPHERICAL
    if (j == jend+1){
      grid->A[JDIR](i,j,k) = DIM_EXPAND(fabs(x1[i])*dx1[i], *fabs(sin(x2p[jend])), *dx3[k]); // = r*dr*sin(thp)*dphi 
    }else{
      grid->A[JDIR](i,j,k) = DIM_EXPAND(fabs(x1[i])*dx1[i], *fabs(sin(x2m[j]))   , *dx3[k]); // = r*dr*sin(thp)*dphi 
    }
    #endif
  }}}

  // 3c. Compute area in the x3-direction

  for (k =  0; k <= kend+1; k++){
  for (j =  0; j <= jend;   j++){
  for (i =  0; i <= iend;   i++){
    #if GEOMETRY == CARTESIAN
    grid->A[KDIR](i,j,k) = DIM_EXPAND(dx1[i], *dx2[j], *1.0);  // = dx*dy 
    #elif GEOMETRY == CYLINDRICAL
    grid->A[KDIR](i,j,k) = DIM_EXPAND(fabs(x1[i])*dx1[i], *dx2[j], *1.0);   // = |r|*dr*dphi 
    #elif GEOMETRY == SPHERICAL
    grid->A[KDIR](i,j,k) = DIM_EXPAND(fabs(x1[i])*dx1[i], *dx2[j], *1.0);   // = r*dr*dth 
    #endif
  }}}

  // 4. Compute and store the reciprocal of cell spacing
  //    between interfaces (inv_dx) and cell centers (inv_dxi)

  DIM_LOOP(idim){
    for (i = 0; i < grid->ntot[idim]; i++) {
      grid->inv_dx[idim][i] = 1.0/(grid->dx[idim][i]);
    }
  
    for (i = 0; i < grid->ntot[idim]-1; i++) {
      grid->inv_dxi[idim][i] = 2.0/(grid->dx[idim][i] + grid->dx[idim][i+1]);
    }
  }
}
