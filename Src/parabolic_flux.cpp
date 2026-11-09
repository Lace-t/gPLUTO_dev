////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Compute parabolic fluxes

  \author D. Crocco \n
          A. Mignone (andrea.mignone@unito.it)
  \date   Apr 19, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void ParabolicFlux (Data *data, timeStep *Dts, Grid *grid, RBox *box);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void ParabolicFlux_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) ParabolicFlux<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) ParabolicFlux<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) ParabolicFlux<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void ParabolicFlux (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  #pragma acc update host(data, grid,box)  // Ensure that the host has the 
                                           // latest data before starting the computation 
 
  mynvtxstart_("PARABOLIC",RED);
  int i, j, k, ti, tj, tk, nv;
  int ni,nj;
  int nbeg, nend, tbeg, tend, bbeg, bend;  

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  RBoxPermuteIndices<dir>(box, nbeg, nend, tbeg, tend, bbeg, bend);
 
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
 
  // ------------------------------------------------------
  // 0. Set location for current and resistivity arrays
  // ------------------------------------------------------

  data->J[IDIR].location_ = CELL_CENTER;
  data->J[JDIR].location_ = CELL_CENTER;
  data->J[KDIR].location_ = CELL_CENTER;
  data->eta[IDIR].location_ = CELL_CENTER;
  data->eta[JDIR].location_ = CELL_CENTER;
  data->eta[KDIR].location_ = CELL_CENTER;

  #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)  
    data->J[IDIR].location_ = X1_EDGE; 
    data->J[JDIR].location_ = X2_EDGE;
    data->J[KDIR].location_ = X3_EDGE;
    data->eta[IDIR].location_ = X1_EDGE;
    data->eta[JDIR].location_ = X2_EDGE;
    data->eta[KDIR].location_ = X3_EDGE;
  #endif

  #pragma acc update device(data, Dts, grid)  // Ensure that the device has the 
                                              // latest data before starting the computation 

  // ------------------------------------------------------
  // 1. Pre-compute resistivity and current density if needed
  // ------------------------------------------------------

  #if RESISTIVITY == EXPLICIT

  RBox box_b, box_j;

  RBoxDefine (0, nx1_tot-INCLUDE_IDIR, 
              0, nx2_tot-INCLUDE_JDIR, 
              0, nx3_tot-INCLUDE_KDIR, &box_b);

  RBoxDefine (1, nx1_tot-2, 
              1, nx2_tot-2, 
              1, nx3_tot-2, &box_j);
  
  #pragma acc enter data copyin (box_b, box_j)
    
  Ary3D B[3];
  
  B[0].create(nx1_tot,nx2_tot,nx3_tot);
  B[1].create(nx1_tot,nx2_tot,nx3_tot);
  B[2].create(nx1_tot,nx2_tot,nx3_tot);

  #pragma omp target enter data map(to: B[:3])
  #pragma acc enter data copyin(B[:3])

  #pragma acc parallel loop collapse(3) present(data,data->Vc)
  BOX_LOOP(box_b,i,j,k){
    B[0](i,j,k) = data->Vc(i,j,k,BX1);
    B[1](i,j,k) = data->Vc(i,j,k,BX2);
    B[2](i,j,k) = data->Vc(i,j,k,BX3);
  }

  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
    #if INCLUDE_IDIR
      #pragma acc parallel loop collapse(3) present (data,box_b)
      BOX_LOOP(box_b,i,j,k){B[0](i,j,k) = data->Vs[BX1s](i,j,k);}
    #endif
    #if INCLUDE_JDIR
      #pragma acc parallel loop collapse(3) present (data,box_b)
      BOX_LOOP(box_b,i,j,k){B[1](i,j,k) = data->Vs[BX2s](i,j,k);}
    #endif
    #if INCLUDE_KDIR
      #pragma acc parallel loop collapse(3) present (data,box_b)
      BOX_LOOP(box_b,i,j,k){B[2](i,j,k) = data->Vs[BX3s](i,j,k);}
    #endif
  #endif  // DIVB_CONTROL == CONSTRAINED_TRANSPORT
  
  VectorCurl(B, data->J, grid, &box_j);  //compute current density J = curl(B)  
  #pragma acc exit data delete(B[:3])  // Free the device memory for B 
    
  // the computestaggeredEta is computed on host because the indexes are interleaved a lot
  // can be improved to be computed on device if needed, but for now it is not a bottleneck
  #pragma acc update host (data,data->eta)
  ComputeStaggeredEta(data, grid);
  #pragma acc update device (data,data->eta)


  //    #pragma omp target exit data map(delete: box_b, box_j)
   #pragma acc exit data delete(box_b, box_j) 

  #endif  // RESISTIVITY == EXPLICIT

  // --------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  // in the direction given by g_dir.
  // --------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(data, Dts, grid)
  for (k = bbeg; k <= bend; k++){
  for (j = tbeg; j <= tend; j++){
  for (i = nbeg; i <= nend+1; i++){

    if constexpr(dir == IDIR)	{ti = i; tj = j; tk = k;}
    if constexpr(dir == JDIR)	{ti = j; tj = i; tk = k;}
    if constexpr(dir == KDIR)	{ti = j; tj = k; tk = i;}

    // ------------------------------------------------------
    // 2a. Compute offsets and assign local pointers
    // ------------------------------------------------------

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*offset;
    long int offset2 = dir*offset;

    Ary1D Dmax(&Dts->Dmax[offset2],ni);
    Ary2D vc(&data->sweep->vC[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    // ------------------------------------------------------
    // 2b. Loop over normal direction
    // ------------------------------------------------------

    double Jp[NVAR], etas[NVAR], vp[NVAR];

    for (nv = 0; nv < NVAR; nv++){
      vp[nv] = 0.5*(vc(i,nv) + vc(i-1,nv));
    }

    #if (DIVB_CONTROL != CONSTRAINED_TRANSPORT)
    if (dir == IDIR) {
      etas[DIRt] = 0.5*(data->eta[DIRt](ti,tj,tk)+data->eta[DIRt](ti-1,tj,tk));
      etas[DIRb] = 0.5*(data->eta[DIRb](ti,tj,tk)+data->eta[DIRb](ti-1,tj,tk));
      Jp[DIRt] = 0.5*(data->J[DIRt](ti,tj,tk)+data->J[DIRt](ti-1,tj,tk));
      Jp[DIRb] = 0.5*(data->J[DIRb](ti,tj,tk)+data->J[DIRb](ti-1,tj,tk));
    } if (dir == JDIR) {
      etas[DIRt] = 0.5*(data->eta[DIRt](ti,tj,tk)+data->eta[DIRt](ti,tj-1,tk));
      etas[DIRb] = 0.5*(data->eta[DIRb](ti,tj,tk)+data->eta[DIRb](ti,tj-1,tk));
      Jp[DIRt] = 0.5*(data->J[DIRt](ti,tj,tk)+data->J[DIRt](ti,tj-1,tk));
      Jp[DIRb] = 0.5*(data->J[DIRb](ti,tj,tk)+data->J[DIRb](ti,tj-1,tk));
    } if (dir == KDIR) {
      etas[DIRt] = 0.5*(data->eta[DIRt](ti,tj,tk)+data->eta[DIRt](ti,tj,tk-1));
      etas[DIRb] = 0.5*(data->eta[DIRb](ti,tj,tk)+data->eta[DIRb](ti,tj,tk-1));
      Jp[DIRt] = 0.5*(data->J[DIRt](ti,tj,tk)+data->J[DIRt](ti,tj,tk-1));
      Jp[DIRb] = 0.5*(data->J[DIRb](ti,tj,tk)+data->J[DIRb](ti,tj,tk-1));       
    }
    #endif

    #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
    if (dir == IDIR) {
      etas[DIRt] = 0.5*(data->eta[DIRt](ti,tj,tk)+data->eta[DIRt](ti,tj,tk+INCLUDE_KDIR));
      etas[DIRb] = 0.5*(data->eta[DIRb](ti,tj,tk)+data->eta[DIRb](ti,tj+INCLUDE_JDIR,tk));
      Jp[DIRt] = 0.5*(data->J[DIRt](ti,tj,tk)+data->J[DIRt](ti,tj,tk+INCLUDE_KDIR));
      Jp[DIRb] = 0.5*(data->J[DIRb](ti,tj,tk)+data->J[DIRb](ti,tj+INCLUDE_JDIR,tk));
    } if (dir == JDIR) {
      etas[DIRt] = 0.5*(data->eta[DIRt](ti,tj,tk)+data->eta[DIRt](ti+INCLUDE_IDIR,tj,tk));
      etas[DIRb] = 0.5*(data->eta[DIRb](ti,tj,tk)+data->eta[DIRb](ti,tj,tk+INCLUDE_KDIR));
      Jp[DIRt] = 0.5*(data->J[DIRt](ti,tj,tk)+data->J[DIRt](ti+INCLUDE_IDIR,tj,tk));
      Jp[DIRb] = 0.5*(data->J[DIRb](ti,tj,tk)+data->J[DIRb](ti,tj,tk+INCLUDE_KDIR));
    } if (dir == KDIR) {
      etas[DIRt] = 0.5*(data->eta[DIRt](ti,tj,tk)+data->eta[DIRt](ti,tj+INCLUDE_JDIR,tk));
      etas[DIRb] = 0.5*(data->eta[DIRb](ti,tj,tk)+data->eta[DIRb](ti+INCLUDE_IDIR,tj,tk));
      Jp[DIRt] = 0.5*(data->J[DIRt](ti,tj,tk)+data->J[DIRt](ti,tj+INCLUDE_JDIR,tk));
      Jp[DIRb] = 0.5*(data->J[DIRb](ti,tj,tk)+data->J[DIRb](ti+INCLUDE_IDIR,tj,tk));
    }
    #endif

    Dmax(i) = etas[DIRt] + etas[DIRb];

    flux(i,BXn) +=   0.0; //fluxes computation
    flux(i,BXt) += - etas[DIRb]*Jp[DIRb];
    flux(i,BXb) +=   etas[DIRt]*Jp[DIRt];

    #if HAVE_ENERGY
    flux(i,ENG) +=  vp[BXt]*(- etas[DIRb]*Jp[DIRb]) 
                  + vp[BXb]*(  etas[DIRt]*Jp[DIRt]);
    #endif

  }}} // End loop on (n,t,b)

  #pragma acc update host (data, grid)
  mynvtxstop_();
} 
