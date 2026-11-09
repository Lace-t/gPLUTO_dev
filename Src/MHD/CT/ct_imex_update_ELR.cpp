////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Update staggered electric field in resistive RMHD

  \b References

  \author A. Mignone  (andrea.mignone@unito.it) \n
          M. Rossazza (marco.rossazza@unito.it)     
  \date   Sep 29, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if (PHYSICS == ResRMHD)
////////////////////////////////////////////////////////////////////////
/*!
  Update face electric fields with the stiff part of the current,
  \f[
  E_f -= dt*j_f
  \f]
  where \f$ j_f \f$ is computed as arithmetic average of the
  cell-centered current obtained during the IMEX source step.
*/
////////////////////////////////////////////////////////////////////////
void CT_IMEXImplicitUpdate(Data *data, Ary4D &J, double dt1, Grid *grid)
{
  int    i,j,k;
  int    ni, nj;
  int    rec = RECONSTRUCTION;
  double gamma;

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  int nx1_tot=grid->ntot[IDIR];
  int nx2_tot=grid->ntot[JDIR];
  int nx3_tot=grid->ntot[KDIR];
  long int size_tot = nx3_tot*nx2_tot*nx1_tot;

  double inv_dt1 = 1.0/dt1;

  static Ary3D Eu, ux1, ux2, ux3, Bx1, Bx2, Bx3;
  RBox tot_box;
  RBox box_rec;
  State state;
  RunConfig &run_config = *(data->run_config);

  // ------------------------------------------------------
  // 0a. Allocate memory for local arrays
  // ------------------------------------------------------

  if (Eu.isAllocated() == 0) {
    Eu.create(nx1_tot,nx2_tot,nx3_tot);
    ux1.create(nx1_tot,nx2_tot,nx3_tot);
    ux2.create(nx1_tot,nx2_tot,nx3_tot);
    ux3.create(nx1_tot,nx2_tot,nx3_tot);
    Bx1.create(nx1_tot,nx2_tot,nx3_tot);
    Bx2.create(nx1_tot,nx2_tot,nx3_tot);
    Bx3.create(nx1_tot,nx2_tot,nx3_tot);

    #pragma acc enter data copyin(Eu[:nx1_tot*nx2_tot*nx3_tot], \
       ux1[:nx1_tot*nx2_tot*nx3_tot], ux2[:nx1_tot*nx2_tot*nx3_tot], ux3[:nx1_tot*nx2_tot*nx3_tot], \
       Bx1[:nx1_tot*nx2_tot*nx3_tot], Bx2[:nx1_tot*nx2_tot*nx3_tot], Bx3[:nx1_tot*nx2_tot*nx3_tot])

  }

  RBoxDefine (0, grid->ntot[IDIR]-1,
              0, grid->ntot[JDIR]-1,
              0, grid->ntot[KDIR]-1, &tot_box);

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box_rec);

  #pragma acc enter data copyin(tot_box, box_rec)
  // #pragma omp target enter data map(tot_box[:1], box_rec[:1])

  // ------------------------------------------------------
  // 0b. Compute 3D arrays to be interpolated (4-vel, mag.
  //     field, E.U)
  // ------------------------------------------------------

mynvtxstart_("IMEXImplicitUpdate",YELLOW);
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, tot_box) \
                                        present(Eu, ux1, ux2, ux3, Bx1, Bx2, Bx3)
  BOX_LOOP(tot_box, i,j,k){

    #ifdef HIGH_ORDER
    gamma =   data->Vpc(i,j,k,VX1)*data->Vpc(i,j,k,VX1)
            + data->Vpc(i,j,k,VX2)*data->Vpc(i,j,k,VX2)
            + data->Vpc(i,j,k,VX3)*data->Vpc(i,j,k,VX3);
    gamma = 1.0/sqrt(1.0 - gamma);

    ux1(i,j,k) = data->Vpc(i,j,k,VX1)*gamma;
    ux2(i,j,k) = data->Vpc(i,j,k,VX2)*gamma;
    ux3(i,j,k) = data->Vpc(i,j,k,VX3)*gamma;

    Eu(i,j,k) =   data->Vpc(i,j,k,EX1)*ux1(i,j,k)
                + data->Vpc(i,j,k,EX2)*ux2(i,j,k)
                + data->Vpc(i,j,k,EX3)*ux3(i,j,k);

    Bx1(i,j,k) = data->Vpc(i,j,k,BX1);
    Bx2(i,j,k) = data->Vpc(i,j,k,BX2);
    Bx3(i,j,k) = data->Vpc(i,j,k,BX3);
    #else
    gamma =   data->Vc(i,j,k,VX1)*data->Vc(i,j,k,VX1)
            + data->Vc(i,j,k,VX2)*data->Vc(i,j,k,VX2)
            + data->Vc(i,j,k,VX3)*data->Vc(i,j,k,VX3);
    gamma = 1.0/sqrt(1.0 - gamma);

    ux1(i,j,k) = data->Vc(i,j,k,VX1)*gamma;
    ux2(i,j,k) = data->Vc(i,j,k,VX2)*gamma;
    ux3(i,j,k) = data->Vc(i,j,k,VX3)*gamma;

    Eu(i,j,k) =   data->Vc(i,j,k,EX1)*ux1(i,j,k)
                + data->Vc(i,j,k,EX2)*ux2(i,j,k)
                + data->Vc(i,j,k,EX3)*ux3(i,j,k);

    Bx1(i,j,k) = data->Vc(i,j,k,BX1);
    Bx2(i,j,k) = data->Vc(i,j,k,BX2);
    Bx3(i,j,k) = data->Vc(i,j,k,BX3);
    #endif
  }

  // ------------------------------------------------------
  // 1. Update x1-face electric field and compute
  //    interface current from
  //
  //         Ex(new) = Ex(old) - dt1*Jx
  //     --> Jx = [ Ex(old) - Ex(new) ]/dt1
  // ------------------------------------------------------

  #if INCLUDE_IDIR
  ARRAY_STRIDES (grid, IDIR, ni, nj);

  ArrayReconstruct<IDIR> (data, Eu,  IDIR, rec, grid, &box_rec, 0); 
  ArrayReconstruct<IDIR> (data, ux1, IDIR, rec, grid, &box_rec, 1); 
  ArrayReconstruct<IDIR> (data, ux2, IDIR, rec, grid, &box_rec, 2); 
  ArrayReconstruct<IDIR> (data, ux3, IDIR, rec, grid, &box_rec, 3); 
  ArrayReconstruct<IDIR> (data, Bx2, IDIR, rec, grid, &box_rec, 4); 
  ArrayReconstruct<IDIR> (data, Bx3, IDIR, rec, grid, &box_rec, 5); 

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3)  present(data, J) 
  for (k = kbeg; k <= kend;   k++){
  for (j = jbeg; j <= jend;   j++){
  for (i = ibeg; i <= iend+1; i++){

    long int offset  = ni*(j + nj*k);
    long int offset1 = ni*(j + nj*k) + size_tot;
    long int offset2 = ni*(j + nj*k) + size_tot*2;
    long int offset3 = ni*(j + nj*k) + size_tot*3;
    long int offset4 = ni*(j + nj*k) + size_tot*4;
    long int offset5 = ni*(j + nj*k) + size_tot*5;
    Ary1D EuL(&data->sweep->vL[offset],ni);
    Ary1D EuR(&data->sweep->vR[offset],ni);
    Ary1D ux1L(&data->sweep->vL[offset1],ni);
    Ary1D ux1R(&data->sweep->vR[offset1],ni);
    Ary1D ux2L(&data->sweep->vL[offset2],ni);
    Ary1D ux2R(&data->sweep->vR[offset2],ni);
    Ary1D ux3L(&data->sweep->vL[offset3],ni);
    Ary1D ux3R(&data->sweep->vR[offset3],ni);
    Ary1D Bx2L(&data->sweep->vL[offset4],ni);
    Ary1D Bx2R(&data->sweep->vR[offset4],ni);
    Ary1D Bx3L(&data->sweep->vL[offset5],ni);
    Ary1D Bx3R(&data->sweep->vR[offset5],ni);

    double eta = Resistive_eta(state, run_config);
    double const uxL = ux1L(i), uyL = ux2L(i), uzL = ux3L(i);
    double const uxR = ux1R(i), uyR = ux2R(i), uzR = ux3R(i);

    double const ByL = Bx2L(i), BzL = Bx3L(i);
    double const ByR = Bx2R(i), BzR = Bx3R(i);

    double const Ex  = data->Vs[EX1s](i,j,k);

    double const gammaL = sqrt(1.0 + uxL*uxL + uyL*uyL + uzL*uzL);
    double const gammaR = sqrt(1.0 + uxR*uxR + uyR*uyR + uzR*uzR);

    double const uXBL = uyL*BzL - uzL*ByL;
    double const uXBR = uyR*BzR - uzR*ByR;

#if 0  // IN CASE WE NEED A VARIABLE-DEPENDENT RESISTIVITY [TODO...]
    double v[256];
    v[RHO] = state->Vc(i,j,k,RHO);    
    v[VX1] = 0.5*(ux1L(i)/gammaL + ux1R(i)/gammaR);
    v[VX2] = 0.5*(ux1L(i)/gammaL + ux1R(i)/gammaR);
    v[VX3] = 0.5*(ux1L(i)/gammaL + ux1R(i)/gammaR);

    STATE_SET(state, vl, ul, bgf, x1L[i], x2C[j], x3C[k]);
    double eta = Resistive_eta(state, run_config);
#endif

    double const ExL = (eta*Ex - dt1*(uXBL - EuL(i)*uxL/gammaL))/(eta + dt1*gammaL);
    double const ExR = (eta*Ex - dt1*(uXBR - EuR(i)*uxR/gammaR))/(eta + dt1*gammaR);

    #ifdef HIGH_ORDER
    data->Vps[EX1s](i,j,k) = 0.5*(ExL + ExR);
    #else
    data->Vs [EX1s](i,j,k) = 0.5*(ExL + ExR);
    #endif

    J(i,j,k,IDIR) = (Ex - data->Vs[EX1s](i,j,k))*inv_dt1;
  }}}
  #endif

  // ------------------------------------------------------
  // 2. Update x2-face electric field
  // ------------------------------------------------------

  #if INCLUDE_JDIR
  ARRAY_STRIDES (grid, JDIR, ni, nj);

  ArrayReconstruct<JDIR> (data, Eu,  JDIR, rec, grid, &box_rec, 0); 
  ArrayReconstruct<JDIR> (data, ux1, JDIR, rec, grid, &box_rec, 1); 
  ArrayReconstruct<JDIR> (data, ux2, JDIR, rec, grid, &box_rec, 2); 
  ArrayReconstruct<JDIR> (data, ux3, JDIR, rec, grid, &box_rec, 3); 
  ArrayReconstruct<JDIR> (data, Bx1, JDIR, rec, grid, &box_rec, 4); 
  ArrayReconstruct<JDIR> (data, Bx3, JDIR, rec, grid, &box_rec, 5); 

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3)  present(data, J) 
  for (k = kbeg; k <= kend; k++){
  for (i = ibeg; i <= iend; i++){
  for (j = jbeg; j <= jend+1; j++){

    long int offset = ni*(i + nj*k);
    long int offset1 = ni*(i + nj*k) + size_tot;
    long int offset2 = ni*(i + nj*k) + size_tot*2;
    long int offset3 = ni*(i + nj*k) + size_tot*3;
    long int offset4 = ni*(i + nj*k) + size_tot*4;
    long int offset5 = ni*(i + nj*k) + size_tot*5;
    Ary1D EuL(&data->sweep->vL[offset],ni);
    Ary1D EuR(&data->sweep->vR[offset],ni);
    Ary1D ux1L(&data->sweep->vL[offset1],ni);
    Ary1D ux1R(&data->sweep->vR[offset1],ni);
    Ary1D ux2L(&data->sweep->vL[offset2],ni);
    Ary1D ux2R(&data->sweep->vR[offset2],ni);
    Ary1D ux3L(&data->sweep->vL[offset3],ni);
    Ary1D ux3R(&data->sweep->vR[offset3],ni);
    Ary1D Bx1L(&data->sweep->vL[offset4],ni);
    Ary1D Bx1R(&data->sweep->vR[offset4],ni);
    Ary1D Bx3L(&data->sweep->vL[offset5],ni);
    Ary1D Bx3R(&data->sweep->vR[offset5],ni);
 
    double eta = Resistive_eta(state, run_config);
    double const uxL = ux1L(j), uyL = ux2L(j), uzL = ux3L(j);
    double const uxR = ux1R(j), uyR = ux2R(j), uzR = ux3R(j);

    double const BxL = Bx1L(j), BzL = Bx3L(j);
    double const BxR = Bx1R(j), BzR = Bx3R(j);

    double const Ey  = data->Vs[EX2s](i,j,k);

    double const gammaL = sqrt(1.0 + uxL*uxL + uyL*uyL + uzL*uzL);
    double const gammaR = sqrt(1.0 + uxR*uxR + uyR*uyR + uzR*uzR);

    double const uXBL = uzL*BxL - uxL*BzL;
    double const uXBR = uzR*BxR - uxR*BzR;

    double const EyL = (eta*Ey - dt1*(uXBL - EuL(j)*uyL/gammaL))/(eta + dt1*gammaL);
    double const EyR = (eta*Ey - dt1*(uXBR - EuR(j)*uyR/gammaR))/(eta + dt1*gammaR);

    #ifdef HIGH_ORDER
    data->Vps[EX2s](i,j,k) = 0.5*(EyL + EyR);
    #else
    data->Vs [EX2s](i,j,k) = 0.5*(EyL + EyR);
    #endif

    J(i,j,k,JDIR) = (Ey - data->Vs[EX2s](i,j,k))*inv_dt1;
  }}}
  #endif

  // ------------------------------------------------------
  // 3. Update x3-face electric field
  // ------------------------------------------------------ 
  
  #if INCLUDE_KDIR
  ARRAY_STRIDES (grid, KDIR, ni, nj);

  ArrayReconstruct<KDIR> (data, Eu,  KDIR, rec, grid, &box_rec, 0); 
  ArrayReconstruct<KDIR> (data, ux1, KDIR, rec, grid, &box_rec, 1); 
  ArrayReconstruct<KDIR> (data, ux2, KDIR, rec, grid, &box_rec, 2); 
  ArrayReconstruct<KDIR> (data, ux3, KDIR, rec, grid, &box_rec, 3); 
  ArrayReconstruct<KDIR> (data, Bx1, KDIR, rec, grid, &box_rec, 4); 
  ArrayReconstruct<KDIR> (data, Bx2, KDIR, rec, grid, &box_rec, 5); 

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3)  present(data,J)
  for (i = ibeg; i <= iend; i++){
  for (j = jbeg; j <= jend; j++){
  for (k = kbeg; k <= kend+1; k++){
    
    long int offset = ni*(i + nj*j);
    long int offset1 = ni*(i + nj*j) + size_tot;
    long int offset2 = ni*(i + nj*j) + size_tot*2;
    long int offset3 = ni*(i + nj*j) + size_tot*3;
    long int offset4 = ni*(i + nj*j) + size_tot*4;
    long int offset5 = ni*(i + nj*j) + size_tot*5;
    Ary1D EuL(&data->sweep->vL[offset],ni);
    Ary1D EuR(&data->sweep->vR[offset],ni);
    Ary1D ux1L(&data->sweep->vL[offset1],ni);
    Ary1D ux1R(&data->sweep->vR[offset1],ni);
    Ary1D ux2L(&data->sweep->vL[offset2],ni);
    Ary1D ux2R(&data->sweep->vR[offset2],ni);
    Ary1D ux3L(&data->sweep->vL[offset3],ni);
    Ary1D ux3R(&data->sweep->vR[offset3],ni);
    Ary1D Bx1L(&data->sweep->vL[offset4],ni);
    Ary1D Bx1R(&data->sweep->vR[offset4],ni);
    Ary1D Bx2L(&data->sweep->vL[offset5],ni);
    Ary1D Bx2R(&data->sweep->vR[offset5],ni);
    
    double eta = Resistive_eta(state, run_config);

    double const uxL = ux1L(k), uyL = ux2L(k), uzL = ux3L(k);
    double const uxR = ux1R(k), uyR = ux2R(k), uzR = ux3R(k);

    double const BxL = Bx1L(k), ByL = Bx2L(k);
    double const BxR = Bx1R(k), ByR = Bx2R(k);

    double const Ez  = data->Vs[EX3s](i,j,k);

    double const gammaL = sqrt(1.0 + uxL*uxL + uyL*uyL + uzL*uzL);
    double const gammaR = sqrt(1.0 + uxR*uxR + uyR*uyR + uzR*uzR);

    double const uXBL = uxL*ByL - uyL*BxL;
    double const uXBR = uxR*ByR - uyR*BxR;

    double const EzL = (eta*Ez - dt1*(uXBL - EuL(k)*uzL/gammaL))/(eta + dt1*gammaL);
    double const EzR = (eta*Ez - dt1*(uXBR - EuR(k)*uzR/gammaR))/(eta + dt1*gammaR);

    #ifdef HIGH_ORDER
    data->Vps[EX3s](i,j,k) = 0.5*(EzL + EzR);
    #else
    data->Vs [EX3s](i,j,k) = 0.5*(EzL + EzR);
    #endif

    J(i,j,k,KDIR) = (Ez - data->Vs[EX3s](i,j,k))*inv_dt1;
  }}}
  #endif

  // ------------------------------------------------------
  // 4. Average face values to cell center
  // ------------------------------------------------------

  #ifndef HIGH_ORDER
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for(k = kbeg; k <= kend; k++){
  for(j = jbeg; j <= jend; j++){
  for(i = ibeg; i <= iend; i++){
    DIM_EXPAND(
      data->Uc(i,j,k,EX1) = 0.5*(data->Vs[EX1s](i,j,k) + data->Vs[EX1s](i+1,j,k));  ,
      data->Uc(i,j,k,EX2) = 0.5*(data->Vs[EX2s](i,j,k) + data->Vs[EX2s](i,j+1,k));  ,
      data->Uc(i,j,k,EX3) = 0.5*(data->Vs[EX3s](i,j,k) + data->Vs[EX3s](i,j,k+1));
    )
  }}}
  #endif

  #pragma acc exit data delete(tot_box, box_rec)
  #pragma omp target exit data map(delete:tot_box[:1], box_rec[:1])
mynvtxstop_();
}

////////////////////////////////////////////////////////////////////////
void CT_ComputeCharge(Data *data, RBox *box, Grid *grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int    i,j,k;
  double dx, dy, dz;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  // ------------------------------------------------------
  // 1. Compute charge
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
    DIM_EXPAND(dx = grid->dx[IDIR][i];,
               dy = grid->dx[JDIR][j];,
               dz = grid->dx[KDIR][k];)
    data->q(i,j,k) = DIM_EXPAND((data->Vs[EX1s](i+1,j,k) - data->Vs[EX1s](i,j,k))/dx,
                              + (data->Vs[EX2s](i,j+1,k) - data->Vs[EX2s](i,j,k))/dy,
                              + (data->Vs[EX3s](i,j,k+1) - data->Vs[EX3s](i,j,k))/dz );
  }}}
}
#endif // (PHYSICS == ResRMHD)  

