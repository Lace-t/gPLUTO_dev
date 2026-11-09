////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Add Powell source term to MHD equastions

  Add Powell's 8 wave source terms to the right hand side of MHD equations.
    
  \b Reference
    - "A conservative orbital advection scheme for simulations 
       of magnetized shear flows with the PLUTO code"\n
       Mignone et al, A&A (2012) 545, A152 (--> Appendix A.1.1)

  \author A. Mignone (andrea.mignone@unito.it)\n
          V. Berta   (vittoria.berta@unito.it)
  \date   Jan 18, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void RightHandSideSource (Data *data, Grid *grid, double dt, RBox *box)
/*!
*/
////////////////////////////////////////////////////////////////////////
 {
  int  i, j, k, nv;
  int  ni, nj;
  int  dir = g_dir;
  int  ibeg = box->ibeg, iend = box->iend;
  int  jbeg = box->jbeg, jend = box->jend;
  int  kbeg = box->kbeg, kend = box->kend;

  double *x1  = grid->xC[IDIR];
  double *x2  = grid->xC[JDIR];
  double *x3  = grid->xC[KDIR];
  double **Bg0;
  double src[NVAR];

  #if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
  int  BXn = BX1 + dir;
  #endif

  ARRAY_OFFSET (grid, ni, nj);

  NVAR_LOOP(nv) src[nv] = 0.0;

  if (g_dir == IDIR){

    #pragma acc parallel loop collapse(2) present(data, x1, x2, x3) \
                                          private (src[:NVAR])
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
      long int offset = NVAR*ni*(j + nj*k);
      Ary2D vc(&data->sweep.vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep.vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep.vR[offset],ni,NVAR);
      Ary2D rhs(&data->sweep.rhs[offset],ni,NVAR);
      Ary2D flux(&data->sweep.flux[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary1D Bn(&data->sweep.Bn[offset],ni);
      #endif

      #pragma acc loop vector
      for (i = ibeg; i <= iend; i++){

        double inv_dx = grid->inv_dx[IDIR][i];

        #if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
        double Bp, Bm, divB;

        Bp = 0.5*(vL(i+1,BXn) + vR(i+1,BXn));
        Bm = 0.5*(vL(i,BXn)   + vR(i,BXn));
        #if GEOMETRY == CARTESIAN
        divB = (Bp - Bm)*inv_dx;
        #else
Bp = flux(i+1,RHO) > 0.0 ? vL(i+1,BX1):vR(i+1,BX1);
Bm = flux(i,RHO)   > 0.0 ? vL(i,BX1)  :vR(i,BX1);
divB = (grid->A[IDIR](i+1,j,k)*Bp - grid->A[IDIR](i,j,k)*Bm)/grid->dV(i,j,k);
        #endif

        double factor = divB*dt;

        src[MX1] = factor*vc(i,BX1);
        src[MX2] = factor*vc(i,BX2);
        src[MX3] = factor*vc(i,BX3);
        src[BX1] = factor*vc(i,VX1);
        src[BX2] = factor*vc(i,VX2);
        src[BX3] = factor*vc(i,VX3);
        #if HAVE_ENERGY
        src[ENG] = factor*(  vc(i,VX1)*vc(i,BX1)
                           + vc(i,VX2)*vc(i,BX2)
                           + vc(i,VX3)*vc(i,BX3));
        #endif
 
        OrthoToCons(src, x1[i], x2[j], x3[k]);

        #pragma acc loop seq 
        NVAR_LOOP(nv) rhs(i,nv) -= src[nv];      

        #endif  // (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
    }}}

  }else if (g_dir == JDIR){

    #pragma acc parallel loop collapse(2) present(data, x1, x2, x3)  \
                                          private (src[:NVAR])
    for (k = kbeg; k <= kend; k++){
    for (i = ibeg; i <= iend; i++){

      long int offset = NVAR*ni*(i + nj*k);
      Ary2D vc(&data->sweep.vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep.vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep.vR[offset],ni,NVAR);
      Ary2D rhs(&data->sweep.rhs[offset],ni,NVAR);
      Ary2D flux(&data->sweep.flux[offset],ni,NVAR);

      #if GEOMETRY != CARTESIAN
      Ary1D Bn(&data->sweep.Bn[offset],ni);
      #endif

      #pragma acc loop vector
      for (j = jbeg; j <= jend; j++){

        double inv_dx = grid->inv_dx[JDIR][j];

        #if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
        double Bp, Bm, divB;

        Bp = 0.5*(vL(j+1,BXn)   + vR(j+1,BXn));
        Bm = 0.5*(vL(j,BXn) + vR(j,BXn));
        #if GEOMETRY == CARTESIAN
        divB = (Bp - Bm)*inv_dx;
        #else
        divB = (grid->A[JDIR](i,j+1,k)*Bn(j+1) - grid->A[JDIR](i,j,k)*Bn(j))/grid->dV(i,j,k);
Bp = flux(j+1,RHO) > 0.0 ? vL(j+1,BX2):vR(j+1,BX2);
Bm = flux(j,RHO)   > 0.0 ? vL(j,BX2)  :vR(j,BX2);
divB = (grid->A[JDIR](i,j+1,k)*Bp - grid->A[JDIR](i,j,k)*Bm)/grid->dV(i,j,k);

        #endif
        
        double factor = divB*dt;
        src[MX1] = factor*vc(j,BX1);
        src[MX2] = factor*vc(j,BX2);
        src[MX3] = factor*vc(j,BX3);
        src[BX1] = factor*vc(j,VX1);
        src[BX2] = factor*vc(j,VX2);
        src[BX3] = factor*vc(j,VX3);
        #if HAVE_ENERGY
        src[ENG] = factor*(  vc(j,VX1)*vc(j,BX1)
                           + vc(j,VX2)*vc(j,BX2)
                           + vc(j,VX3)*vc(j,BX3));
        #endif

        OrthoToCons(src, x1[i], x2[j], x3[k]);

        #pragma acc loop seq 
        NVAR_LOOP(nv) rhs(j,nv) -= src[nv];      

        #endif // (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)

    }}}

  }else if (g_dir == KDIR){

    #pragma acc parallel loop collapse(2) present(data, x1, x2, x3) \
                                          private(src[:NVAR])
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      long int offset = NVAR*ni*(i + nj*j);
      Ary2D vc(&data->sweep.vC[offset],ni,NVAR);
      Ary2D vL(&data->sweep.vL[offset],ni,NVAR);
      Ary2D vR(&data->sweep.vR[offset],ni,NVAR);
      Ary2D rhs(&data->sweep.rhs[offset],ni,NVAR);
      #if GEOMETRY != CARTESIAN
      Ary1D Bn(&data->sweep.Bn[offset],ni);
      #endif

      #pragma acc loop vector
      for (k = kbeg; k <= kend; k++){

        double inv_dx = grid->inv_dx[KDIR][k];
        double scrh  = dt;

        #if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
        double Bp, Bm, divB;
        Bp = 0.5*(vL(k+1,BXn) + vR(k+1,BXn));
        Bm = 0.5*(vL(k,BXn)   + vR(k,BXn));
        #if GEOMETRY == CARTESIAN
        divB = (Bp - Bm)*inv_dx;
        #else
        divB = (grid->A[KDIR](i,j,k+1)*Bn(k+1) - grid->A[KDIR](i,j,k)*Bn(k))/grid->dV(i,j,k);
        #endif

        double factor = dt*divB;
        src[MX1] = factor*vc(k,BX1);
        src[MX2] = factor*vc(k,BX2);
        src[MX3] = factor*vc(k,BX3);
        src[BX1] = factor*vc(k,VX1);
        src[BX2] = factor*vc(k,VX2);
        src[BX3] = factor*vc(k,VX3);
        #if HAVE_ENERGY
        src[ENG] = factor*(  vc(k,VX1)*vc(k,BX1)
                           + vc(k,VX2)*vc(k,BX2)
                           + vc(k,VX3)*vc(k,BX3));
        #endif

        OrthoToCons(src, x1[i], x2[j], x3[k]);

        #pragma acc loop seq 
        NVAR_LOOP(nv) rhs(k,nv) -= src[nv];      

        #endif
    }}}
  }
}
