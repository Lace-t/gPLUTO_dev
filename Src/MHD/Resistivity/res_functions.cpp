////////////////////////////////////////////////////////////////////////
/*!
  //to be commented by D.Crocco
  \file
  \brief  Compute resistivity.

  Compute resistivity using the same staggering of the current:

  -  eta_x1 at   (i, j+1/2, k+1/2)
  -  eta_x2 at   (i+1/2, j, k+1/2)
  -  eta_x3 at   (i+1/2, j+1/2, k)
 
  Since each eta_xn may depend on the total current, we also need
  to interpolate different current components at the same place.

  \c eta is stored in this function as a static array.

  \authors D. Crocco \n
           A. Mignone (andrea.mignone@unito.it)
           
  \date   Apr 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void ComputeStaggeredEta (const Data *data, Grid *grid)
/*!
  
  \param [in,out]  d       pointer to the PLUTO data structure
  \param [in]      grid    pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
  int i, j, k, nv;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  double Js[3], etas[3], vs[NVAR];
  double *x1, *x2, *x3;
  double *x1r, *x2r, *x3r;
  RBox box_red;
  RunConfig &run_config = *(data->run_config);

  x1 = grid->xC[IDIR]; x1r = grid->xR[IDIR];
  x2 = grid->xC[JDIR]; x2r = grid->xR[JDIR];
  x3 = grid->xC[KDIR]; x3r = grid->xR[KDIR];

  RBoxDefine(1, nx1_tot-1, 1, nx2_tot-1, 1, nx3_tot-1,&box_red);

  BOX_LOOP (box_red, i, j, k) {
    Js[IDIR] = data->J[IDIR](i,j,k);
    Js[JDIR] = AVERAGE_XY(data->J[JDIR], i-1, j, k);
    Js[KDIR] = AVERAGE_XZ(data->J[KDIR], i-1, j, k);
    NVAR_LOOP(nv){
      vs[nv] = 0.25*(data->Vc(i,j,k,nv)+data->Vc(i,j-1,k,nv)
                    +data->Vc(i,j,k-1,nv)+data->Vc(i,j-1,k-1,nv));
    }
    Resistive_eta (vs, x1[i], x2r[j], x3r[k], Js, run_config, etas);
    data->eta[IDIR](i,j,k) = etas[IDIR];

    Js[IDIR] = AVERAGE_XY(data->J[IDIR], i, j-1, k);
    Js[JDIR] = data->J[JDIR](i,j,k);
    Js[KDIR] = AVERAGE_YZ(data->J[KDIR], i, j-1, k);
    NVAR_LOOP(nv){
      vs[nv] = 0.25*(data->Vc(i,j,k,nv)+data->Vc(i-1,j,k,nv)
                    +data->Vc(i,j,k-1,nv)+data->Vc(i-1,j,k-1,nv));
    }
    Resistive_eta (vs, x1r[i], x2[j], x3r[k], Js, run_config, etas);
    data->eta[JDIR](i,j,k) = etas[JDIR];

    Js[IDIR] = AVERAGE_XZ(data->J[IDIR], i, j, k-1);
    Js[JDIR] = AVERAGE_YZ(data->J[JDIR], i, j, k-1);
    Js[KDIR] = data->J[KDIR](i,j,k);
    NVAR_LOOP(nv){
      vs[nv] = 0.25*(data->Vc(i,j,k,nv)+data->Vc(i-1,j,k,nv)
                    +data->Vc(i,j-1,k,nv)+data->Vc(i-1,j-1,k,nv));
    }
    Resistive_eta (vs, x1r[i], x2r[j], x3[k], Js, run_config, etas);
    data->eta[KDIR](i,j,k) = etas[KDIR];
  }
}