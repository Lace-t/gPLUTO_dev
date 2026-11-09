#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void UserDefOutput (const Data *, Grid *)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
/*
  int i,j,k;
  int *ntot = grid ->ntot;
  double *dx = grid->dx[IDIR]; // shortcut to dx 
  double *dy = grid->dx[JDIR]; // shortcut to dx 
  static Ary3D T, vortz;
  RBox box;

  RBoxDefine (grid->lbeg, grid->lend, &box); 
  if (T.isAllocated() == NO){
    T.create (ntot[IDIR], ntot[JDIR], ntot[KDIR]);
    vortz.create (ntot[IDIR], ntot[JDIR], ntot[KDIR]); 
  }

  LinkUserVar(&T, "T", DBL_OUTPUT); 
  LinkUserVar(&vortz, "vortz", DBL_OUTPUT);

  BOX_LOOP(box,i,j,k){
    T(i,j,k)     = data->Vc(i,j,k,PRS)/data->Vc(i,j,k,RHO);
    vortz(i,j,k) = 0.5*(data->Vc(i+1,j,k,VX2) - data->Vc(i-1,j,k,VX2))/dx[i]
  }

  SetOutputVar("prs", PARTICLES_FLT_OUTPUT, YES);
*/
}


