#include "gpluto.hpp"

/* *************************************************************** */
void UserDefOutput (const Data *d, Grid *grid)
/*
 *
 ***************************************************************** */
{
  int i, j, k;  
  int *ntot = grid->ntot;
  static Ary3D flag;
  RBox box;

  RBoxDefine (grid->lbeg, grid->lend, &box);

  if (flag.isAllocated() == NO){
    flag.create (ntot[IDIR], ntot[JDIR], ntot[KDIR]);
  }

  LinkUserVar(&flag, "flag", DBL_OUTPUT);

  BOX_LOOP(box, i,j,k){
    flag(i,j,k) = d->flag(i,j,k);
  }

}
