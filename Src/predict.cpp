#include "gpluto.hpp"

void Predict (Data *data, Grid *grid)
{
  int i,j,k;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR]; 

  static Ary3D rho0, rho1;
  static Ary3D prs0, prs1;
  
  RBox box, box_tot;

  // 0. Define boxes and allocate memory

  RBoxDefine (ibeg, iend, jbeg, jend, kbeg, kend, &box);
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &box_tot);
  #pragma acc enter data copyin(box, box_tot)

  if (rho0.isAllocated() == 0){
    rho0.create (nx1_tot,nx2_tot,nx3_tot);
    rho1.create (nx1_tot,nx2_tot,nx3_tot);
    prs0.create (nx1_tot,nx2_tot,nx3_tot);
    prs0.create (nx1_tot,nx2_tot,nx3_tot);
  }

  // 1. 
  static int first_call = 1;
  FILE *fp;
  if (first_call) fp = fopen("rho.dat","w");
  else            fp = fopen("rho.dat","a");

  i = 51; j = 0; k = 9;
  fprintf (fp, "%d   %12.6e   %12.6e  %12.6e  %12.6e  %12.6e\n", 
           g_stepNumber-49100,  data->Dts->time, 
           data->Vc(i,j,k,RHO), data->Vc(i,j,k,PRS), 
           data->Uc(i,j,k,RHO), data->Uc(i,j,k,ENG));

  fclose(fp);
  first_call = 0; 
}