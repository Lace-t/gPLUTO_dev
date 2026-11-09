#include "gpluto.hpp"

template<int dir>
inline long int SweepIndex (const Grid *grid, int, int, int, int);

////////////////////////////////////////////////////////////////////////
/* Check whether symmetry is preserved across a symmetry axis         */
////////////////////////////////////////////////////////////////////////

void CheckFluxSymmetry(const Data *data, RBox *box, Grid *grid)
{
  int i,j,k, nv;
  double eps = 1.e-12;
  double fzm, fzp;

long int ni,nj;
ARRAY_STRIDES (grid, KDIR, ni, nj);

// Create Ary4D()

  static Ary4D Flux;
  if (Flux.isAllocated() == 0) {
    Flux.create(grid->ntot[IDIR], grid->ntot[JDIR], grid->ntot[KDIR],NVAR);
  }

// Copy flux to Ary4D()

  RBox box_tot;
  RBoxDefine(0, grid->ntot[IDIR]-1, 0, grid->ntot[JDIR]-1, 0, grid->ntot[KDIR]-1, 
             &box_tot);  

  for (nv = 0; nv < NVAR; nv++){
  BOX_LOOP(box_tot,i,j,k){
    long int offset1 = NVAR*ni*(i + nj*j);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Flux(i,j,k,nv) = flux(k,nv);
  }}

  // Check data equivalence (variables)

BOX_LOOP(*box,i,j,k){   
  NVAR_LOOP(nv){
  long int k1 = SweepIndex<KDIR>(grid, i, j, k, nv);

   double dv = data->sweep->vC[k1] - data->Vc(i,j,k,nv);
   if (fabs(dv) > 1.e-12){
     printLog ("! Variables are not the same at (%d %d %d)\n", i,j,k);
     QuitPLUTO("Stop in CheckFluxSymmetry\n",1);
   }
  }
}


  // Check flux symmetry at the axis

  // Are we at the axis ? 
  if (grid->lbound[IDIR] == AXISYMMETRIC){

    // Compare fluxes at (i,j,k) and (i+1,j,k)
    j = 0;
    
    for (k = box->kbeg; k <= box->kend; k++){
      int  im = grid->lbeg[IDIR]-1;
      int  ip = grid->lbeg[IDIR];

      NVAR_LOOP(nv){

        long int km = SweepIndex<KDIR>(grid, im, j, k, nv);
        long int kp = SweepIndex<KDIR>(grid, ip, j, k, nv);

        int sign = 1;
        if (nv == MX1) sign = -1;

        double fzm = data->sweep->flux[km];
        double fzp = data->sweep->flux[kp];

// fzm = Flux(im,j,k,nv);
// fzp = Flux(ip,j,k,nv);

        double dflx = fzp - fzm*sign;
        if (fabs(dflx) > eps){
          printLog ("! CheckFluxSymmetry(): flux not symmetric (i,j,k = %d, %d, %d, nv = %d)\n",
                     im,j,k, nv);
          printLog ("  Fz(i-1) = %12.6e; Fz(i) = %12.6e\n",fzm, fzp);
          printLog ("  Vc      = %12.6e, %12.6e\n",data->Vc(im,j,k,nv), data->Vc(ip,j,k,nv));
          QuitPLUTO("Stop in CheckFluxSymmetry\n",1);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
/* Map grid indices (i,j,k,nv) to the linear index of the sweep buffer
   for the given sweep direction `dir`.
   Same logic used in the ARRAY_STRIDES macro                         */
////////////////////////////////////////////////////////////////////////
template<int dir>
inline long int SweepIndex(const Grid *grid, int i, int j, int k, int nv)
{
  long int idx;

  if constexpr (dir == IDIR){

    const int ni = grid->ntot[IDIR];
    const int nj = grid->ntot[JDIR];
    idx = NVAR*(ni*(j + nj*k)) + nv*ni + i;

  } else if constexpr (dir == JDIR){

    const int ni = grid->ntot[JDIR];
    const int nj = grid->ntot[IDIR];
    idx = NVAR*(ni*(i + nj*k)) + nv*ni + j;

  } else { // KDIR

    const int ni = grid->ntot[KDIR];
    const int nj = grid->ntot[IDIR];
    idx = NVAR*(ni*(i + nj*j)) + nv*ni + k;
  }

  return idx;
}

