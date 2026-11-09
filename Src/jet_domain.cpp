////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Adjust the effective domain size when the -xNjet option is given.

  The SetJetDomain() function shortens the domain integration index in the 
  direction specified by either '\c -x1jet', '\c -x2jet' or '\c -x3jet' to save 
  computational time.
  This is done by setting the final integration index in the direction
  of jet propagation to that of the first zone (counting from the top) with 
  non-zero pressure gradient (GetRightmostIndex()).
  The function UnsetJetDomain() restores the initial computational domain
  size.
  Useful for problems involving jet propagation.

  \note In parallel, the domain is \e not decomposed along the
        propagation direction (see ParseCmdLineArgs()).
  
  \author A. Mignone (andrea.mignone@unito.it)
  \date   Sep 14, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

struct Jet{
  int jd_nbeg;
  int jd_nend;
  int jd_npt;
  int jd_ntot;
  int rbound;
};

static Jet jet;
#pragma acc declare create(jet)

int GetRightmostIndex(int, Data const *, Grid*, int var);

////////////////////////////////////////////////////////////////////////
/*!
  
  Adjust the size of computational domain by reducing the final
  computationa index in the direction of propagation.
  
  \param [in]         d   pointer to Data structure
  \param [in]       dir   the direction of propagation
  \param [in]  log_freq   the output log frequency
  \param [in,out]  grid   pointer to Grid structure
*/
////////////////////////////////////////////////////////////////////////
void SetJetDomain (Data const *data, int dir, int log_freq, Grid *grid)
{
  mynvtxstart_("SetJetDomain",GREEN);
  int ngh;
  int n;
  static int first_call = 1;

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  #if HAVE_ENERGY
  constexpr int var = PRS;
  #else
  constexpr int var = RHO;
  #endif

  int  nx1_tot = grid->ntot[IDIR];
  int  nx2_tot = grid->ntot[JDIR];
  int  nx3_tot = grid->ntot[KDIR];
  int  nx1_int = grid->nint[IDIR];
  int  nx2_int = grid->nint[JDIR];
  int  nx3_int = grid->nint[KDIR];

  ngh = grid->nghost[dir];

  // 0. Save original domain offsets 
  //    and return to call for the first time

  if (first_call){
    if (dir == IDIR) {
      jet.jd_nbeg = ibeg; jet.jd_nend = iend; 
      jet.jd_npt  = nx1_int;  jet.jd_ntot = nx1_tot;
    }else if (dir == JDIR){
      jet.jd_nbeg = jbeg; jet.jd_nend = jend; 
      jet.jd_npt  = nx2_int;  jet.jd_ntot = nx2_tot;
    }else if (dir == KDIR){
      jet.jd_nbeg = kbeg; jet.jd_nend = kend; 
      jet.jd_npt  = nx3_int;  jet.jd_ntot = nx3_tot;
    }

    jet.rbound = grid->rbound[dir];

    #pragma acc update device(jet)
    first_call = 0;
    mynvtxstop_();

    return;
  }

  // 1. Find where grad(p) is not zero and 
  //  add safety guard cells 

  n = GetRightmostIndex(dir, data, grid, var) + 2*ngh;

  if (n < jet.jd_nbeg + ngh) n = jet.jd_nbeg + ngh;
  if (n > jet.jd_nend)       n = jet.jd_nend;

  #ifdef PARALLEL
  int n_glob;
  MPI_Allreduce (&n, &n_glob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  n = n_glob;
  #endif

  if (g_stepNumber%log_freq==0){
    print ("- SetJetDomain: dir: %d, index %d / %d\n", dir, n,jet.jd_nend);
  }

  // 2. Change global integer variables giving the rightmost integration 
  //    indexes. Also, suppress the rightmost physics boundary (rbound) 
  //    if solution has not yet reached it. 

  if (dir == IDIR) {
    iend    = grid->lend[IDIR] = n;
    nx1_tot = grid->ntot[IDIR] = iend + ngh;
    nx1_int = grid->nint[IDIR] = iend - ngh + 1;

    grid->rbound[IDIR] = jet.rbound;
    if (jet.jd_nend != iend) grid->rbound[IDIR] = 0;
    
  } else if (dir == JDIR){


    jend    = grid->lend[JDIR] = n;
    nx2_tot = grid->ntot[JDIR] = jend + ngh;
    nx2_int = grid->nint[JDIR] = jend - ngh + 1;

    grid->rbound[JDIR] = jet.rbound;
    if (jet.jd_nend != jend) grid->rbound[JDIR] = 0;

  }else if (dir == KDIR){

    kend    = grid->lend[KDIR] = n;
    nx3_tot = grid->ntot[KDIR] = kend + ngh;
    nx3_int = grid->nint[KDIR] = kend - ngh + 1;

    grid->rbound[KDIR] = jet.rbound;
    if (jet.jd_nend != kend) grid->rbound[KDIR] = 0;

  }
  #pragma acc update device(grid)
  mynvtxstop_();

}

////////////////////////////////////////////////////////////////////////
/*!
  Restore original (full domain) indexes.
*/
////////////////////////////////////////////////////////////////////////
void UnsetJetDomain (int dir, Grid *grid)
{
  mynvtxstart_("UnSetJetDomain",GREEN);

  grid->lbeg[dir] = jet.jd_nbeg;
  grid->lend[dir] = jet.jd_nend;
  grid->ntot[dir] = jet.jd_ntot;
  grid->nint[dir] = jet.jd_npt;
  grid->rbound[dir] = jet.rbound;

/*
  if (dir == IDIR){
    grid->lbeg[IDIR] = jet.jd_nbeg;
    grid->lend[IDIR] = jet.jd_nend;
    grid->ntot[IDIR] = jet.jd_ntot;
    grid->nint[IDIR] = jet.jd_npt;
    grid->rbound[IDIR] = jet.rbound;
  }

  if (dir == JDIR){
    grid->lbeg[JDIR] = jet.jd_nbeg;
    grid->lend[JDIR] = jet.jd_nend;
    grid->ntot[JDIR] = jet.jd_ntot;
    grid->nint[JDIR] = jet.jd_npt;
    grid->rbound[JDIR] = jet.rbound;
  }

  if (dir == KDIR){
    grid->lbeg[KDIR] = jet.jd_nbeg;
    grid->lend[KDIR] = jet.jd_nend;
    grid->ntot[KDIR] = jet.jd_ntot;
    grid->nint[KDIR] = jet.jd_npt;
    grid->rbound[KDIR] = jet.rbound;
  }
*/
  mynvtxstop_();
}

////////////////////////////////////////////////////////////////////////
/*!
  
  Find the local index j where grad(p) exceeds
  a cerating threshold
*/
////////////////////////////////////////////////////////////////////////
int GetRightmostIndex (int dir, Data const *data, Grid *grid, int var)
{

// -- backward sweep on x-axis -- 

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  int beg = -1;
  if (dir == IDIR){
    int iJetBeg = jet.jd_nbeg;
    #pragma acc parallel loop collapse(2) present(data, grid, jet) copy(iJetBeg) reduction(max:iJetBeg)
    for (int k = kbeg; k <= kend; k++){
    for (int j = jbeg; j <= jend; j++){
      int ii = jet.jd_nbeg;
      bool done = false;
      #pragma acc loop seq
      for (int i = jet.jd_nend; i >= jet.jd_nbeg; i--){
        if(done) continue;

        double dp = data->Vc(i+1,j,k,var) - data->Vc(i-1,j,k,var);
        dp /= data->Vc(i+1,j,k,var) + data->Vc(i-1,j,k,var);

        if (fabs(dp) > 1.e-5){
          done = true;
          ii = i;
        }
      }
      iJetBeg = std::max(ii, iJetBeg);
    
    }}
    beg = (iJetBeg);
  }

  if (dir == JDIR){

    int jJetBeg = jet.jd_nbeg;
    #pragma acc parallel loop collapse(2) present(data, grid, jet) copy(jJetBeg) reduction(max:jJetBeg)
    for (int k = kbeg; k <= kend; k++){
    for (int i = ibeg; i <= iend; i++){
      int jj = jet.jd_nbeg;
      bool done = false;

      #pragma acc loop seq
      for (int j = jet.jd_nend; j >= jet.jd_nbeg; j--){
        if(done) continue;

        double dp = data->Vc(i,j+1,k,var) - data->Vc(i,j-1,k,var);
        dp /= data->Vc(i,j+1,k,var) + data->Vc(i,j-1,k,var);
        if (fabs(dp) > 1.e-5){
          done = true;
          jj = j;
        }
      }
      jJetBeg = std::max(jj, jJetBeg);
    
    }}
    beg = (jJetBeg);
  }

  if (dir == KDIR){
    int kJetBeg = jet.jd_nbeg;
    #pragma acc parallel loop collapse(2) present(data, grid, jet) copy(kJetBeg) reduction(max:kJetBeg)
    for (int j = jbeg; j <= jend; j++){
    for (int i = ibeg; i <= iend; i++){
      int kk = jet.jd_nbeg;
      bool done = false;
      #pragma acc loop seq
      for (int k = jet.jd_nend; k >= jet.jd_nbeg; k--){
        if(done) continue;

        double dp = data->Vc(i,j,k+1,var) - data->Vc(i,j,k-1,var);
        dp /= data->Vc(i,j,k+1,var) + data->Vc(i,j,k-1,var);
        if (fabs(dp) > 1.e-5){

          done = true;
          kk = k;
        }
      }
      kJetBeg = std::max(kk, kJetBeg);
    
    }}
    beg = (kJetBeg);
  }

  return beg + static_cast<int>(0.05*static_cast<double>(jet.jd_ntot));
}
