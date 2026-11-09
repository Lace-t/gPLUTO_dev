////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Exchange boundary conditions for a flag() array, with buffers

  This function exchanges and merges potentially overlapping flags 
  in the ghost zones and using a logical or (|) operator.

  Here \c op can be one of:

  - <tt>FLAG_COPY</tt>: values in ghost zones are copied from adjacent 
                        processors active zones (default bounady exchange)
  - <tt>FLAG_OR</tt>:   values in ghost zones are merged with adjacent 
                        processors active zones using logical \c or operator
  - <tt>FLAG_AND</tt>:  values in ghost zones are merged with adjacent 
                        processors active zones using logical \c and operator

  In other words, this function does a symmetric merge because either 
  rank might already have flags set in its ghost zones after the 
  "enlarge-by-buffer" step, and you want the union (OR) across the interface.

  \author A. Mignone  (andrea.mignone@unito.it)
  \date   Apr 17, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define DEBUG NO

// int IsFlagged(int, int, const Grid *);

#ifdef PARALLEL

// Explicit template instantiations
template void ParallelExchangeFlag<FLAG_COPY>(Ary3Duint32_t &, const Grid *);
template void ParallelExchangeFlag<FLAG_OR>  (Ary3Duint32_t &, const Grid *);
template void ParallelExchangeFlag<FLAG_AND> (Ary3Duint32_t &, const Grid *);

////////////////////////////////////////////////////////////////////////
/*!
  \param [in,out] flag  an unit32_t array filled also in the ghost 
                        zones.
  \param [in]     grid  pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
template<int flag_op> 
void ParallelExchangeFlag (Ary3Duint32_t &flag, const Grid *grid)
{
  int  i,j,k;
  int  par_dim[3] = {0, 0, 0};

  const int  *ntot = grid->ntot;
  const int  *ngh  = grid->nghost;

  MPI_Status status;

  static int **nrnks;
  static int ***nrnks_edg;
  static int ***nrnks_crn;

  static uint32_t *send_bufL;
  static uint32_t *recv_bufL;
  static uint32_t *send_bufR;
  static uint32_t *recv_bufR;

  int nmax = ntot[IDIR]*ntot[JDIR];
  nmax     = MAX(nmax, ntot[JDIR]*ntot[KDIR]);
  nmax     = MAX(nmax, ntot[IDIR]*ntot[KDIR]);

  int nmaxg = ngh[IDIR];
  nmaxg     = MAX(nmaxg,ngh[JDIR]);
  nmaxg     = MAX(nmaxg,ngh[KDIR]);
  nmaxg    *= 2;

  int size  = nmaxg*nmax;

  // ------------------------------------------------------
  // 0a. Allocate arrays
  // ------------------------------------------------------

  if (nrnks == NULL){
    nrnks     = ARRAY_2D(3, 2, int);
    nrnks_edg = ARRAY_3D(3, 2, 2, int);
    nrnks_crn = ARRAY_3D(2, 2, 2, int);

    send_bufL = ARRAY_1D(size, uint32_t);  
    recv_bufL = ARRAY_1D(size, uint32_t);
    send_bufR = ARRAY_1D(size, uint32_t);
    recv_bufR = ARRAY_1D(size, uint32_t);

    #pragma acc enter data create(send_bufL[:size], recv_bufL[:size])
    #pragma acc enter data create(send_bufR[:size], recv_bufR[:size])

    GetNeighbourRanks(grid, nrnks, nrnks_edg, nrnks_crn);
  }

  // ------------------------------------------------------
  // 0b. Check the number of processors in each direction
  // ------------------------------------------------------

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

  // ------------------------------------------------------
  // 1. Loop over directions.
  //    Box naming convention:
  //           X         X         X                   
  //  ........ |---------|---------|---------|---------|......... 
  // |recv_boxL|send_boxL|                   |send_boxR|recv_boxR
  //
  // Note that recv_box'es are the same as the send_box'es
  // ------------------------------------------------------


// if (g_rank == 0){
// RBox boxt, box;
// RBoxDefine(0,ntot[IDIR]-1,0,ntot[JDIR]-1,0,ntot[KDIR]-1, &boxt);
// RBoxDefine(grid->lbeg, grid->lend, &box);
//   int count = 0;
//   BOX_LOOP(box, i,j,k){
//     if (IsFlagged(i,j,k,grid) && flag(i,j,k) != 0){
//       printLog ("Flag(%d, %d, %d) = %d\n",i,j,k,flag(i,j,k));
//       count++;
//     }
//   }
//   printLog ("Total count (before) = %d\n",count);
// }


  int dir;
  for (dir = 0; dir < 3; dir++){
    
    if (!par_dim[dir]) continue;  // Direction is not parallel, skip
    #if DEBUG == YES
    printLog ("----> dir = %d\n",dir);
    #endif

    int  l;
    int  beg[3], end[3];
    long ind;
    int  countL, countR;
    int  rankL = nrnks[dir][0];
    int  rankR = nrnks[dir][1];
    RBox send_boxL, send_boxR;
    RBox recv_boxL, recv_boxR; 

    // ----------------------------------------------------
    // 1A. Construct send recv boxes which coincide on each
    //     side of the domain and they are: 
    //
    //     i.   nghost long (when flag_op == FLAG_COPY) 
    //     ii. 2*nghost long (in FLAG_OR / FLAG_AND mode)
    // 
    //     in the exchange direction:
    //
    //      |...|...|---|---|    . . .      |---|---|...|...|
    //        0   1   2   3                  e-3 e-2 e-1  e
    //     
    //     Here e = end[dir]
    // ----------------------------------------------------

    for (l = 0; l < 3; l++) {  // Initialize to whole domain
      beg[l] = 0;
      end[l] = ntot[l] - 1;
    }

    // 1A1. Define send_boxL, recv_boxL, send_boxR and recv_boxR
    //      in FLAG_OR or FLAG_AND mode

    if constexpr (flag_op == FLAG_OR || flag_op == FLAG_AND){
      beg[dir] = 0;
      end[dir] = 2*ngh[dir] - 1;

      RBoxDefine (beg, end, &send_boxL);
      RBoxDefine (beg, end, &recv_boxL);
      
      beg[dir] = (ntot[dir] - 1) - (2*ngh[dir] - 1);
      end[dir] = (ntot[dir] - 1);

      RBoxDefine (beg, end, &send_boxR);
      RBoxDefine (beg, end, &recv_boxR); 
    }
    
    // 1A2. Define send_boxL, recv_boxL, send_boxR and recv_boxR
    //      in FLAG_COPY

    if constexpr (flag_op == FLAG_COPY) {

      beg[dir] = grid->lbeg[dir];
      end[dir] = beg[dir] + grid->nghost[dir] - 1;

      RBoxDefine (beg, end, &send_boxL);
      
      // Define recv_boxR 

      beg[dir] += grid->nint[dir];
      end[dir] += grid->nint[dir];

      RBoxDefine (beg, end, &recv_boxR);

      // Define send_boxR

      beg[dir] -= grid->nghost[dir];
      end[dir] -= grid->nghost[dir];

      RBoxDefine (beg, end, &send_boxR);

      // Define recv_boxL

      beg[dir] = 0;
      end[dir] = grid->nghost[dir]-1;

      RBoxDefine (beg, end, &recv_boxL);
    }

    countL = (send_boxL.nx1)*(send_boxL.nx2)*(send_boxL.nx3); 
    countR = (send_boxR.nx1)*(send_boxR.nx2)*(send_boxR.nx3); 

if (countL != countR){
  printLog("! Error in ParallelExchangeFlag: countL = %d, countR = %d\n", 
                 countL, countR);
  QuitPLUTO(" ",1);
}
    #pragma acc enter data copyin(send_boxL, send_boxR, recv_boxL, recv_boxR)

    // ----------------------------------------------------
    // 1B. Pack data (fill send buffers)
    // ----------------------------------------------------

    #pragma acc parallel loop collapse(3) \
                present(flag, send_bufL[:size], send_boxL)
    BOX_LOOP(send_boxL,i,j,k) {
      ind = BOX_INDEX(send_boxL,i,j,k);
      send_bufL[ind] = flag(i,j,k);
    }

    #pragma acc parallel loop collapse(3) \
                present(flag, send_bufR[:size], send_boxR)
    BOX_LOOP(send_boxR,i,j,k) {
      ind = BOX_INDEX(send_boxR,i,j,k);
      send_bufR[ind] = flag(i,j,k);
    }

    #if DEBUG == YES
    printLog ("rankL = %d, rankR = %d\n", rankL, rankR);
    printLog("send_boxL = "); RBoxShow(&send_boxL);
    printLog("recv_boxR = "); RBoxShow(&recv_boxR);
    printLog("send_boxR = "); RBoxShow(&send_boxR);
    printLog("recv_boxL = "); RBoxShow(&recv_boxL);
    printLog ("countL, R = %d, %d\n", countL, countR);
    #endif

    // ----------------------------------------------------
    // 1C. Send / recv data
    // ----------------------------------------------------

    #pragma acc host_data use_device(send_bufL, recv_bufL, send_bufR, recv_bufR)
    {
      int tag1 = 10 + 2*dir;
      int tag2 = 10 + 2*dir + 1;
      MPI_Sendrecv (send_bufL, countL, MPI_UINT32_T, rankL, tag1,
                    recv_bufR, countL, MPI_UINT32_T, rankR, tag1,  
                    MPI_COMM_WORLD, &status);
      MPI_Sendrecv (send_bufR, countR, MPI_UINT32_T, rankR, tag2,
                    recv_bufL, countR, MPI_UINT32_T, rankL, tag2,  
                    MPI_COMM_WORLD, &status);
    }

    // ----------------------------------------------------
    // 1D. Unpack data, and "or" them 
    // ----------------------------------------------------

    if (rankL != MPI_PROC_NULL){
      #pragma acc parallel loop collapse(3) \
                  present(flag, recv_bufL[:size], recv_boxL)
      BOX_LOOP(recv_boxL,i,j,k){
        ind = BOX_INDEX(recv_boxL,i,j,k);
        if constexpr (flag_op == FLAG_OR) 
          flag(i,j,k) |= recv_bufL[ind];
        else if constexpr (flag_op == FLAG_AND)
          flag(i,j,k) &= recv_bufL[ind];
        else if constexpr (flag_op == FLAG_COPY)
          flag(i,j,k)  = recv_bufL[ind];

      }
    }

    if (rankR != MPI_PROC_NULL){
      #pragma acc parallel loop collapse(3) \
                  present(flag, recv_bufR[:size], recv_boxR)
      BOX_LOOP(recv_boxR,i,j,k){
        ind = BOX_INDEX(recv_boxR,i,j,k);
        if constexpr (flag_op == FLAG_OR) 
          flag(i,j,k) |= recv_bufR[ind];
        else if constexpr (flag_op == FLAG_AND)
          flag(i,j,k)  &= recv_bufR[ind];
        else if constexpr (flag_op == FLAG_COPY)
          flag(i,j,k)  = recv_bufR[ind];
      }
    }

    // 1E. Delete data on GPU 

    #pragma acc exit data delete (send_boxL, send_boxR, recv_boxL, recv_boxR)

  } // End loop on dir

}
#endif // PARALLEL

// int IsFlagged(int i, int k, const Grid *grid)
// {
//   double gradp;
//   double *x1 = grid->xC[IDIR];
//   double *x3 = grid->xC[KDIR];
//   double r  = sqrt(x1[i]*x1[i] + x3[k]*x3[k]);
//   double eps = grid->dx[IDIR][0];
  
//   if (r < (2.0 + eps) && r > (2.0 - eps)) gradp = 1.e3;
//   else gradp = 0.0;

//   return gradp > 0;
// }
