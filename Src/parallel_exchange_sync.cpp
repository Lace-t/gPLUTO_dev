////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Exchange boundary conditions with neighbour processor.

  Exchange data values between processors using blocking MPI_Sendrecv()

TODO: 
Allow also for staggered electric field components.
Add NCCL implementation 

  
  \author A. Mignone  (andrea.mignone@unito.it)
          M. Rossazza (marco.rossazza@unito.it)
  \date   Feb 02, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define DEBUG  NO

#if (defined PARALLEL) && (PARALLEL_COMM_ASYNC == NO)
////////////////////////////////////////////////////////////////////////
void ParallelExchange (const Ary4D &Vc, const Ary3D (&Vs)[6], Grid *grid)
/*!
  Set boundary conditions on one or more sides of the computational
  domain.
  
  \param [in,out] d  pointer to PLUTO Data structure containing the
  solution array (including centered and staggered
  fields)
  \param [in]  grid   pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
{
  int  i,j,k,nv;
  int  par_dim[3] = {0, 0, 0};

  int  nx1_tot = grid->ntot[IDIR];
  int  nx2_tot = grid->ntot[JDIR];
  int  nx3_tot = grid->ntot[KDIR];
  int  n1g = grid->nghost[IDIR];
  int  n2g = grid->nghost[JDIR];
  int  n3g = grid->nghost[KDIR];

  MPI_Status status;

  static int **nrnks;
  static int ***nrnks_edg;
  static int ***nrnks_crn;

  static double *send_bufL;
  static double *recv_bufL;
  static double *send_bufR;
  static double *recv_bufR;

  int nmax = nx1_tot*nx2_tot;
  nmax=MAX(nmax,nx2_tot*nx3_tot);
  nmax=MAX(nmax,nx1_tot*nx3_tot);

  int nmaxg = n1g;
  nmaxg = MAX(nmaxg,n2g);
  nmaxg = MAX(nmaxg,n3g);

  int size = NVAR*nmaxg*nmax;

// --------------------------------------------------------
// 0a. Allocate arrays
// -------------------------------------------------------- 

  if (nrnks == NULL){
    nrnks     = ARRAY_2D(3, 2, int);
    nrnks_edg = ARRAY_3D(3, 2, 2, int);
    nrnks_crn = ARRAY_3D(2, 2, 2, int);

    send_bufL = ARRAY_1D(size, double);  // Arrays for cell-centered buffers 
    recv_bufL = ARRAY_1D(size, double);
    send_bufR = ARRAY_1D(size, double);
    recv_bufR = ARRAY_1D(size, double);

    #pragma acc enter data create(send_bufL[:size], recv_bufL[:size])
    #pragma acc enter data create(send_bufR[:size], recv_bufR[:size])

    GetNeighbourRanks(grid, nrnks, nrnks_edg, nrnks_crn);
  }

  // --------------------------------------------------------
  // 0b. Check the number of processors in each direction
  // --------------------------------------------------------

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

  // ------------------------------------------------------------------
  // 1. Loop over directions.
  // Box naming convention:
  // X         X         X                   
  // ........ |---------|---------|---------|---------|......... 
  // |recv_boxL|send_boxL|                   |send_boxR|recv_boxR
  // 
  // Note that recv_box'es are obtained by simply offsetting the
  // send_box'es.
  // ------------------------------------------------------------------


  int dir;
  for (dir = 0; dir < 3; dir++){
    
    if (!par_dim[dir]) continue;  // Direction is not parallel, skip
    if (DEBUG) printLog ("----> dir = %d\n",dir);

    int l;
    int beg[3], end[3];
    long ind, countL, countR;
    int rankL = nrnks[dir][0];
    int rankR = nrnks[dir][1];
    RBox send_boxL, send_boxR;
    RBox recv_boxL, recv_boxR; 
    RBox *box;

    // --------------------------------------------
    // 1. Construct send / recv boxes 
    // for cell-centered data.
    // --------------------------------------------

  // -- 1a. Define send_boxL box -- 

    for (l = 0; l < 3; l++) {  // Initialize to whole domain
      beg[l] = 0;
      end[l] = grid->ntot[l] - 1;
    }
    beg[dir] = grid->lbeg[dir] ;
    end[dir] = beg[dir] + grid->nghost[dir] - 1;

    RBoxDefine (beg, end, &send_boxL);
    
  // -- 1b. Define recv_boxR -- 

    beg[dir] += grid->nint[dir];
    end[dir] += grid->nint[dir];

    RBoxDefine (beg, end, &recv_boxR);

  // -- 1c. Define send_boxR -- 

    beg[dir] -= grid->nghost[dir];
    end[dir] -= grid->nghost[dir];

    RBoxDefine (beg, end, &send_boxR);

  // -- 1d. Define recv_boxL -- 

    beg[dir] = 0;
    end[dir] = grid->nghost[dir]-1;

    RBoxDefine (beg, end, &recv_boxL);

    #pragma acc enter data copyin(send_boxL, send_boxR, recv_boxL, recv_boxR)

    // ------------------------------------------
    // 2. Pack cell-centered data 
    // ------------------------------------------

    // if (DEBUG) {
    //   printLog ("rankL = %d, rankR = %d\n", rankL, rankR);
    //   printLog("send_boxL = "); RBoxShow(&send_boxL);
    //   printLog("recv_boxR = "); RBoxShow(box);
    //   printLog("send_boxR = "); RBoxShow(&send_boxR);
    //   printLog("recv_boxL = "); RBoxShow(&recv_boxL);
    //   printLog ("countL, R = %d, %d\n", countL, countR);
    // }

    for (nv = 0; nv < NVAR; nv++){

      box = &send_boxL;
      #pragma acc parallel loop collapse(3) present(Vc, send_bufL[:size], box)
      BOX_LOOP(*box,i,j,k) {
        ind = BOX_INDEX(*box,i,j,k);
        send_bufL[ind] = Vc(i,j,k, nv);
      }
      countL = box->nx1*box->nx2*box->nx3; 

      box = &send_boxR;
      #pragma acc parallel loop collapse(3) present(Vc, send_bufR[:size], box)
      BOX_LOOP(*box,i,j,k) {
        ind = BOX_INDEX(*box,i,j,k);
        send_bufR[ind] = Vc(i,j,k, nv);
      }
      countR = box->nx1*box->nx2*box->nx3; 

      // ------------------------------------------
      // 3. Send / recv data
      // ------------------------------------------

      #pragma acc host_data use_device(send_bufL, recv_bufL, send_bufR, recv_bufR)
      {
      MPI_Sendrecv (send_bufL, countL, MPI_DOUBLE, rankL, 1,
                    recv_bufR, countL, MPI_DOUBLE, rankR, 1,  MPI_COMM_WORLD, &status);
      MPI_Sendrecv (send_bufR, countR, MPI_DOUBLE, rankR, 2,
                    recv_bufL, countR, MPI_DOUBLE, rankL, 2,  MPI_COMM_WORLD, &status);
      }

      // ------------------------------------------
      // 4. Unpack cell-centered data 
      // ------------------------------------------

      if (rankL != MPI_PROC_NULL){
        box = &recv_boxL;
        #pragma acc parallel loop collapse(3) present(Vc, recv_bufL[:size], box)
        BOX_LOOP(*box,i,j,k){
          ind = BOX_INDEX(*box,i,j,k);
          Vc(i,j,k,nv) = recv_bufL[ind];
        }
      }

      if (rankR != MPI_PROC_NULL){
        box = &recv_boxR;
        #pragma acc parallel loop collapse(3) present(Vc, recv_bufR[:size], box)
        BOX_LOOP(*box,i,j,k){
          ind = BOX_INDEX(*box,i,j,k);
          Vc(i,j,k,nv) = recv_bufR[ind];
        }
      }
    } // End Loop on NVAR 

    // --------------------------------------------------- 
    // 5. Staggered Fields
    // ----------------------------------------------------
 
    #ifdef STAGGERED_MHD

    // -------------------------------------------------------------------
    // 6. Define staggered boxes (one per field component).
    // When passing staggered field component in the direction of 
    // communication (e.g. BX1s when dir = IDIR), since the same 
    // interface is shared by two neighbouring processors, we 
    // adopt the convention that the processor "owns" its
    // leftmost active face.
    // 
    // When communicating between processors, for instance,
    // 
    // [RankL]
    // |----|----|....|....|            
    // o    o    o
    // ^    ^    ^
    // *    *    *
    // |....|....|----|----|----|----|----|----|....|....|
    // *    *
    // v    v
    // o    o
    // |....|....|----|----|  
    // [RankR]                              
    // 
    // values denoted by a '*' are sent to the processor to the left 
    // or to the right and received in the corresponding 'o' faces.
    // Thus, the sendboxL includes one extra point
    // (3 in the diagram) while sendBoxR has 2 points. 
    // For transverse components, the sendbuf always contains
    // one extra values.
    // Summarizing, the sendboxes always extend to end+1 except for 
    // the sendBoxR in the normal direction which ends at end.
    // -------------------------------------------------------------------

    RBox send_stag_boxL[3], send_stag_boxR[3];
    RBox recv_stag_boxL[3], recv_stag_boxR[3];

    #pragma acc enter data create(send_stag_boxL[:3], send_stag_boxR[:3], \
                                  recv_stag_boxL[:3], recv_stag_boxR[:3])

    // -------------------------------------------- 
    // 7. Loop over staggered components 
    // --------------------------------------------

    int staggered_fields = 3; // DIMENSIONS;
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    staggered_fields += 3; // DIMENSIONS;
    #endif
    
    constexpr bool include_dirs[] = {INCLUDE_IDIR, INCLUDE_JDIR, INCLUDE_KDIR};
    for (int ksf = 0; ksf < staggered_fields; ksf++){

      int ks = ksf%3; 
      if (!include_dirs[ks]) continue;

    // -- 7a. Define boxes: send_stag_bufL / recv_stag_bufR -- 

      RBoxGetIndexes(&send_boxL, beg, end);
      end[ks] += 1;
      RBoxDefine (beg, end, &(send_stag_boxL[ks]));  

      beg[dir] += grid->nint[dir];
      end[dir] += grid->nint[dir];
      RBoxDefine (beg, end, &(recv_stag_boxR[ks]));  

    // -- 7b. Define boxes: send_stag_bufR /recv_stag_bufL -- 

      RBoxGetIndexes(&send_boxR, beg, end);
      // end[ks]  += 1; // (ks != dir);  We need += 1 since at restart time, 
      //                // the leftmost interface is not correctly assigned.
      if (ks != dir) end[ks] += 1;
      RBoxDefine (beg, end, &(send_stag_boxR[ks]));  

      beg[dir] -= grid->nint[dir];
      end[dir] -= grid->nint[dir];
      RBoxDefine (beg, end, &(recv_stag_boxL[ks]));  

      if (DEBUG){
        printLog ("dir = %d, ks = %d, send_stag_boxL:", dir, ks);
        RBoxShow(&(send_stag_boxL[ks]));
        printLog ("dir = %d, ks = %d, recv_stag_boxR:", dir, ks);
        RBoxShow(&(recv_stag_boxR[ks]));

        printLog ("dir = %d, ks = %d, recv_stag_boxL:", dir, ks);
        RBoxShow(&(recv_stag_boxL[ks]));
        printLog ("dir = %d, ks = %d, send_stag_boxR:", dir, ks);
        RBoxShow(&(send_stag_boxR[ks]));
      }

      #pragma acc update device (send_stag_boxL[:3], send_stag_boxR[:3], \
                                 recv_stag_boxL[:3], recv_stag_boxR[:3])

    // -- 7c. Pack data -- 

      box = &(send_stag_boxL[ks]);
      #pragma acc parallel loop collapse(3) present(Vs, send_bufL[:size], box)
      BOX_LOOP(*box,i,j,k) {
        ind = BOX_INDEX(*box,i,j,k);
        send_bufL[ind] = Vs[ksf](i,j,k);
      }
      countL = box->nx1*box->nx2*box->nx3; 

      box = &(send_stag_boxR[ks]);
      #pragma acc parallel loop collapse(3) present(Vs, send_bufR[:size], box)
      BOX_LOOP(*box,i,j,k) {
        ind = BOX_INDEX(*box,i,j,k);
        send_bufR[ind] = Vs[ksf](i,j,k);
      }
      countR = box->nx1*box->nx2*box->nx3; 

    // -- 7d. Send / recv data -- 

      #pragma acc host_data use_device(send_bufL, recv_bufL, send_bufR, recv_bufR)
      {
      MPI_Sendrecv (send_bufL, countL, MPI_DOUBLE, rankL, 1,
                    recv_bufR, countL, MPI_DOUBLE, rankR, 1, MPI_COMM_WORLD, &status);
      MPI_Sendrecv (send_bufR, countR, MPI_DOUBLE, rankR, 2,
                    recv_bufL, countR, MPI_DOUBLE, rankL, 2, MPI_COMM_WORLD, &status);
      }

    // -- 7e. Unpack data -- 

      if (rankL != MPI_PROC_NULL){
        box = &(recv_stag_boxL[ks]);
        #pragma acc parallel loop collapse(3) present(Vs, recv_bufL[:size], box)
        BOX_LOOP(*box,i,j,k){
          ind = BOX_INDEX(*box,i,j,k);
          Vs[ksf](i,j,k) = recv_bufL[ind];
        }
      }

      if (rankR != MPI_PROC_NULL){
        box = &(recv_stag_boxR[ks]);
        #pragma acc parallel loop collapse(3) present(Vs, recv_bufR[:size], box)
        BOX_LOOP(*box,i,j,k){
          ind = BOX_INDEX(*box,i,j,k);
          Vs[ksf](i,j,k) = recv_bufR[ind];
        }
      }
    }  // End Loop on ksf 
    #endif // STAGGERED_MHD 

  // -- Delete data on GPU -- 

    #pragma acc exit data delete (send_boxL, send_boxR, recv_boxL, recv_boxR)
    #ifdef STAGGERED_MHD
    #pragma acc exit data delete (send_stag_boxL[:3], send_stag_boxR[:3], \
                                  recv_stag_boxL[:3], recv_stag_boxR[:3])
    #endif

  } // End loop on dir    
}
#endif  // (defined PARALLEL) && (PARALLEL_COMM_ASYNC == NO) 

#ifdef PARALLEL
////////////////////////////////////////////////////////////////////////
void ParallelExchangeDouble (Ary3D &flag, Grid *grid)
/*!
  Set boundary conditions on one or more sides of the computational
  domain.
  
  \param [in,out] d  pointer to PLUTO Data structure containing the
  solution array (including centered and staggered
  fields)
  \param [in]  grid   pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
{
  int  i,j,k;
  int  par_dim[3] = {0, 0, 0};

  int  nx1_tot = grid->ntot[IDIR];
  int  nx2_tot = grid->ntot[JDIR];
  int  nx3_tot = grid->ntot[KDIR];
  int  n1g = grid->nghost[IDIR];
  int  n2g = grid->nghost[JDIR];
  int  n3g = grid->nghost[KDIR];

  MPI_Status status;

  static int **nrnks;
  static int ***nrnks_edg;
  static int ***nrnks_crn;

  static double *send_bufL;
  static double *recv_bufL;
  static double *send_bufR;
  static double *recv_bufR;

  int nmax = nx1_tot*nx2_tot;
  nmax=MAX(nmax,nx2_tot*nx3_tot);
  nmax=MAX(nmax,nx1_tot*nx3_tot);

  int nmaxg = n1g;
  nmaxg = MAX(nmaxg,n2g);
  nmaxg = MAX(nmaxg,n3g);

  int size  = nmaxg*nmax;

  // ------------------------------------------------------
  // 0a. Allocate arrays
  // ------------------------------------------------------

  if (nrnks == NULL){
    nrnks    = ARRAY_2D(3, 2, int);
    nrnks_edg = ARRAY_3D(3, 2, 2, int);
    nrnks_crn = ARRAY_3D(2, 2, 2, int);

    send_bufL = ARRAY_1D(size, double);  // Arrays for cell-centered buffers 
    recv_bufL = ARRAY_1D(size, double);
    send_bufR = ARRAY_1D(size, double);
    recv_bufR = ARRAY_1D(size, double);

    #pragma omp target enter data map(alloc: send_bufL[:size], recv_bufL[:size], \
                                      send_bufR[:size], recv_bufR[:size])
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
  //    Note that recv_box'es are obtained by simply offsetting the
  //    send_box'es.
  // ------------------------------------------------------

  int dir;
  for (dir = 0; dir < 3; dir++){
    
    if (!par_dim[dir]) continue;  // Direction is not parallel, skip
    if (DEBUG) printLog ("----> dir = %d\n",dir);

    int l;
    int beg[3], end[3];
    long countL, countR;
    int rankL = nrnks[dir][0];
    int rankR = nrnks[dir][1];
    RBox send_boxL, send_boxR;
    RBox recv_boxL, recv_boxR; 
    RBox *box;

    // 1a. Construct send / recv boxes for cell-centered data:
    //     Define send_boxL 

    for (l = 0; l < 3; l++) {  // Initialize to whole domain
      beg[l] = 0;
      end[l] = grid->ntot[l] - 1;
    }
    beg[dir] = grid->lbeg[dir] ;
    end[dir] = beg[dir] + grid->nghost[dir] - 1;

    RBoxDefine (beg, end, &send_boxL);
    
    // 1b. Define recv_boxR

    beg[dir] += grid->nint[dir];
    end[dir] += grid->nint[dir];

    RBoxDefine (beg, end, &recv_boxR);

    // 1c. Define send_boxR

    beg[dir] -= grid->nghost[dir];
    end[dir] -= grid->nghost[dir];

    RBoxDefine (beg, end, &send_boxR);

    // 1d. Define recv_boxL

    beg[dir] = 0;
    end[dir] = grid->nghost[dir]-1;

    RBoxDefine (beg, end, &recv_boxL);
    
    #pragma omp target enter data map(to: send_boxL, send_boxR, recv_boxL, recv_boxR)
    #pragma acc enter data copyin(send_boxL, send_boxR, recv_boxL, recv_boxR)

    // ----------------------------------------------------
    // 1e. Pack cell-centered data 
    // ----------------------------------------------------

    // if (DEBUG) {
    //   printLog ("rankL = %d, rankR = %d\n", rankL, rankR);
    //   printLog("send_boxL = "); RBoxShow(&send_boxL);
    //   printLog("recv_boxR = "); RBoxShow(box);
    //   printLog("send_boxR = "); RBoxShow(&send_boxR);
    //   printLog("recv_boxL = "); RBoxShow(&recv_boxL);
    //   printLog ("countL, R = %d, %d\n", countL, countR);
    // }

    box = &send_boxL;
    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(flag, send_bufL[:size], box)
    BOX_LOOP(*box,i,j,k) {
      long int ind = BOX_INDEX(*box,i,j,k);
      send_bufL[ind] = flag(i,j,k);
    }
    countL = box->nx1*box->nx2*box->nx3; 

    box = &send_boxR;
    #pragma omp target teams distribute parallel for collapse(3)
    #pragma acc parallel loop collapse(3) present(flag, send_bufR[:size], box)
    BOX_LOOP(*box,i,j,k) {
      long int ind = BOX_INDEX(*box,i,j,k);
      send_bufR[ind] = flag(i,j,k);
    }
    countR = box->nx1*box->nx2*box->nx3; 

    // ----------------------------------------------------
    // 1f. Send / recv data
    // ----------------------------------------------------

    #pragma omp target data use_device_addr(send_bufL, recv_bufL, send_bufR, recv_bufR)
    #pragma acc host_data use_device(send_bufL, recv_bufL, send_bufR, recv_bufR)
    {
    MPI_Sendrecv (send_bufL, countL, MPI_DOUBLE, rankL, 1,
                  recv_bufR, countL, MPI_DOUBLE, rankR, 1,  MPI_COMM_WORLD, &status);
    MPI_Sendrecv (send_bufR, countR, MPI_DOUBLE, rankR, 2,
                  recv_bufL, countR, MPI_DOUBLE, rankL, 2,  MPI_COMM_WORLD, &status);
    }

    // ----------------------------------------------------
    // 1g. Unpack cell-centered data 
    // ----------------------------------------------------

    if (rankL != MPI_PROC_NULL){
      box = &recv_boxL;
      #pragma omp target teams distribute parallel for collapse(3)
      #pragma acc parallel loop collapse(3) present(flag, recv_bufL[:size], box)
      BOX_LOOP(*box,i,j,k){
        long int ind = BOX_INDEX(*box,i,j,k);
        flag(i,j,k) = recv_bufL[ind];
      }
    }

    if (rankR != MPI_PROC_NULL){
      box = &recv_boxR;
      #pragma omp target teams distribute parallel for collapse(3)
      #pragma acc parallel loop collapse(3) present(flag, recv_bufR[:size], box)
      BOX_LOOP(*box,i,j,k){
        long int ind = BOX_INDEX(*box,i,j,k);
        flag(i,j,k) = recv_bufR[ind];
      }
    }

   // Delete data on GPU

    #pragma omp target exit data map(delete: send_boxL, send_boxR, recv_boxL, recv_boxR)
    #pragma acc exit data delete (send_boxL, send_boxR, recv_boxL, recv_boxR)

  } // End loop on dir    
}
#endif // PARALLEL
