////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Exchange boundary conditions with neighbour processor.

  Exchange data values between processors using non blocking 
  MPI_Isend/recv()

TODO: 
Add NCCL implementation 
  
  \author A. Mignone  (andrea.mignone@unito.it)
          M. Rossazza (marco.rossazza@unito.it)
  \date   Nov 28, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if (defined PARALLEL) && (PARALLEL_COMM_ASYNC == YES)
#define DEBUG  NO

////////////////////////////////////////////////////////////////////////
/*!
  Set boundary conditions on one or more sides of the computational
  domain.
  
  \param [in,out] d  pointer to PLUTO Data structure containing the
  solution array (including centered and staggered
  fields)
  \param [in]  grid   pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
void ParallelExchange (const Ary4D &Vc, const Ary3D (&Vs)[6], Grid *grid)
{
  int i, j, k, nv, nvs;
  int ip, jp, kp;
  int dir;
  static int first_call = 1;
  #ifdef _OPENACC
  int count = 0; // async counter
  #endif

  static RBox send_box[3][3][3];
  static RBox recv_box[3][3][3];

  static RBox send_boxs[3][3][3][3];
  static RBox recv_boxs[3][3][3][3];

  static double *send_buf[3][3][3];  
  static double *recv_buf[3][3][3];  

  static int send_rnk[3][3][3];
  static int recv_rnk[3][3][3];

  MPI_Request send_req[3][3][3];
  MPI_Request recv_req[3][3][3];
  
  static int beg[3], end[3];

  int coords[3], rnk;
  MPI_Comm cartcomm;
#if USE_AL == NO
cartcomm = new_SZ.cart_comm;
#else
  AL_Get_cart_comm(SZ, &cartcomm); 
#endif

  // ------------------------------------------------------
  //  0. Define the domain of neghbouring processes:
  //     for every direction we will loop from -1 (if it exists 
  //     a process to the left) to +1 (if it exists a process 
  //     to the right).
  // ------------------------------------------------------

  if (first_call){
    DIM_LOOP(dir){
      beg[dir] = 0;
      end[dir] = 0;
      if (PROC_EXIST(grid->lbound[dir], grid->nproc[dir])) beg[dir] = -1;
      if (PROC_EXIST(grid->rbound[dir], grid->nproc[dir])) end[dir] = +1;
    }
  }

  // ------------------------------------------------------
  // 1. Define send / recv boxes.
  //    Box dimension is uniquely set by the triplet that 
  //    identifies:
  //
  //    (-1, 0, 0)    [ibeg, ibeg + ngh][jbeg, jend][kbeg,kend]
  //    (+1, 0, 0)    [iend-nght, iend]
  //
  //    (-1,-1, 0)    (-1,0,0), jend -= ny-ngh 
  //    (-1,+1, 0)    (-1,0,0), jbeg += ny-ngh 
  //    (+1,-1, 0)    (+1,0,0), jend -= ny-ngh 
  //    (+1,+1, 0)    (+1,0,0), jbeg += ny-ngh
  //
  //    (-1,-1,-1)    (-1,-1,0), kend -= nz-ngh
  //    (-1,-1,+1)    (-1,-1,0), kbeg += nz-ngh
  //    (-1,+1,-1)    (-1,+1,0), kend -= nz-ngh
  //    (-1,+1,+1)    (-1,+1,0), kbeg += nz-ngh
  // ------------------------------------------------------
  
  if (first_call){
    for (kp = beg[KDIR]; kp <= end[KDIR]; kp++){
    for (jp = beg[JDIR]; jp <= end[JDIR]; jp++){
    for (ip = beg[IDIR]; ip <= end[IDIR]; ip++){
      
      int iv = ip+1, jv = jp+1 , kv = kp+1 ; // Map from ip to iv, etc...
      if (ip == 0 && jp == 0 && kp == 0) continue; // No need to allocate this buffer
      
      int cbeg[3], cend[3], icoord[3];

      icoord[IDIR] = ip; // coordinate increments
      icoord[JDIR] = jp;
      icoord[KDIR] = kp;

      // --------------------------------------------------
      // 1A. Define send boxes for central variables
      // --------------------------------------------------

      DIM_LOOP(dir) {  
        int lbeg = grid->lbeg[dir];
        int lend = grid->lend[dir];
        int nbuf = grid->nghost[dir]-1;

        cbeg[dir] = (icoord[dir] == +1 ? lend-nbuf:lbeg); 
        cend[dir] = (icoord[dir] == -1 ? lbeg+nbuf:lend); 

      }
      
      RBoxDefine (cbeg, cend, &(send_box[iv][jv][kv]));
    
      // --------------------------------------------------
      //  1B. Define staggered boxes (one per field component).
      //      When passing staggered field component in the direction of 
      //      communication (e.g. BX1s when dir = IDIR), since the same 
      //      interface is shared by two neighbouring processors, we 
      //      adopt the convention that the processor "owns" its
      //      leftmost active face.
            
      //      When communicating between processors, for instance,
            
      //      [RankL]
      //    |----|----|....|....|            
      //              o    o    o
      //              ^    ^    ^
      //              *    *    *
      //    |....|....|----|----|----|----|----|----|....|....|
      //                                  *    *
      //                                  v    v
      //                                  o    o
      //                                  |....|....|----|----|  
      //                                              [RankR]                              

      //      values denoted by a '*' are sent to the processor to the left 
      //      or to the right and received in the corresponding 'o' faces.
      //      Thus, the sendboxL includes one extra point
      //      (3 in the diagram) while sendBoxR has 2 points in this example. 
      //      Fields transverse to the direction of exchange gets 1 additional 
      //      points, in the spatial component where they are staggered. 
      //      i.e. X2-component of BX2s has 1 additional point in this case.
      //      In the case of DIVE_CONTROL == CONSTRAINED_TRANSPORT the boxes of
      //      the electric fields are  
      // --------------------------------------------------

      #ifdef STAGGERED_MHD
      DIM_LOOP(nvs){
        
        int sbeg[3], send[3];
        RBoxGetIndexes(&(send_box[iv][jv][kv]), sbeg, send);

        DIM_LOOP(dir) send[dir] += (icoord[dir] == -1 && nvs == dir ? 1:0);           
        
        DIM_LOOP(dir) send[dir] += (icoord[dir] == 0 && nvs == dir ? 1:0); 

        RBoxDefine(sbeg, send, &(send_boxs[nvs][iv][jv][kv]));
      }
      #endif
      
      // --------------------------------------------------
      // 1C. Get rank of processor to which we send
      // --------------------------------------------------

      coords[IDIR] = grid->rank_coord[IDIR];
      coords[JDIR] = grid->rank_coord[JDIR];
      coords[KDIR] = grid->rank_coord[KDIR];

      send_rnk[iv][jv][kv] = MPI_PROC_NULL;
      
      DIM_LOOP(dir){
        if(icoord[dir]==-1){
          if (PROC_EXIST(grid->lbound[dir], grid->nproc[dir])){
            coords[dir] = grid->rank_coord[dir] - 1;
            MPI_Cart_rank(cartcomm, coords, &rnk);
            send_rnk[iv][jv][kv] = rnk;
          }
        }

        if(icoord[dir]==+1){
          if (PROC_EXIST(grid->rbound[dir], grid->nproc[dir])){
            coords[dir] = grid->rank_coord[dir] + 1;
            MPI_Cart_rank(cartcomm, coords, &rnk);
            send_rnk[iv][jv][kv] = rnk;
          }
        }
      }

      // --------------------------------------------------
      // 1D. Define recv boxes for central variables
      // --------------------------------------------------

      DIM_LOOP(dir) { 
        int lbeg = grid->lbeg[dir];
        int lend = grid->lend[dir];
        int nbuf = grid->nghost[dir]-1;

        if (icoord[dir] == -1) {
          cbeg[dir] = 0;
          cend[dir] = lbeg-1;
        }
        if (icoord[dir] == +1) {
          cbeg[dir] = lend+1;
          cend[dir] = lend+1+nbuf;
        }
      }

      RBoxDefine (cbeg, cend, &(recv_box[iv][jv][kv]));
    
      // --------------------------------------------------
      // 1E. Define recv boxes for staggered variables
      // --------------------------------------------------

      #ifdef STAGGERED_MHD
      int sbeg[3], send[3];
      
      DIM_LOOP(nvs){
        RBoxGetIndexes(&(recv_box[iv][jv][kv]), sbeg, send);

        DIM_LOOP(dir) send[dir] += (icoord[dir] == +1 && nvs == dir ? 1:0); 

        DIM_LOOP(dir) send[dir] += (icoord[dir] == 0 && nvs == dir ? 1:0); 
        
        RBoxDefine(sbeg, send, &(recv_boxs[nvs][iv][jv][kv]));
      }
      #endif
      
      // --------------------------------------------------
      // 1F. Get rank of processor from which we receive
      // --------------------------------------------------

      recv_rnk[iv][jv][kv] = send_rnk[iv][jv][kv];

    }}}
  }

  #pragma omp target enter data map(to: send_box[:3][:3][:3], recv_box[:3][:3][:3], \
                                    send_boxs[:3][:3][:3][:3], recv_boxs[:3][:3][:3][:3])

  #pragma acc enter data copyin(send_box)
  #pragma acc enter data copyin(recv_box)
  #pragma acc enter data copyin(send_boxs)
  #pragma acc enter data copyin(recv_boxs)

  // ------------------------------------------------------
  // 1G. Allocate buffers
  // ------------------------------------------------------

  if (first_call){
    for (kp = beg[KDIR]; kp <= end[KDIR]; kp++){
    for (jp = beg[JDIR]; jp <= end[JDIR]; jp++){
    for (ip = beg[IDIR]; ip <= end[IDIR]; ip++){
      int iv = ip+1, jv = jp+1 , kv = kp+1 ; 
      if (ip == 0 && jp == 0 && kp == 0) continue; 
      
      long int send_size = NVAR*RBoxGetSize(&(send_box[iv][jv][kv])); 
      DIM_LOOP(nvs) send_size += RBoxGetSize(&(send_boxs[nvs][iv][jv][kv]));
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_LOOP(nvs) send_size += RBoxGetSize(&(send_boxs[nvs][iv][jv][kv]));
      #endif
      send_buf[iv][jv][kv] = MALLOC_DOUBLE(send_size);

      long int recv_size = NVAR*RBoxGetSize(&(recv_box[iv][jv][kv])); 
      DIM_LOOP(nvs) recv_size += RBoxGetSize(&(recv_boxs[nvs][iv][jv][kv]));
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_LOOP(nvs) recv_size += RBoxGetSize(&(recv_boxs[nvs][iv][jv][kv]));
      #endif
      recv_buf[iv][jv][kv] = MALLOC_DOUBLE(recv_size);

    }}}
  }

  // ------------------------------------------------------
  // 2. Begin non-blocking receives
  // ------------------------------------------------------

  for (kp = end[KDIR]; kp >= beg[KDIR]; kp--){
  for (jp = end[JDIR]; jp >= beg[JDIR]; jp--){
  for (ip = end[IDIR]; ip >= beg[IDIR]; ip--){
    if (ip == 0 && jp == 0 && kp == 0) continue; 
    int iv = ip+1, jv = jp+1 , kv = kp+1 ; 

    long int recv_size = NVAR*RBoxGetSize(&(recv_box[iv][jv][kv])); 
    DIM_LOOP(nvs) recv_size += RBoxGetSize(&(recv_boxs[nvs][iv][jv][kv]));
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_LOOP(nvs) recv_size += RBoxGetSize(&(recv_boxs[nvs][iv][jv][kv]));
    #endif

    MPI_Irecv(recv_buf[iv][jv][kv], recv_size, MPI_DOUBLE, recv_rnk[iv][jv][kv], 
              100, MPI_COMM_WORLD, &(recv_req[iv][jv][kv]));
  }}}

  // ------------------------------------------------------
  // 3. Pack buffers 
  // ------------------------------------------------------

  for (kp = beg[KDIR]; kp <= end[KDIR]; kp++){
  for (jp = beg[JDIR]; jp <= end[JDIR]; jp++){
  for (ip = beg[IDIR]; ip <= end[IDIR]; ip++){
    if (ip == 0 && jp == 0 && kp == 0) continue; 
    int iv = ip+1, jv = jp+1 , kv = kp+1 ; 
    long int ind;

    double *buf = send_buf[iv][jv][kv];
    RBox   &box = send_box[iv][jv][kv];
    long int size = RBoxGetSize(&box);

    #pragma omp target teams distribute parallel for collapse(4) is_device_ptr(buf) 
    #pragma acc parallel loop collapse(4) deviceptr(buf) present(Vc, box) \
                                          async(count++)
    NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
      ind = nv*size + BOX_INDEX(box,i,j,k);
      buf[ind] = Vc(i,j,k,nv);
    }

    #ifdef STAGGERED_MHD
    long int offset = NVAR*size;
    DIM_LOOP(nvs){
      RBox &boxs = send_boxs[nvs][iv][jv][kv];
      long int sizes = RBoxGetSize(&boxs);

      #pragma omp target teams distribute parallel for collapse(3) is_device_ptr(buf)
      #pragma acc parallel loop collapse(3) deviceptr(buf) present(Vs, boxs) \
                                            async(count++)
      BOX_LOOP(boxs,i,j,k){
        ind = offset + BOX_INDEX(boxs,i,j,k);
        buf[ind] = Vs[nvs](i,j,k);
        #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
        // after BX*s we fill the buf correspondent to EX*s 
        ind = offset +  sizes + BOX_INDEX(boxs,i,j,k);
        buf[ind] = Vs[nvs+DIMENSIONS](i,j,k);
        #endif
      }
      offset += sizes;
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      offset += sizes;
      #endif
    }
    #endif
  }}}
  #pragma acc wait

  // ------------------------------------------------------
  // 4a. Send buffers 
  // ------------------------------------------------------

  for (kp = beg[KDIR]; kp <= end[KDIR]; kp++){
  for (jp = beg[JDIR]; jp <= end[JDIR]; jp++){
  for (ip = beg[IDIR]; ip <= end[IDIR]; ip++){
    if (ip == 0 && jp == 0 && kp == 0) continue; 
    int iv = ip+1, jv = jp+1 , kv = kp+1 ; 

    long int send_size = NVAR*RBoxGetSize(&(send_box[iv][jv][kv])); 
    DIM_LOOP(nvs) send_size += RBoxGetSize(&(send_boxs[nvs][iv][jv][kv]));
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_LOOP(nvs) send_size += RBoxGetSize(&(send_boxs[nvs][iv][jv][kv]));
    #endif

    MPI_Isend(send_buf[iv][jv][kv], send_size, MPI_DOUBLE, send_rnk[iv][jv][kv], 
              100, MPI_COMM_WORLD, &(send_req[iv][jv][kv]));
  }}}

  // ------------------------------------------------------
  // 4b. MPI_Waitall on MPI_Irecv
  // ------------------------------------------------------

  for (kp = end[KDIR]; kp >= beg[KDIR]; kp--){
  for (jp = end[JDIR]; jp >= beg[JDIR]; jp--){
  for (ip = end[IDIR]; ip >= beg[IDIR]; ip--){
    if (ip == 0 && jp == 0 && kp == 0) continue; 
    int iv = ip+1, jv = jp+1 , kv = kp+1 ; 
  
    MPI_Waitall(1, &(recv_req[iv][jv][kv]), MPI_STATUSES_IGNORE); 
  }}}

  // ------------------------------------------------------
  // 4c. Unpack buffers
  // ------------------------------------------------------
  
  #ifdef _OPENACC
  count = 0;
  #endif
  for (kp = beg[KDIR]; kp <= end[KDIR]; kp++){
  for (jp = beg[JDIR]; jp <= end[JDIR]; jp++){
  for (ip = beg[IDIR]; ip <= end[IDIR]; ip++){
    if (ip == 0 && jp == 0 && kp == 0) continue; 
    int iv = ip+1, jv = jp+1 , kv = kp+1 ; 
    long int ind;

    double *buf = recv_buf[iv][jv][kv];
    RBox   &box = recv_box[iv][jv][kv];
    long int size = RBoxGetSize(&box);
    
    #pragma omp target teams distribute parallel for collapse(4) is_device_ptr(buf)
    #pragma acc parallel loop collapse(4) deviceptr(buf) present(Vc, box) \
                                          async(count++)
    NVAR_LOOP(nv) BOX_LOOP(box,i,j,k){
      ind = nv*size + BOX_INDEX(box,i,j,k);
      Vc(i,j,k,nv) = buf[ind]; 
    }

    #ifdef STAGGERED_MHD
    long int offset = NVAR*size;
    DIM_LOOP(nvs){
      RBox &boxs = recv_boxs[nvs][iv][jv][kv];
      long int sizes = RBoxGetSize(&boxs);

      #pragma omp target teams distribute parallel for collapse(3) is_device_ptr(buf)
      #pragma acc parallel loop collapse(3) deviceptr(buf) present(Vs, boxs) \
                                            async(count++)
      BOX_LOOP(boxs,i,j,k){
        ind = offset + BOX_INDEX(boxs,i,j,k);
        Vs[nvs](i,j,k) = buf[ind];
        #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
        ind = offset +  sizes + BOX_INDEX(boxs,i,j,k);
        Vs[nvs+DIMENSIONS](i,j,k) = buf[ind];
        #endif
      }
      offset += sizes;
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      offset += sizes;
      #endif
    }
    #endif
  }}}
  #pragma acc wait

  // ------------------------------------------------------
  // 5. MPI_Waitall on MPI_Isend
  // ------------------------------------------------------

  for (kp = beg[KDIR]; kp <= end[KDIR]; kp++){
  for (jp = beg[JDIR]; jp <= end[JDIR]; jp++){
  for (ip = beg[IDIR]; ip <= end[IDIR]; ip++){
    if (ip == 0 && jp == 0 && kp == 0) continue; 
    int iv = ip+1, jv = jp+1 , kv = kp+1 ; 
    
    MPI_Waitall(1, &(send_req[iv][jv][kv]), MPI_STATUSES_IGNORE); 
  }}}

  first_call = 0;

  #ifndef STAGGERED_MHD
  #endif
}
#endif  // (defined PARALLEL) && (PARALLEL_COMM_ASYNC == YES) 
