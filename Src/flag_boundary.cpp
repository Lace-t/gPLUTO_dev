////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief   Apply boundary condition on flag

  Symmetrize flag at physical boundaries to preserve symmetries, 
  periodicity.
  In order to correctly preserve serial / parallel equality, 
  two strategies are possible, depending on the value of parallel_op

      1. Double Parallel exchange (parallel_op = FLAG_COPY).
        This option requires two ParallelExchange():

          i) Set root flag on domain (D):
          proc 0    |...|...|---|-F-|-F-|---|...|...|     
          proc 1                    |...|...|---|---|---|---|...|...|

          ii) Apply Boundary/Parallel Exchange() to root flag: 
          proc 0    |...|...|---|-F-|-F-|---|...|...|     
          proc 1                    |.F.|...|---|---|---|---|...|...|

          iii) Add buffer zones on domain + boundary (D+B):
          proc 0    |...|.B.|-B-|-F-|-F-|-B-|.B.|...|     
          proc 1                    |.F.|.B.|-B-|---|---|---|...|...|

          iv) Apply again Boundary/Parallel Exchange() to buffered flag: 
          proc 0    |...|.B.|-B-|-F-|-F-|-B-|.B.|...|     
          proc 1                    |.F.|.B.|-B-|---|---|---|...

          (notice that the B in the 1st ghost zone on the right on proc 0 
          is being copied from proc 1)

      2. Single parallel exchange using logical or (parallel_op = FLAG_OR).
        This option take care of the buffered zones that may extend to 
        the ghots zones in one single communication:

          i) Set root flag on domain (D), add buffer:
          proc 0    |...|.B.|-B-|-F-|-F-|-B-|.B.|...|     
          proc 1                    |...|...|---|---|---|---|...|...|

          ii) Apply Boundary/Parallel Exchange() to buffered flag using logical or: 
          proc 0    |...|.B.|-B-|-F-|-F-|-B-|.B.|...|     
          proc 1                    |.F.|.B.|-B-|---|---|---|...|...|


  \authors A. Mignone (andrea.mignone@unito.it)\n
           
  \date    Jan 10, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Set physical boundary condition on flag.
*/
////////////////////////////////////////////////////////////////////////
void FlagBoundary(Ary3Duint32_t &fl, const int parallel_op, Grid *grid)
{
  int i,  j,  k;
  int ib, jb, kb;
  int ie, je, ke;
  int type[6], is;
  int sbeg = 0; 
  int send = 5;
  int side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int par_dim[3] = {0, 0, 0};
  int  *lbound = grid->lbound;
  int  *rbound = grid->rbound;
  RBox box;
  #pragma acc enter data create(box)

  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  int nx1 = grid->nint[IDIR];
  int nx2 = grid->nint[JDIR];
  int nx3 = grid->nint[KDIR];

  // ------------------------------------------------------
  // 0. Determine if a direction is parallel or not
  // ------------------------------------------------------

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

  // ------------------------------------------------------
  // 1. Exchange data between processors
  // ------------------------------------------------------

// #ifdef PARALLEL  !! HIGH ORDER ONLY
//   ParallelExchangeFlag (fl, grid);
// #endif
  #ifdef PARALLEL
  if (parallel_op == FLAG_OR)   ParallelExchangeFlag<FLAG_OR>(fl, grid);
  if (parallel_op == FLAG_AND)  ParallelExchangeFlag<FLAG_AND>(fl, grid);
  if (parallel_op == FLAG_COPY) ParallelExchangeFlag<FLAG_COPY>(fl, grid);
  #else
  #endif

  // ------------------------------------------------------
  // 2. Set physical b.c.
  // ------------------------------------------------------

  type[0] = lbound[IDIR]*INCLUDE_IDIR; type[1] = rbound[IDIR]*INCLUDE_IDIR;
  type[2] = lbound[JDIR]*INCLUDE_JDIR; type[3] = rbound[JDIR]*INCLUDE_JDIR;
  type[4] = lbound[KDIR]*INCLUDE_KDIR; type[5] = rbound[KDIR]*INCLUDE_KDIR;

  for (is = sbeg; is <= send; is++){

    if (type[is] == 0) continue;
    
    // ------------------------------------------
    // 2A. Define bc boxes, sweeping direction
    // ------------------------------------------

    ib = 0; ie = nx1_tot-1;
    jb = 0; je = nx2_tot-1;
    kb = 0; ke = nx3_tot-1;

    if (side[is] == X1_BEG) {
      ib  = 0;    ie  = ibeg-1;
    } else if (side[is] == X1_END) {
      ib  = iend+1; ie  = nx1_tot-1;
    } else if (side[is] == X2_BEG) {
      jb  = 0;    je  = jbeg-1;
    } else if (side[is] == X2_END) {
      jb  = jend+1; je  = nx2_tot-1;
    } else if (side[is] == X3_BEG) {
      kb  = 0; ke  = kbeg-1;
    } else if (side[is] == X3_END) {
      kb  = kend+1; ke  = nx3_tot-1;
    }

    RBoxDefine (ib, ie, jb, je, kb, ke, &box);

    #pragma acc update device(box)

    // ------------------------------------------
    // 2B. Apply boundary conditions.
    // ------------------------------------------

    if (type[is] == OUTFLOW) {

      #ifdef HIGH_ORDER
      if (side[is] == X1_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = 1; // fl[k][j][ibeg];
      }else if (side[is] == X1_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = 1; // fl[k][j][iend];
      }else if (side[is] == X2_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = 1; // fl[k][jbeg][i];
      }else if (side[is] == X2_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = 1; // fl[k][jend][i];
      }else if (side[is] == X3_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j,kbeg);
      }else if (side[is] == X3_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j,kend);
      }
      #endif

    } else if (  (type[is] == REFLECTIVE) 
              || (type[is] == AXISYMMETRIC)
              || (type[is] == EQTSYMMETRIC)){ 

      if (side[is] == X1_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(2*ibeg - 1 - i,j,k);
      }else if (side[is] == X1_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(2*iend + 1 - i,j,k);
      }else if (side[is] == X2_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,2*jbeg - 1 - j,k);
      }else if (side[is] == X2_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,2*jend + 1 - j,k);
      }else if (side[is] == X3_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j,2*kbeg - 1 - k);
      }else if (side[is] == X3_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j,2*kend + 1 - k);
      }
  
    } else if (type[is] == PERIODIC && !par_dim[is/2]) {

      if (side[is] == X1_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i + nx1,j,k);
      }else if (side[is] == X1_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i - nx1,j,k);
      }else if (side[is] == X2_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j + nx2,k);
      }else if (side[is] == X2_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j - nx2,k);
      }else if (side[is] == X3_BEG){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j,k + nx3);
      }else if (side[is] == X3_END){
        #pragma acc parallel loop collapse(3) present(fl, box)
        BOX_LOOP(box,i,j,k) fl(i,j,k) = fl(i,j,k - nx3);
      }

    } else if (type[is] == USERDEF){
      // Any Idea ?? 
    } else {
      // printLog ("! FlagBoundary(): unknown boundary condition \n");
      // printLog ("                  type[%d] = %d\n",is,type[is]);
    }

  }  // End Loop on sides 

  #pragma acc exit data delete(box)
}

