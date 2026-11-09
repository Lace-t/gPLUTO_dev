////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Collect Internal Boundary routines.


  \author A. Mignone (andrea.mignone@unito.it)
  \date   Sep 28, 2025

  \b References
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)

#define IB_SET(flag) (flag & FLAG_INTERNAL_BOUNDARY)

// Explicit template instantiations
template void InternalBoundaryReflect<IDIR> (Data *, Grid *, RBox *);
template void InternalBoundaryReflect<JDIR> (Data *, Grid *, RBox *);
template void InternalBoundaryReflect<KDIR> (Data *, Grid *, RBox *);

// ********************************************************************* 
/*! 
 *  Reflect local 1D states across an internal boundary.
 * 
 *  This function should be called prior to any reconstruction 
 *  method in order to correctly symmetrize or anti-symmetrize 
 *  variables (\c sweep->v[]) across a fluid-solid interface boundary.
 *  Zones immediately to the left or to the right of an interface 
 *  are reconstructed in a symmetric way for all variables except 
 *  for the normal velocity which is anti-symmetric.
 *
 *********************************************************************** */
template<int dir>
void InternalBoundaryReflect (Data *data, Grid *grid, RBox *box)
{
  int i, j, k, nv;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int  ni, nj;
  int ngh = grid->nghost[dir];

  constexpr int VXn = VX1 + dir;

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);

  ARRAY_STRIDES (grid, dir, ni, nj);

  // ------------------------------------------------------
  // 1. Loop through box to select 
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(2) present(data, grid)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*offset;

    Ary2D v(&data->sweep->vC[offset1],ni,NVAR);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    #pragma acc loop vector
    for (i = ibeg; i <= iend; i++){
      
      // ------------------------------------------
      // 1a. Fluid (i) | Solid (i+1) interface.
      // Loop on "solid" cells, symmetrize profiles
      // -------------------------------------------
      
      if (!IB_SET(flag(i)) && IB_SET(flag(i+1))){  // Fluid (i) |  Solid (i+1) 
        int sbeg = i+1;
        int send = sbeg + ngh - 1;
        int npoints = grid->nint[dir];

        send = MIN(send, npoints);
        #pragma acc loop seq
        for (int i1 = sbeg; i1 <= send; i1++){
          if (IB_SET(flag(i1))){
            NVAR_LOOP(nv) v(i1,nv) = v(2*sbeg-i1-1,nv);
            v(i1,VXn) *= -1.0;
          }
        } 
      }

      // ------------------------------------------
      // 1b. Solid (i) | Fluid (i+1) interface.
      // Loop on solid cells, symmetrize profiles.
      // ------------------------------------------

      if (dir != IDIR && IB_SET(flag(i)) && !IB_SET(flag(i+1))){   // Solid (i) |  Fluid (i+1) 
        int sbeg = i;
        int send = sbeg - ngh + 1;
        send = MAX(sbeg, 0);
        #pragma acc loop seq
        for (int i1 = sbeg; i1 >= send; i1--){
          if (IB_SET(flag(i1))){
            NVAR_LOOP(nv) v(i1,nv) = v(2*sbeg-i1+1,nv);
            v(i1,VXn) *= -1.0;
          }
        } 
      }
    }
  }}
}      
#endif // (INTERNAL_BOUNDARY == YES) && (INTERNAL_BOUNDARY_REFLECT == YES)
