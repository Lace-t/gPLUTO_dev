////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Collects various functions to operate on the RBox structure.

  - RBoxDefine() is used to set the box extent in terms of six
    indices (beg..end for each direction) and the variable location
    inside the cell.

  - RBoxSeepSet() is used to set normal, tangent and binormal indices
    with respect to the specified (sweeping) direction.
    The convention is here respects the fastest running index,
    <tt> (i,j,k) -> (j,i,k) -> (k,i,j) </tt>.
    Useful for sweeping along different directions during the time stepping
    routines.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   May 16, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Define a box given the beg / end indices in each direction
  
  \param [in]  ib    leftmost  index in the X1 direction
  \param [in]  ie    rightmost index in the X1 direction
  \param [in]  jb    leftmost  index in the X2 direction
  \param [in]  je    rightmost index in the X2 direction
  \param [in]  kb    leftmost  index in the X3 direction
  \param [in]  ke    rightmost index in the X3 direction
  \param [out] box   pointer to a RBox structure
*/
////////////////////////////////////////////////////////////////////////
void RBoxDefine(int ib, int ie, int jb, int je, int kb, int ke, RBox *box)
{
  box->ibeg = ib;
  box->iend = ie;

  box->jbeg = (INCLUDE_JDIR ? jb:0);
  box->jend = (INCLUDE_JDIR ? je:0);

  box->kbeg = (INCLUDE_KDIR ? kb:0);
  box->kend = (INCLUDE_KDIR ? ke:0);

  box->nx1 = box->iend - box->ibeg + 1;
  box->nx2 = box->jend - box->jbeg + 1;
  box->nx3 = box->kend - box->kbeg + 1;
}

////////////////////////////////////////////////////////////////////////
/*!
  Define a box given the beg / end indices in each direction (array version)
  
  \param [in]  beg   array of beg indexes
  \param [in]  end   array of end indexes
  \param [out] box   pointer to a RBox structure
*/
////////////////////////////////////////////////////////////////////////
void RBoxDefine(int *beg, int *end, RBox *box)
{
  box->ibeg = beg[IDIR];
  box->iend = end[IDIR];

  box->jbeg = (INCLUDE_JDIR ? beg[JDIR]:0);
  box->jend = (INCLUDE_JDIR ? end[JDIR]:0);

  box->kbeg = (INCLUDE_KDIR ? beg[KDIR]:0);
  box->kend = (INCLUDE_KDIR ? end[KDIR]:0);

  box->nx1 = box->iend - box->ibeg + 1;
  box->nx2 = box->jend - box->jbeg + 1;
  box->nx3 = box->kend - box->kbeg + 1;
}

////////////////////////////////////////////////////////////////////////
/*!
  Return beg and end indices as arrays of ints
  
  \param [in]   box     pointer to input box
  \param [out]  beg     an array of 3 indices specifying the beginning
  \param [out]  end     an array of 3 indices specifying the end
*/
////////////////////////////////////////////////////////////////////////
void RBoxGetIndexes (RBox *box, int *beg, int *end)
{
  beg[IDIR] = box->ibeg;
  beg[JDIR] = box->jbeg;
  beg[KDIR] = box->kbeg;

  end[IDIR] = box->iend;
  end[JDIR] = box->jend;
  end[KDIR] = box->kend;
}

////////////////////////////////////////////////////////////////////////
/*!
  Return the total box size.
  
  \param [in]   box     pointer to input box
*/
////////////////////////////////////////////////////////////////////////
long int RBoxGetSize (RBox *box)
{
  long int nx = (box->iend - box->ibeg) + 1;
  long int ny = (box->jend - box->jbeg) + 1;
  long int nz = (box->kend - box->kbeg) + 1;

  return nx*ny*nz;
}

////////////////////////////////////////////////////////////////////////
/*!
  Shift the beg, end indices of an RBox
  
  \param [in]   box     pointer to input box
  \param [in]   offset  an array of 3 indices specifying the offset
*/
////////////////////////////////////////////////////////////////////////
void RBoxOffset (RBox *box, int *offset)
{
  box->ibeg += offset[IDIR];
  box->iend += offset[IDIR];

  box->jbeg += offset[JDIR];
  box->jend += offset[JDIR];

  box->kbeg += offset[KDIR];
  box->kend += offset[KDIR];
}

////////////////////////////////////////////////////////////////////////
/*!
  Display info about a RBox
*/
////////////////////////////////////////////////////////////////////////
void RBoxShow(RBox *box)
{
  printLog ("[%d, %d] [%d, %d] [%d, %d]\n",
             box->ibeg, box->iend,box->jbeg, box->jend,box->kbeg, box->kend);
}

// Explicit template instantiations
template void RBoxPermuteIndices<IDIR>(RBox *, int &, int &, 
                                               int &, int &, 
                                               int &, int &);
template void RBoxPermuteIndices<JDIR>(RBox *, int &, int &, 
                                               int &, int &, 
                                               int &, int &);
template void RBoxPermuteIndices<KDIR>(RBox *, int &, int &, 
                                               int &, int &, 
                                               int &, int &);

////////////////////////////////////////////////////////////////////////
/*!
  \brief Permutes the index ranges of an RBox based on the given direction.
  
  This function reorders the `ibeg`, `iend`, `jbeg`, `jend`, `kbeg`, and `kend`
  values according to the specified direction template parameter.
  
  \tparam dir The direction along which to permute. Must be one of IDIR, JDIR, or KDIR.
  \param box  Pointer to the RBox object whose indices will be permuted.
  \param ibeg Start index for the first axis (updated by the function).
  \param iend End index for the first axis (updated by the function).
  \param jbeg Start index for the second axis (updated by the function).
  \param jend End index for the second axis (updated by the function).
  \param kbeg Start index for the third axis (updated by the function).
  \param kend End index for the third axis (updated by the function).
  
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void RBoxPermuteIndices(RBox *box, int &ibeg, int &iend, 
                                   int &jbeg, int &jend, 
                                   int &kbeg, int &kend)
{
  if constexpr (dir == IDIR) {
    ibeg = box->ibeg; iend = box->iend;
    jbeg = box->jbeg; jend = box->jend;
    kbeg = box->kbeg; kend = box->kend;
  } else if constexpr (dir == JDIR) {
    ibeg = box->jbeg; iend = box->jend;
    jbeg = box->ibeg; jend = box->iend;
    kbeg = box->kbeg; kend = box->kend;
  } else if constexpr (dir == KDIR) {
    ibeg = box->kbeg; iend = box->kend;
    jbeg = box->ibeg; jend = box->iend;
    kbeg = box->jbeg; kend = box->jend;
  } else {
    string msg = "! RBoxPermuteIndices(): invalid dir = " + to_string(dir);
    QuitPLUTO(msg, 1);
  }
}
