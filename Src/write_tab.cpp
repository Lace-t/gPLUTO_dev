////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Write tabulated 1D or 2D ascii data files.

  WriteTabArray() provides ascii output for 2-D and 1-D data 
  in tabulated multi-column format:
  \verbatim
      .     .         .         .             .
      .     .         .         .             .
      .     .         .         .             .
    x1(i) x2(j)   var_1(i,j)  var_2(i,j)   var_3(i,j) ...
      .     .         .         .             .
      .     .         .         .             .
      .     .         .         .             .
  \endverbatim
 
  Blocks with different x2 coordinates are separated by a blank row
  (in more than 1D).
  One file (with all variables) per output is written to disk. 
  This format does not work in parallel mode and can be used for simple
  serial runs.
 
  \author A. Mignone (andrea.mignone@unito.it)
  \date   May 29, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void WriteTabArray (Output *output, char *filename, Grid *grid)
/*!
  Write tabulated array.
  
  \param [in] output    a pointer to the output structure 
                        corresponding to the tab format
  \param [in] filename  the output file name
  \param [in] grid      pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
#ifdef PARALLEL
  print ("! WriteTabArray: tab output not supported in parallel\n");
  return;
#else
  #if DIMENSIONS < 3
  int  nv, i, j, k;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  Ary4D *Vary4D = output->Vary4D;
  RBox box;

  RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, &box);

  // 1. Dump arrays in ascii format to disk 

  FILE *fout;
  fout = fopen (filename, "w");
  k = 0;
  for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){
      fprintf (fout, "%f %f ", grid->xC[IDIR][i], grid->xC[JDIR][j]);
      for (nv = 0; nv < output->nvar; nv++) {
        if (output->dump_var[nv]) {
          fprintf (fout, "%12.6e ", (*Vary4D)(i,j,k,nv));
        }
      }
      fprintf (fout, "\n");   // newline 
    }
    #if DIMENSIONS > 1
    fprintf (fout, "\n");   // skip one more empty line in 2D 
    #endif
  }
  fclose (fout);
  return;
  #else
  print ("! WriteTabArray: tab output not supported in 3D\n");
  return;
  #endif // DIMENSIONS < 3
#endif // PARALLEL

}
