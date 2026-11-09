////////////////////////////////////////////////////////////////////////
/* Collect various tools.

   This files collect various tools such as memory allocation, wrapper
   to prim/cons variables and so forth.

   This code is distributed under the BSD 3-Clause License (see below).
   Last modified: 14 May 2020 by A. Mignone

  ----------------------------------------------------------------------
   BSD 3-Clause License
   Copyright (c) [2020], [Andrea Mignone]
   Copyright (c) [2020], NVIDIA CORPORATION
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  ----------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////
   
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Allocate memory for a 1-D array of any basic data type starting
  at 0.
  
  \param [in] nx    number of elements to be allocated
  \param [in] dsize data-type of the array to be allocated
  
  \return A pointer of type (char *) to the allocated memory area.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
char *Array1D (int nx, size_t dsize)
{
  char *v;
  v = (char *) malloc (nx*dsize);
  return v;
}

////////////////////////////////////////////////////////////////////////
/*!
  Free memory allocated by the pointer *v.
*/
////////////////////////////////////////////////////////////////////////
void FreeArray1D (void *v)
{
  free ((char *)v);
}

////////////////////////////////////////////////////////////////////////
/*!
  Allocate memory for a 2-D array of any basic data type.
  
  \param [in] nx    number of elements in the 2nd dimension
  \param [in] ny    number of elements in the 1st dimension
  \param [in] dsize data-type of the array to be allocated
  
  \return A pointer of type (char **) to the allocated memory area
  with index range [0...nx-1][0...ny-1]
*/
////////////////////////////////////////////////////////////////////////
char **Array2D (int nx, int ny, size_t dsize)
{
  int i;
  char **m;

  m    = (char **)malloc ((size_t) nx*sizeof(char *));
  m[0] = (char *) malloc ((size_t) nx*ny*dsize);

  for (i = 1; i < nx; i++) m[i] = m[(i - 1)] + ny*dsize;

  return m;
}

////////////////////////////////////////////////////////////////////////
/*!
  Free memory allocated by the double pointer **v.
*/
////////////////////////////////////////////////////////////////////////
void FreeArray2D (void **m)
{
  free ((char *) m[0]);
  free ((char *) m);
}

////////////////////////////////////////////////////////////////////////
/*!
  Allocate memory for a 3-D array of any basic data type.
  
  \param [in] nx    number of elements in the 3rd dimension
  \param [in] ny    number of elements in the 2nd dimension
  \param [in] nz    number of elements in the 1st dimension
  \param [in] dsize data-type of the array to be allocated
  
  \return A pointer of type (char ***) to the allocated memory area
  with index range [0...nx-1][0...ny-1][0...nz-1]
*/
////////////////////////////////////////////////////////////////////////
char ***Array3D (int nx, int ny, int nz, size_t dsize)
{
  int i, j;
  char ***m;

  m = (char ***) malloc ((size_t) nx*sizeof (char **));
  m[0] = (char **) malloc ((size_t) nx*ny*sizeof(char *));
  m[0][0] = (char *) malloc ((size_t) nx*ny*nz*dsize);

  // ---------------------------
  // single subscript: i
  // ---------------------------

  for (i = 1; i < nx; i++) m[i] = m[i - 1] + ny;

  // ---------------------------
  // double subscript:
  // 
  // (i,0)  (0,j)
  // (i,j)
  // ---------------------------

  for (j = 1; j < ny; j++) m[0][j] = m[0][j - 1] + nz*dsize;
  for (i = 1; i < nx; i++) m[i][0] = m[i - 1][0] + ny*nz*dsize;

  for (i = 1; i < nx; i++) {
  for (j = 1; j < ny; j++) {
    m[i][j] = m[i][j - 1] + nz*dsize;
  }}

  for (j = 0; j < ny; j++){
  for (i = 0; i < nx; i++){
    if (m[i][j] == NULL){
      printf ("! Allocation failure in Array3D\n");
      exit(1);
    }
  }}

  return m;
}

////////////////////////////////////////////////////////////////////////
/*!
  Free memory allocated by the pointer ***v.
*/
////////////////////////////////////////////////////////////////////////
void FreeArray3D (void ***m)
{
  free ((char *) m[0][0]);
  free ((char *) m[0]);
  free ((char *) m);
}

////////////////////////////////////////////////////////////////////////
/*!
  Allocate memory for a 4-D array of any basic data type.
  
  \param [in] nx    number of elements in the 4th dimension
  \param [in] ny    number of elements in the 3rd dimension
  \param [in] nz    number of elements in the 2nd dimension
  \param [in] nv    number of elements in the 1st dimension
  \param [in] dsize data-type of the array to be allocated
  
  \return A pointer of type (char ****) to the allocated memory area
  with index range [0...nx-1][0...ny-1][0...nz-1][0...nv-1]
*/
////////////////////////////////////////////////////////////////////////
char ****Array4D (int nx, int ny, int nz, int nv, size_t dsize)
{
  int i, j, k;
  char ****m;

  m = (char ****) malloc ((size_t) nx*sizeof (char ***));
  m[0] = (char ***) malloc ((size_t) nx*ny*sizeof (char **));
  m[0][0] = (char **) malloc ((size_t) nx*ny*nz*sizeof (char *));
  m[0][0][0] = (char *) malloc ((size_t) nx*ny*nz*nv*dsize);

  // ---------------------------
  // single subscript: i
  // ---------------------------

  for (i = 1; i < nx; i++) m[i] = m[i - 1] + ny;

  // ---------------------------
  // double subscript:
  // 
  // (i,0)  (0,j)
  // (i,j)
  // ---------------------------

  for (i = 1; i < nx; i++) {
    m[i][0] = m[i - 1][0] + ny*nz;
  }
  for (j = 1; j < ny; j++) {
    m[0][j] = m[0][j - 1] + nz;
  }

  for (i = 1; i < nx; i++) {
  for (j = 1; j < ny; j++) {
    m[i][j] = m[i][j - 1] + nz;
  }}

  // ---------------------------
  // triple subscript:
  // 
  // (i,0,0) (0,j,0) (0,0,k)
  // (i,j,0) (i,0,k) (0,j,k)
  // (i,j,k)
  // ---------------------------

  for (i = 1; i < nx; i++) {
    m[i][0][0] = m[i - 1][0][0] + ny*nz*nv*dsize;
  }
  for (j = 1; j < ny; j++) {
    m[0][j][0] = m[0][j - 1][0] + nz*nv*dsize;
  }

  for (k = 1; k < nz; k++) {
    m[0][0][k] = m[0][0][k - 1] + nv*dsize;
  }


  for (i = 1; i < nx; i++) {
  for (j = 1; j < ny; j++) {
    m[i][j][0] = m[i][j - 1][0] + nz*nv*dsize;
  }}

  for (i = 1; i < nx; i++) {
  for (k = 1; k < nz; k++) {
    m[i][0][k] = m[i][0][k - 1] + nv*dsize;
  }}

  for (j = 1; j < ny; j++) {
  for (k = 1; k < nz; k++) {
    m[0][j][k] = m[0][j][k - 1] + nv*dsize;
  }}

  for (i = 1; i < nx; i++) {
    for (j = 1; j < ny; j++) {
      for (k = 1; k < nz; k++) {
        m[i][j][k] = m[i][j][k - 1] + nv*dsize;
      }
    }
  }

  return m;
}

////////////////////////////////////////////////////////////////////////
/*!
  Free memory allocated by the pointer ****v.
*/
////////////////////////////////////////////////////////////////////////
void FreeArray4D (void ****m)
{
  free ((char *) m[0][0][0]);
  free ((char *) m[0][0]);
  free ((char *) m[0]);
  free ((char *) m);
}
