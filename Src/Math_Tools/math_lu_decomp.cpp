////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Functions for LU decomposition and matrix inversion.

  \author A. Mignone (andrea.mignone@unito.it)

  \date   March 7, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include <utility>

////////////////////////////////////////////////////////////////////////
/*!
  Perform LU decomposition. Adapted from Numerical Recipes.
  On output, replaces the matrix a[0..n-1][0..n-1] with its LU
  decomposition. indx[0..n-1] is a row vector that records the row
  permutation while d = 1 (-1) if the number of raw interchanges
  was even (odd).
  This function should be used in combination with LUBackSubst() to
  solve linear system of equations.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
int LUDecompose (double *J, int n, int *indx, double *d)
#define TINY 1.0e-20
{
  int i, imax, j, k;
  int max_size = 0;
  double big, dum, sum, temp;
  double vv[(NIONS+1)];
  double a[(NIONS+1)][(NIONS+1)];
  n= (NIONS + 1);

  for (i = 0; i < n; i++){
  for (j = 0; j < n; j++){
    a[i][j] = *((J+i*n) + j);
  }}

  if (n > max_size) {
    max_size = n;
  }

  *d = 1.0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs(a[i][j])) > big)
        big = temp;
    if (big == 0.0) {
//      printLog ("! Singular matrix in routine LUDecompose - (i=%d, j=%d)",i,j); 
      return (0);
    }
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = a[i][j];
      for (k = 0; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = a[i][j];
      for (k = 0; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs (sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = TINY;
    if (j != n - 1) {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i < n; i++)
        a[i][j] *= dum;
    }
  }

  for (i = 0; i < n; i++){
  for (j = 0; j < n; j++){
    *((J+i*n) + j) = a[i][j];
  }}

  return (1); // -- success -- 
}
#undef TINY

////////////////////////////////////////////////////////////////////////
/*!
  Solve a linear system of the type <tt>Ax = b</tt>, where \c A is
  a square matrix, \c x is the array of unknowns and \c b is given.
  This function must be called after LUDecompose()
  (adapted from Numerical Recipes).
  For the 1st time, use:
  \code
  LUDecompose (A,n,indx,&d);
  LUBackSubst (A,n,indx,b);
  \endcode
  The answer \c x will be given back in \c b and the original matrix \c A
  will be destroyed.
  For all subsequent time, to solve the same system wth a different
  right-hand side \c b, use
  \code
  LUBackSubst(A,n,indx,b);
  \endcode
  where \c A is the matrix that has been decomposed by LUDecompose() and
  \em not the original matrix \c A.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void LUBackSubst (double *J, int n, int *indx, double b[])
{
  int i, ii = 0, ip, j;
  double sum;
  double a[3][3];

  for (i = 0; i < n; i++){
  for (j = 0; j < n; j++){
    a[i][j] = *((J+i*n) + j);
  }}

  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii - 1; j <= i - 1; j++)
        sum -= a[i][j] * b[j];
    else if (sum)
      ii = i + 1;
    b[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = b[i];
    for (j = i + 1; j < n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }

  for (i = 0; i < n; i++){
  for (j = 0; j < n; j++){
    *((J+i*n) + j) = a[i][j];
  }}
}

////////////////////////////////////////////////////////////////////////
/*!
  Solve a system Ax = b by Gaussian elimination (no LU decomp)
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
void GaussianElimination (double A[][3], double *x, double *b, int n)
{
  int i, j, k, imax;
  double Amax, tmp;

  for (k = 0; k < n-1; k++){        // Loop over the Gk's

    // Apply Partial pivoting here

    Amax = fabs(A[k][k]);
    imax = k;
    for (i = k+1; i < n; i++) {
      tmp = fabs(A[i][k]);
      if (tmp > Amax) {
        imax = i;
        Amax = tmp;
      }
    }
    // if (k != imax) SwapRows (A, b, k, imax, n);
    if (k != imax) {
      for (j = 0; j < n; j++)
        std::swap(A[k][j],A[imax][j]);
      std::swap(b[k],b[imax]);
    }

    // Loop over rows

    for (i = k+1; i < n; i++){
      tmp = A[i][k]/A[k][k];
      for (j = k+1; j < n; j++) A[i][j] -= tmp*A[k][j];
      A[i][k] = 0.0;
      b[i]   -= tmp*b[k];
    }
  }

  // Backsubstitution

  for (i = n-1; i >= 0; i--){
    tmp = b[i];
    for (j = n-1; j > i; j--) tmp -= x[j]*A[i][j];
    x[i] = tmp/A[i][i];
  }
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*!
  Swap rows i1 and i2 in the matrix A[][] and vector b[]
  A[][] is assumed to be a square matrix of dimension n
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
#pragma omp declare target
void SwapRows (double A[][3], double *b, int i1, int i2, int n)
{
  int j;
  for (j = 0; j < n; j++){
    std::swap(A[i1][j],A[i2][j]);
  }

  std::swap(b[i1],b[i2]);
}
#pragma omp end declare target
