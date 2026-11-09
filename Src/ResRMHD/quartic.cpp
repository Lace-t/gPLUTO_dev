////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Solve quartic and cubic equations.



#####################################################################
  THIS FILE HAS BEEN DISMISSED ON October 4 2016.
  RELEVANT FUNCTIONS INCORPORATED INSIDE math_root_finders.c
#####################################################################



  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 4, 2016
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"


// #pragma acc routine seq
void QuarticPrintCoeffs(double b, double c, double d, double e);
// #pragma acc routine seq
int CubicSolve2 (double b, double c, double d, double z[]);
#pragma acc routine seq
int CubicSolve (double, double, double, double *z);


static int debug_print=0;

////////////////////////////////////////////////////////////////////////
int QuarticSolve (double b, double c, double d, double e, double *z)
/*!
  Solve a quartic equation in the form
  \f[
  z^4 + bz^3 + cz^2 + dz + e = 0
  \f]
  For its purpose, it is assumed that \b ALL roots are real.
  This makes things faster.
  
  \param [in] b   coefficient of the quartic
  \param [in] c   coefficient of the quartic
  \param [in] d   coefficient of the quartic
  \param [in] e   coefficient of the quartic
  \param [out] z  vector containing the
  (double) roots of the quartic
  
  \b Reference:
  
  http://www.1728.com/quartic2.htm  (????)
  
  \return Return 0 on success, 1 if cubic solver fails, 2 if NaN has
  been found and 3 if quartic is not satisfied.
*/
////////////////////////////////////////////////////////////////////////
{
  int  n, j, ifail;
  double b2, sq;
  double a2, a1, a0, u[4];
  double p, q, r, f;
  const double three_256 = 3.0/256.0;
  const double one_64 = 1.0/64.0;
  double sz1, sz2, sz3, sz4;

  b2 = b*b;

  // --------------------------------------------------------------
  // 1) Compute cubic coefficients using the method outlined in
  // http://eqworld.ipmnet.ru/En/solutions/ae/ae0108.pdf
  // --------------------------------------------------------------

  p = c - b2*0.375;
  q = d + b2*b*0.125 - b*c*0.5;
  r = e - 3.0*b2*b2/256.0 + b2*c/16.0 - 0.25*b*d;

  a2 = 2.0*p;
  a1 = p*p - 4.0*r;
  a0 = -q*q;

  ifail = CubicSolve(a2, a1, a0, u);
  if (ifail != 0) return 1;

  u[0] = MAX(u[0],0.0);
  u[1] = MAX(u[1],0.0);
  u[2] = MAX(u[2],0.0);


  if (u[0] != u[0] || u[1] != u[1] || u[2] != u[2]) return 1;

  sq  = -0.5*DSIGN(q);
  sz1 = sq*sqrt(u[0]);
  sz2 = sq*sqrt(u[1]);
  sz3 = sq*sqrt(u[2]);

  z[0] = -0.25*b + sz1 + sz2 + sz3;
  z[1] = -0.25*b + sz1 - sz2 - sz3;
  z[2] = -0.25*b - sz1 + sz2 - sz3;
  z[3] = -0.25*b - sz1 - sz2 + sz3;
  SortArray(z,4);
/*
  if (debug_print){
    printLog ("Quartic roots = %f  %f  %f  %f; q = %8.3e\n",z[0],z[1],z[2],z[3],q);
    CheckSolutions(b,c,d,e,z);
  }
*/
  return 0;
}

////////////////////////////////////////////////////////////////////////
int CubicSolve (double b, double c, double d, double z[])
/*!
  Solve a cubic equation in the form
  \f[
  z^3 + bz^2 + cz + d = 0
  \f]
  For its purpose, it is assumed that \b ALL roots are real.
  This makes things faster.
  
  \param [in] b coefficient of the cubic
  \param [in] c coefficient of the cubic
  \param [in] d coefficient of the cubic
  \param [out] z vector containing the roots of the cubic.
  Roots should be sorted in increasing order.
  
  \b Reference:
  - Section 5.6 of Numerical Recipe
  
  \return Return 0 on success.
*/
////////////////////////////////////////////////////////////////////////
{
  double b2;
  double Q, R;
  double sQ, arg, theta, cs, sn, p;
  const double one_3  = 1.0/3.0;
  const double one_27 = 1.0/27.0;
  const double sqrt3  = sqrt(3.0);

  b2 = b*b;

  // ----------------------------------------------
  // the expression for f should be
  // Q = c - b*b/3.0; however, to avoid negative
  // round-off making h > 0.0 or g^2/4 - h < 0.0
  // we let c --> c(1- 1.1e-16)
  // ----------------------------------------------

  Q  = b2*one_3 - c*(1.0 - 1.e-16);    // = 3*Q, with Q given by Eq. [5.6.10] 
  R  = b*(2.0*b2 - 9.0*c)*one_27 + d;  // = 2*R, with R given by Eq. [5.6.10] 

Q = MAX(Q, 0.0);
/*
if (fabs(Q) < 1.e-18){
  printLog ("CubicSolve() very small Q = %12.6e\n",Q);
  QUIT_PLUTO(1);
}
if (Q < 0.0){
  printLog ("! CubicSolve(): Q = %8.3 < 0 \n",Q);
  QUIT_PLUTO(1);
}
*/

  // -------------------------------------------------------
  // We assume the cubic *always* has 3 real root for
  // which R^2 > Q^3.
  // It follows that Q is always > 0
  // -------------------------------------------------------

  sQ  = sqrt(Q)/sqrt3;
  arg = -1.5*R/(Q*sQ);

//  this is to prevent unpleseant situation
//  where both g and i are close to zero       

  arg = MAX(-1.0, arg);
  arg = MIN( 1.0, arg);


  theta = acos(arg)*one_3;     // Eq. [5.6.11], note that  pi/3 < theta < 0 

  cs = cos(theta);              //   > 0   
  sn = sqrt3*sin(theta);    //   > 0   
  p  = -b*one_3;

  z[0] = -sQ*(cs + sn) + p;
  z[1] = -sQ*(cs - sn) + p;
  z[2] =  2.0*sQ*cs    + p;

  return(0);
}
