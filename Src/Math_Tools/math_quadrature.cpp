#include "gpluto.hpp"

#pragma acc routine seq
double h1h2h3 (double, double, double);
#pragma acc routine seq
double   h2h3 (double, double, double);
#pragma acc routine seq
double   h1h3 (double, double, double);
#pragma acc routine seq
double   h1h2 (double, double, double);
#pragma acc routine seq
double   h1 (double, double, double);
#pragma acc routine seq
double   h2 (double, double, double);
#pragma acc routine seq
double   h3 (double, double, double);

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
double InitGaussAverage(intList *nv_list, 
                        double xb, double xe, 
                        double yb, double ye, 
                        double zb, double ze, double *uav,
                        int Ng0, int Ng1, int Ng2, RunConfig &run_config)
/*!
  Compute volume- or surface-averaged conserved variables using 
  Gaussian quadrature, starting from the profiles of primitive variables 
  defined in the Init() function.
  The general formula is:
  
  \f[
  <Q> = \frac{\int q(x_1, x_2, x_3) J dx1\,dx2\,dx3}
  {\int J dx1\,dx2\,dx3} 
  \f]
  
  where J(x1,x2,x3) contains one or more scale factors.
  
  When all directions have non-zero integration intervals 
  [xb != xe, yb != ye, zb != ze or, in compact notation (d1, d2, d3)] 
  then we set J = h1*h2*h3 and a volume integral is obtained.
  
  However, when the interval over a given direction is zero, an area 
  integral is performed and the corresponding scale factor is removed.
  In polar coordinates, for instance,
  
  (d1, d2, d3) --> Volume: J() = h1*h2*h3 = R 
  (    d2, d3) --> AR:     J() = h2*h3    = R   (when xb = xe)  
  (d1,   , d3) --> Aphi:   J() = h1*h3    = 1   (when yb = ye)   
  (d1, d2,   ) --> Az:     J() = h1*h2    = R   (when zb = ze)  
  
  Likewise, in spherical coordinates:
  
  (d1, d2, d3) --> Volume:  J() = h1*h2*h3 = r*r*s
  (    d2, d3) --> Ar:      J() = h2*h3    = r*r*s
  (d1,   , d3) --> Ath:     J() = h1*h3    = r*s
  (d1, d2,   ) --> Aphi:    J() = h1*h2    = r
  
  
  \param [in]   nv    useless at the moment
  \param [in]   xb,xe  integration limit for the x1-direction
  \param [in]   yb,ye  integration limit for the x2-direction
  \param [in]   zb,ze  integration limit for the x3-direction
  \param [out]  uav    volume- or surface-averaged conservative array
  \param [in]   Ngx    Number of Gaussian point in the x1 direction
  \param [in]   Ngy    Number of Gaussian point in the x2 direction
  \param [in]   Ngz    Number of Gaussian point in the x3 direction
*/
////////////////////////////////////////////////////////////////////////
{

int Ngx = Ng0; int Ngy = Ng1; int Ngz = Ng2; // there's a bug passing Ngx,Ngy,Ngz directly

  int  nv;
  int  isub, jsub, ksub;
  char  d1, d2, d3;
  double x1s, x2s, x3s;
  double dx1, dx2, dx3;
  double eps = 1.e-8;   // Interval width below eps are considered zero 
  double wdim;
  double Jv, norm;
  double xg[8], wgx[8];
  double yg[8], wgy[8];
  double zg[8], wgz[8];
  double v[256], u[256];
  double (*J)(double, double, double);

  // ------------------------------------------------------
  // 0. Set domain widths
  // ------------------------------------------------------

  dx1 = (xe - xb); d1 = (fabs(dx1) > eps) ? 1:0;
  dx2 = (ye - yb); d2 = (fabs(dx2) > eps) ? 1:0;
  dx3 = (ze - zb); d3 = (fabs(dx3) > eps) ? 1:0;

  // ------------------------------------------------------
  // 1. Automatically detect Jacobian from domain size
  // ------------------------------------------------------

  if (d1 && d2 && d3){
    J = h1h2h3;
  }else if (d1 && d2)  {
    J   = h1h2;
    Ngz = 1;
  }else if (d1 && d3) {
    J = h1h3;
    Ngy = 1;
  }else if (d2 && d3) { 
    J = h2h3;
    Ngx = 1;
  }else if (d1){
    J = h1;
  }else if (d2){
    J = h2;
  }else if (d3){
    J = h3;
  }else{
    // printLog ("! InitGaussAverage: something's wrong with integration limits\n");
    // QUIT_PLUTO(1);
  }

  // ------------------------------------------------------
  // 2. Set Gaussian points and weights
  // ------------------------------------------------------

  GaussianWeights (xg, wgx, Ngx);
  GaussianWeights (yg, wgy, Ngy);
  GaussianWeights (zg, wgz, Ngz);

  // ------------------------------------------------------
  // 3. Integrate
  // ------------------------------------------------------

  norm = 0.0;
  if (nv_list != NULL) FOR_EACH (nv, nv_list) uav[nv] = 0.0;
  else                 NVAR_LOOP(nv) uav[nv] = 0.0;

  NVAR_LOOP(nv) u[nv] = v[nv] = 0.0;

  for (ksub = 0; ksub < Ngz; ksub++){ 
  for (jsub = 0; jsub < Ngy; jsub++){ 
  for (isub = 0; isub < Ngx; isub++){
    x1s  = 0.5*dx1*xg[isub] + 0.5*(xe + xb);
    x2s  = 0.5*dx2*yg[jsub] + 0.5*(ye + yb);
    x3s  = 0.5*dx3*zg[ksub] + 0.5*(ze + zb);
    wdim = wgx[isub]*wgy[jsub]*wgz[ksub]*0.125;

    Jv   = 1.0;
    // Jv   = J;
    
    #pragma acc routine (Init) seq
    Init (v, x1s, x2s, x3s, run_config);

    if (nv_list == NULL){ // Means all variables must be converted 
      PrimToCons (v, u, run_config);
      NVAR_LOOP(nv) uav[nv] += wdim*Jv*u[nv];
    } else {
      FOR_EACH(nv,nv_list) uav[nv] += wdim*Jv*v[nv];
    }

    norm += wdim*Jv;
  }}}
  
  if (nv_list == NULL) NVAR_LOOP(nv)         uav[nv] /= norm;
  else                 FOR_EACH(nv, nv_list) uav[nv] /= norm;

  return 0;
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h1h2h3 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == CARTESIAN
  return 1.0;
  #elif GEOMETRY == CYLINDRICAL
  return x1;
  #elif GEOMETRY == SPHERICAL
  return DIM_EXPAND(x1*x1, *sin(x2), *1.0);
  #endif
}
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h1h2 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == CARTESIAN
  return 1.0;
  #elif GEOMETRY == CYLINDRICAL
  return x1;
  #elif GEOMETRY == SPHERICAL
  return x1;
  #endif
}
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h1h3 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == SPHERICAL
  return DIM_EXPAND(x1, *sin(x2), *1.0);
  #else
  return 1.0;
  #endif
}
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h2h3 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == CARTESIAN
  return 1.0;
  #elif GEOMETRY == CYLINDRICAL
  return x1;
  #elif GEOMETRY == SPHERICAL
  return DIM_EXPAND(x1*x1, *sin(x2), *1.0);
  #endif
}


////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h1 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  return 1.0;
}
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h2 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == CARTESIAN
  return 1.0;
  #elif GEOMETRY == CYLINDRICAL
  return x1;
  #elif GEOMETRY == SPHERICAL
  return x1;
  #endif
}
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
double h3 (double x1, double x2, double x3)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  #if GEOMETRY == SPHERICAL
  return x1*sin(x2);
  #else
  return 1.0;
  #endif
}
