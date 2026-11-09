////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Grid generation using Mapped Grids.

  Provide default grid mapping and grid generation tools.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Apr 11, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if MAPPED_GRID == YES
////////////////////////////////////////////////////////////////////////
/*! 
  Exponential map
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_Exp(Grid *grid, int dir, double xi)
{
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];

  double x = xbeg*pow(xend/xbeg, xi);
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Sinusoidal grid.
  
  \param [in]   r  stretching ratio
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_Sin(Grid *grid, int dir, double xi)
{
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];
  double A = (xend - xbeg) / 2.0;
  double B = xbeg;
  double x = A * (2.0*xi + sin(2.0*M_PI*xi) / (2.0*M_PI)) + B;
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Hyperbolic Sin grid.
  
  \param [in]   xc  the clustering center
  \param [in]   a   the clustering strength
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_Sinh(Grid *grid, int dir, double xi, double a)
{
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];
  double xc   = 0.5*(xbeg + xend);
  double A = (xend - xbeg) / 2.0;
  double x = xc + A*sinh(a*(2*xi-1.0))/sinh(a);
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Single stretched grid.
  
  \param [in]   r  stretching ratio
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_Stretched(Grid *grid, int dir, double xi, double r)
{
  int    N    = grid->nint_glob[dir];
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];

  double x = xbeg + (xend-xbeg)*(pow(r,N*xi) - 1.0)/(pow(r,N) - 1.0);
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Composite grid formed by 2 patches: a linear uniform grid followed
  by a parabolically stretched one.
  
  \param[in]  q       beg of parabolic map (logical space)
  \param[in]  xq      beg of parabolic map (physical space)
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_UP(Grid *grid, int dir, double xi, double q, double xq)
{
  double x;
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];

  if (xi < q){
    x = xbeg + (xq - xbeg)*xi/q;
  }else{
    double L   = xend - xbeg;
    double den = (q-1.0)*(q-1.0);
    double a   = (q*L - (xq - xbeg))/(q*den);
    double b   = -(q*q*(2.0*xend - xbeg - xq) - (xq - xbeg))/(q*den);
    double c   =  (q*(q*xend - xbeg - xq) + xbeg)/den;
  
    x = a*xi*xi + b*xi + c;
  }
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Composite grid formed by 2 patches: a linear uniform grid followed
  by a geometrically stretched one.
  
  \param[in]  q       beg of stretched map (logical space)
  \param[in]  xq      beg of stretched map (physical space)
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_US(Grid *grid, int dir, double xi, double q, double xq)
{
  double x;
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];
  static int k[3] = {0,0,0};
  static double r_dir[3];
  double r;

  if (xi < q){
    x = xbeg + (xq - xbeg)*xi/q;
  }else{

    if (k[dir] == 0){ 
      double dx_dxi = (xq - xbeg)/q;
      double C = dx_dxi*(1.0 - q)/(xend - xq);
      r_dir[dir] = MappedGrid_StretchRatio (C);
      k[dir] = 1;
    }
    
    r = r_dir[dir];
    double esp = (xi - q)/(1.0 - q);
    if (r != 1.0) x = xq + (xend - xq)*(pow(r,esp) - 1.0)/(r - 1.0);
    else          x = xq + (xend - xq)*(xi-q)/(1.0-q);
  }
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Composite grid formed by 3 patches: a stretched one followed by a
  linear uniform grid followed by a geometrically stretched one.
  
  \param[in]  q1       beg of stretched map on the left  (logical space)
  \param[in]  xq1      beg of stretched map on the left (physical space)
  \param[in]  q2       beg of stretched map on the right (logical space)
  \param[in]  xq2      beg of stretched map on the right (physical space)
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_SUS(Grid *grid, int dir, double xi, 
                      double q1, double xq1, double q2, double xq2)
{
  double x;
  double xb = grid->xbeg_glob[dir];
  double xe = grid->xend_glob[dir];
  static int k1[3] = {0,0,0};
  static int k2[3] = {0,0,0};
  static double r1_dir[3];
  static double r2_dir[3];
  double r1, r2;

  if (xi < q1){      

    // --------------------------------------------
    // Left map (xi < q1) 
    // --------------------------------------------

    if (k1[dir] == 0){ 
      double dx_dxi = (xq2 - xq1)/(q2 - q1);
      double C = dx_dxi*q1/(xq1 - xb);
      r1_dir[dir] = MappedGrid_StretchRatio (C); 
      k1[dir]  = 1;
    }
    
    r1 = r1_dir[dir];
    double esp = (q1 - xi)/q1;
    if (r1 != 1.0) x = xq1 - (xq1 - xb)*(pow(r1,esp) - 1.0)/(r1 - 1.0);
    else           x = xb  + (xq1 - xb)*xi/q1;

  } else if (xi < q2){

    // --------------------------------------------
    // Central  map (q1 < xi < q2) 
    // --------------------------------------------

    x = xq1 + (xq2 - xq1)*(xi - q1)/(q2 - q1);

  }else{

    // --------------------------------------------
    // Right  map (q2 < xi < 1) 
    // --------------------------------------------

    if (k2[dir] == 0){ 
      double dx_dxi = (xq2 - xq1)/(q2 - q1);
      double C      = dx_dxi*(1.0 - q2)/(xe - xq2);
      r2_dir[dir] = MappedGrid_StretchRatio (C);
      k2[dir] = 1;
    }
    
    r2 = r2_dir[dir];
    double esp = (xi - q2)/(1.0 - q2);
    if (r2 != 1.0) x = xq2 + (xe - xq2)*(pow(r2,esp) - 1.0)/(r2 - 1.0);
    else           x = xq2 + (xe - xq2)*(xi - q2)/(1.0 - q2);
  }
  return x;
}

////////////////////////////////////////////////////////////////////////
/*! 
  Uniform (linear, constant spacing macp) 
*/
////////////////////////////////////////////////////////////////////////
double MappedGrid_Uniform(Grid *grid, int dir, double xi)
{
  double xbeg = grid->xbeg_glob[dir];
  double xend = grid->xend_glob[dir];
  double x = xbeg + xi*(xend - xbeg);
  return x;
}

////////////////////////////////////////////////////////////////////////
/*!
  Use a combination of Bisection / Newton-Raphson methods to 
  compute the stretch ratio by solving the equation
  
  \f[
  \frac{\log r}{1-r} = C
  \f]
*/
////////////////////////////////////////////////////////////////////////
double StretchRatioFunc (double r);

double MappedGrid_StretchRatio(double C)
{
  #define MAX_ITER 50
  int k = 0, newton_fail = 1;
  int kb = 0, kn = 0;
  double ra, rb, rm; 
  double fa, fb, fm;
  double acc  = 1.e-8;
  double dr   = 1.0;

  // -------------------------------------------------------- 
  // 0. Determine lower and upper bounds for root solver. 
  // Solution is in 0 < r < 1 when lhs > 1 and
  // r > 1     when lhs < 1
  // --------------------------------------------------------

  if (C > 1.0) {
    ra = 1.e-12; // 0.001/(C*C);
    rb = 1.0;
  }else{
    ra = 1.0;
    rb = 10.0/C;  // Empirical estimate of upper bound
  }

  fa = StretchRatioFunc(ra) - C;
  fb = StretchRatioFunc(rb) - C;

  if (fa*fb > 0.0){
    print ("! MappedGrid_StretchRatio(): no root\n");
    QuitPLUTO("",1);
  }

  if (C <= 0.0){
    print ("! MappedGrid_StretchRatio(): negative values of C not allowed\n");
    QuitPLUTO("",1);
  }

// printf ("C = %12.6e, fa = %12.6e, fb = %12.6e\n", C, fa, fb);

  // --------------------------------------------------------
  // 1. Special case: this equation has the trivial solution 
  // r = 1 when lhs = 1.0;
  // --------------------------------------------------------

  if (fabs(C-1.0) < 1.e-12) {
    rm = 1.0;
    dr = 0.0;
  }

  // --------------------------------------------------------
  // 2. Apply Newton when possible or revert to 
  // bisection algorithm 
  // --------------------------------------------------------

  while(fabs(dr) > acc){

    if (newton_fail) rm  = 0.5*(ra+rb);

    // --------------------------------------------
    // 2a. Try Newton first 
    // -------------------------------------------- 

    fm  = StretchRatioFunc(rm) - C;
    double dfm = (-log(rm)*rm + rm - 1.0)/(rm*(rm-1.0)*(rm-1.0));
    dr  = fm/dfm;
    rm -= dr;

    if (rm < ra || rm > rb){ // Newton failed
      newton_fail = 1;

      // ------------------------------------------
      // 2b. Apply Bisection 
      // ------------------------------------------ 

      rm  = 0.5*(ra+rb);
      fm  = StretchRatioFunc(rm) - C;

      dr = rb-ra;

      if (fm*fa < 0.0) {
        rb = rm;  
        fb = fm;
      }else{
        ra = rm;  
        fa = fm;
      }
      kb++;
    }else {
      newton_fail = 0;
      kn++;
    }

    k++;

    if (k == 40) {
      print ("! MappedGrid_StretchRatio(): too many iterations\n");
      QuitPLUTO("",1);
    }
  }

  print ("  MappedGrid_StretchRatio(): ");
  print ("found r = %12.6e [# iterations = %d; kNewt = %d, kBisect = %d]\n",
             rm, k, kn, kb);
  return rm;
}

double StretchRatioFunc (double r)
{
  double f;

  // Use Taylor expansion around r = 1 to avoid singularity.
  if (fabs(r-1.0) < 1.e-4) f = 1.0 - 0.5*(r-1.0) + (r-1.0)*(r-1.0)/3.0;
  else                     f = log(r)/(r - 1.0);

  return f;
}


#endif // MAPPED_GRID == YES
