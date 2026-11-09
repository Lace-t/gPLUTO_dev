////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief 

  \author M. Rossazza (marco.rossazza@unito.it) 
          A. Mignone  (andrea.mignone@unito.it)
  \date   Nov 17, 2025
*/
////////////////////////////////////////////////////////////////////////

#ifndef MP5_ALPHA
  #define MP5_ALPHA 4.0
#endif

////////////////////////////////////////////////////////////////////////
/* Inline functions                                                   */
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/*!
  Monotonicity-preserving reconstruction of Suresh-Huynh, JCP (1997) 136.
  
  \param [in]  F    a 1D array to be reconstructed
  \param [in] sign  an integer with value +1/-1 giving the orientation
  
  \return On output it returns the right-edge value.
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline double MP5_States(double *F, int sign)
{
  double f, d2, d2p, d2m;
  double dMMm, dMMp;
  double scrh1, scrh2, fmin, fmax;
  double fAV, fMD, fLC, fUL, fMP;

  const double alpha = MP5_ALPHA;
  const double epsm = 1.e-12;

  int ipp = 0+2*sign;
  int ip  = 0+1*sign;
  int im  = 0-1*sign;
  int imm = 0-2*sign;

  double dFp = F[ip] - F[0];
  double dFm = F[0]  - F[im];

  // Reconstruction from punctual or average values 

  #ifdef HIGH_ORDER
  f  = 3.0*F[imm] - 20.0*F[im] + 90.0*F[0] + 60.0*F[ip] - 5.0*F[ipp];
  f /= 128.0;
  #else
  f  = 2.0*F[imm] - 13.0*F[im] + 47.0*F[0] + 27.0*F[ip] - 3.0*F[ipp]; // Eq. (2.1) 
  f /= 60.0;
  #endif

  fMP = F[0] + MINMOD_LIMITER(dFp, alpha*dFm);  // -- Eq. (2.12) -- 

  if ((f - F[0])*(f - fMP) <= epsm) return f;

  d2m = F[imm] + F[0] - 2.0*F[im];    // -- Eq. (2.19) -- 
  d2  = F[im]  + F[ip ] - 2.0*F[0];
  d2p = F[0] + F[ipp] - 2.0*F[ip];    // -- Eq. (2.19) -- 

  scrh1 = MINMOD_LIMITER(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD_LIMITER(d2, d2p);
  dMMp  = MINMOD_LIMITER(scrh1,scrh2);   // -- Eq. (2.27) -- 

  scrh1 = MINMOD_LIMITER(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD_LIMITER(d2, d2m);
  dMMm  = MINMOD_LIMITER(scrh1,scrh2);  // -- Eq. (2.27) -- 

  fUL = F[0] + alpha*dFm;              // -- Eq. (2.8) -- 
  fAV = 0.5*(F[0] + F[ip]);            // -- Eq. (2.16) -- 
  fMD = fAV - 0.5*dMMp;                // -- Eq. (2.28) -- 
  fLC = 0.5*(3.0*F[0] - F[im]) + 4.0/3.0*dMMm;  // -- Eq. (2.29) -- 

  scrh1 = MIN(F[0], F[ip]); scrh1 = MIN(scrh1, fMD);
  scrh2 = MIN(F[0], fUL);   scrh2 = MIN(scrh2, fLC);
  fmin  = MAX(scrh1, scrh2);  // -- Eq. (2.24a) -- 

  scrh1 = MAX(F[0], F[ip]); scrh1 = MAX(scrh1, fMD);
  scrh2 = MAX(F[0], fUL);   scrh2 = MAX(scrh2, fLC);
  fmax  = MIN(scrh1, scrh2);  // -- Eq. (2.24b) -- 

//  f = Median(f, fmin, fmax);            /* -- Eq. (2.26) -- */
  f = f + MINMOD_LIMITER(fmin-f,fmax-f);  // -- Eq. (2.26) -- 

  return f;
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*!
  Monotonicity-preserving reconstruction of Suresh-Huynh, JCP (1997) 136.
  
  \param [in]  F   a 1D array to be reconstructed
  \param [in]  j   the point at which reconstruction is needed.
  
  \return On output it returns the right-edge value.
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline double WENOZ_States(double *F, int dir)
{
  int jmm = 0 - 2*dir;
  int jm  = 0 - 1*dir;
  int jp  = 0 + 1*dir;
  int jpp = 0 + 2*dir;

  double a0, a1, a2;
  double f0, f1, f2;
  double w0, w1, w2;
  double sum_a, f;
  double thirteen_12 = 13.0/12.0;
  double b0, b1, b2, t5;

  a0 = F[jmm] - 2.0*F[jm] +     F[0];
  a1 = F[jmm] - 4.0*F[jm] + 3.0*F[0];
  b0 = thirteen_12*a0*a0 + 0.25*a1*a1;

  a0 = F[jm] - 2.0*F[0] + F[jp];
  a1 = F[jm] - F[jp];
  b1 = thirteen_12*a0*a0 + 0.25*a1*a1;

  a0 =     F[0] - 2.0*F[jp] + F[jpp];
  a1 = 3.0*F[0] - 4.0*F[jp] + F[jpp];
  b2 = thirteen_12*a0*a0 + 0.25*a1*a1;

  t5 = fabs(b0-b2);

  #ifdef HIGH_ORDER
  a0 =  1.0*(1.0 + t5/(b0 + 1.e-40));
  a1 = 10.0*(1.0 + t5/(b1 + 1.e-40));
  a2 =  5.0*(1.0 + t5/(b2 + 1.e-40));
  #else
  a0 = 1.0*(1.0 + t5/(b0 + 1.e-40));
  a1 = 6.0*(1.0 + t5/(b1 + 1.e-40));
  a2 = 3.0*(1.0 + t5/(b2 + 1.e-40));
  #endif

  sum_a = 1.0/(a0 + a1 + a2);
  w0 = a0*sum_a;
  w1 = a1*sum_a;
  w2 = a2*sum_a;

  #ifdef HIGH_ORDER
  f0 = (3.0*F[jmm] - 10.0*F[jm] + 15.0*F[0]);
  f1 = (   -F[jm]  + 6.0*F[0]   +  3.0*F[jp]);
  f2 = (3.0*F[0]   + 6.0*F[jp]  -      F[jpp]);
  f = (w0*f0 + w1*f1 + w2*f2)/8.0;
  #else
  f0 = (2.0*F[jmm] - 7.0*F[jm] + 11.0*F[0]);
  f1 = (   -F[jm]  + 5.0*F[0]  +  2.0*F[jp]);
  f2 = (2.0*F[0]   + 5.0*F[jp] -      F[jpp]);
  f = (w0*f0 + w1*f1 + w2*f2)/6.0;
  #endif

  return(f);
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*!  
  Provide interface values F_{i+1/2} using 3rd-order WENO.
  We implement two formulations: 
  
  - "Classical" [Jiang \& Shu, JCP 126 (1996) 202]
  - "Improved"  [Yamaleev and Carpenter, JCP 228 (2009) 3025] 
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline double WENO3_States(double *F, int dir, double dx)
{
  int jp = 0 + dir;
  int jm = 0 - dir;
  double a0, b0, w0, f0, a1, b1, w1, f1;
  double tau, dx2;
  // const double eps = 1.e-6;

  dx2 = dx*dx;

  b0 = F[jp] - F[0];
  b1 = F[jm] - F[0];

  b0 = b0*b0;
  b1 = b1*b1;

/* -- "Classical" version -- */
/*
  a0 = eps + b0;
  a1 = eps + b1;
  
  #ifdef HIGH_ORDER
  a0 = 3.0/(a0*a0);
  a1 = 1.0/(a1*a1);
  #else
  a0 = 2.0/(a0*a0);
  a1 = 1.0/(a1*a1);
  #endif
*/
/* -- "Improved" version -- */

  tau = (F[jp] - 2.0*F[0] + F[jm]);
  tau = tau*tau;

  #ifdef HIGH_ORDER
  a0 = 3.0*(1.0 + tau/(dx2 + b0));
  a1 = 1.0*(1.0 + tau/(dx2 + b1));
  #else
  a0 = 2.0*(1.0 + tau/(dx2 + b0));
  a1 = 1.0*(1.0 + tau/(dx2 + b1));
  #endif

// -- Complete the reconstruction -- 

  w0 = a0/(a0 + a1);
  w1 = a1/(a0 + a1);

  f0 =  0.5*( F[jp] +     F[0]);
  f1 =  0.5*(-F[jm] + 3.0*F[0]);
  return (w0*f0 + w1*f1);
}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*!
  Implement the 3rd-order limiter function, 
  Eq. [3.22]
*/
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
#pragma acc routine seq
inline double LimO3_States (double *v, int dir, double dx)
{
  double r = 0.1;
  double a,b,c,q, q1, lim;
  double eta, psi, eps = 1.e-12;
  int ip = 0 + dir;
  int im = 0 - dir;

  double dvp = v[ip] - v[0];
  double dvm = v[0] - v[im];

  double sp   = (dvp >= 0.0 ? 1.0:-1.0); 
  double advp = fabs(dvp);
  double dvm1 = sp*dvm; 

  q  = (2.0*dvp + dvm)/3.0;
  q1 = (2.0*fabs(dvp) + sp*dvm)/3.0;
  
  a = MIN(1.5*advp, 2.0*dvm1);
  a = MIN(q1,a);
  b = MAX(-0.5*dvm1,a);
  c = MIN( q1,b);
  psi = sp*MAX(0.0,c);

  eta = r*dx;
  eta = (dvm*dvm + dvp*dvp)/(eta*eta);
  if ( eta <= 1.0 - eps) {
    lim = q;
  }else if (eta >= 1.0 + eps){
    lim = psi;
  }else{
    psi =   (1.0 - (eta - 1.0)/eps)*q
          + (1.0 + (eta - 1.0)/eps)*psi;
    lim = 0.5*psi;
  }

  return v[0] + 0.5*lim;
}
#pragma omp end declare target




////////////////////////////////////////////////////////////////////////
/*!
    Implement PPM reconstruction
  
   ///////////////////////////////////////////////////////////////////// */
#pragma omp declare target
#pragma acc routine seq
inline void PPM_States(double *v, double &vm, double &vp)
{
  const int i = 0;  
  double dv, dvp, dvm;
  const double cm = 2.0, cp = 2.0;

#if 1
  // 1. Initialize the left and right states to unique interface value.

  vp =  (-v[i-1] + 7.0*(v[i]   + v[i+1]) - v[i+2])/12.0;
  vm =  (-v[i-2] + 7.0*(v[i-1] + v[i])   - v[i+1])/12.0;

  // 2. Constraint values to lie between adjacent cell averages.
  //     Same as imposing 
  //     MIN[v(i),v(i+1)] < vp < MAX[v(i), v(i+1)]

  dv  = v[i+1] - v[i];
  dvp = vp - v[i];
  vp  = v[i] + MINMOD_LIMITER(dvp, dv); 

  dv  = v[i] - v[i-1];
  dvm = vm - v[i-1];
  vm  = v[i-1] + MINMOD_LIMITER(dvm, dv);

  // 3. Apply parabolic limiter

  dvp = vp - v[i];
  dvm = vm - v[i];

  if (dvp*dvm >= 0.0) dvp = dvm = 0.0;
  else{
    if      (fabs(dvp) >= cm*fabs(dvm)) dvp = -cm*dvm;
    else if (fabs(dvm) >= cp*fabs(dvp)) dvm = -cp*dvp;
  }

  // 4. Construct left and right states

  vp = v[i] + dvp;
  vm = v[i] + dvm;
#else

  // 1. Initialize the left and right states to unique interface value.

  vp =  (-v[i-1] + 7.0*(v[i]   + v[i+1]) - v[i+2])/12.0;
  vm =  (-v[i-2] + 7.0*(v[i-1] + v[i])   - v[i+1])/12.0;

  double vmin = MIN(v[i-1], v[i+1]);
  double vmax = MAX(v[i-1], v[i+1]);
  vmin = MIN(vmin, v[i]);
  vmax = MAX(vmax, v[i]);
  if (vp > vmax) vp = vmax;
  if (vp < vmin) vp = vmin;
  if (vm > vmax) vm = vmax;
  if (vm < vmin) vm = vmin;
  

double vL = vm;   // v_{i-1/2} seen from cell i
double vR = vp;   // v_{i+1/2} seen from cell i

// Parabola curvature term
dv  = vR - vL;
double v6  = 6.0*v[i] - 3.0*(vL + vR);

// Parabolic limiting
if ( (vR - v[i])*(v[i] - vL) <= 0.0 ) {
  vL = vm = v[i];
  vR = vp = v[i];
  v6 = 0.0;
} else {
  double lim = 2.0 * MIN(fabs(vR - v[i]), fabs(v[i] - vL));
  if (fabs(v6) > lim) v6 = copysign(lim, v6);

  // (optional but common) recompute edges consistent with limited curvature
  vL = vm = v[i] - 0.5*dv - (1.0/6.0)*v6;
  vR = vp = v[i] + 0.5*dv - (1.0/6.0)*v6;
}

#endif

// dvp = v[i+1] - v[i];
// dvm = v[i] - v[i-1];
// vp = v[i] + 0.5*MC_LIMITER(dvp, dvm);
// vm = v[i] - 0.5*MC_LIMITER(dvp, dvm);

}
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*!
    Implement PLM reconstruction
  
   ///////////////////////////////////////////////////////////////////// */
#pragma omp declare target
#pragma acc routine seq
inline void PLM_States(double *v, double &vm, double &vp)
{
  const int i = 0;  

  double dvp = v[i+1] - v[i];
  double dvm = v[i] - v[i-1];
  double dv  = MC_LIMITER(dvp, dvm);

  vp = v[i] + 0.5*dv;
  vm = v[i] - 0.5*dv;

}
#pragma omp end declare target