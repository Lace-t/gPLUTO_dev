////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Header file for the High-Order module.

  \author A. Mignone  (andrea.mignone@unito.it)
          V. Berta    (vittoria.berta@unito.it)
  \date   Sep 30, 2024
*/
////////////////////////////////////////////////////////////////////////
#define HO_RDER_LIM            6
#define HO_JAMESON_LIM         7

/*! The HO_AVERAGE and DE_AVARAGE macros select how the averageing 
    and de-averaging are handled throughout the code.
    In particular, the following values are possible:

      0 = No Laplacian (d2Q = 0)
      1 = traditional Laplacian (3 points, unlimited) 
      2 = -
      3 = arbitrary coefficient version (5 points ) 
      4 = WENO3 parabola (3-points stencil) 
      5 = Split-limited version (5 points stencil) 
*/
#ifndef HO_AVERAGE   
  #define HO_AVERAGE  1   
#endif

#ifndef HO_DEAVERAGE   
  #define HO_DEAVERAGE  1
#endif

  // -------------------------------------------------------- 
  // The HO_FIELD_AVERAGE macro sets how central point-value 
  // and volume-average magnetic field are computed inside the
  // cell from staggered components.
  // 1 = yields the Felker & Stone 
  // 2 = new formulation  
  // --------------------------------------------------------
#ifndef HO_FIELD_AVERAGE
  #define HO_FIELD_AVERAGE    1  /**< Choose the magnetic field averaging scheme */
#endif

#ifndef HO_LAP_LIMITER
  #define HO_LAP_LIMITER     NO
#endif

#ifndef HO_ORDER_REDUCTION
  #define HO_ORDER_REDUCTION  NO
#endif

#ifndef HO_STARTUP_NGAUSS
  #define HO_STARTUP_NGAUSS 4
#endif

// #ifndef HO_WENO_LAPLACIAN
//   #define HO_WENO_LAPLACIAN   NO
// #endif

// --------------------------------------------------------
// Function prototyping
// --------------------------------------------------------

void DiscontinuityDetect(Data *, Grid *, RBox *);
void FlipSign (int, int, int *);

template<int dir>
void HO_FluxAverage(Data *, Grid *, RBox *, RBox *);

// void OutflowBoundary    (const Ary4D &, const Ary3D (&)[6], const Ary3D &, RBox *, int);
// void PeriodicBoundary   (const Ary4D &,const Ary3D (&)[6], RBox *, int, Grid *);
// void ReflectiveBoundary (const Ary4D &,const Ary3D (&)[6], int *, RBox *, int, Grid *);
void PointValue (Data *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/* Inline functions                                                   */
////////////////////////////////////////////////////////////////////////

#define WENO3_METHOD  1

static const double thirteen_12 = 13.0/12.0;

static inline double WENO3_Average(double Fm, double Fc, double Fp, 
                                   double d0, double d1, double w)
{
  double b0, C0, a0, w0;
  double b1, C1, a1, w1, tau;
  double dx = 1.0/RuntimeGet()->npoint[IDIR];
  double dx2 = 1.e-4*dx*dx;
  double eps = 1.e-18;

  b0 = Fc - Fm; 
  b0 = b0*b0;

  b1 = Fp - Fc; 
  b1 = b1*b1; 

// Reconstruct L --> R

#if WENO3_METHOD == 0  // "Classical" WENO3
  a0 = eps + b0;
  a1 = eps + b1;
  C0 = 1.0/(a0*a0);
  C1 = 1.0/(a1*a1);
#elif WENO3_METHOD == 1  // "Improved" (YC) WENO3
  tau = Fp - 2.0*Fc + Fm; 
  tau = tau*tau;
  C0  = 1.0 + tau/(dx2 + b0);
  C1  = 1.0 + tau/(dx2 + b1);
#elif WENO3_METHOD == 4  //  "WENO3-Z" of Li et al with p = 3/2
  tau = fabs(b0 - b1);
  tau = tau*sqrt(tau);
  C0  = 1.0 + tau/(eps + b0);
  C1  = 1.0 + tau/(eps + b1);
#endif

  a0  = d0*C0;
  a1  = d1*C1;

  w0 = a0/(a0 + a1);
  w1 = a1/(a0 + a1);
  double Qp = 0.5*(w0*(-Fm + 3.0*Fc) + w1*(Fp + Fc));  

// Reconstruct L <-- R
  a0  = d0*C1;
  a1  = d1*C0;

  w0 = a0/(a0 + a1);
  w1 = a1/(a0 + a1);
  double Qm =  0.5*(w0*(-Fp + 3.0*Fc) + w1*( Fm + Fc));  

  return w*(Qp - 2.0*Fc + Qm);

/*
L->R:   d0*C0,  d1*C1,   a0 + a1 = d0*C0 + d1*C1
R->L:   d0*C1,  d1*C0,   a0 + a1 = d0*C1 + d1*C0


Qp + Qm =   [ d0*C0*(-Fm + 3*Fc) + d1*C1*(Fp+Fc) ] / (d0*C0 + d1*C1)  
          + [ d0*C1*(-Fp + 3*Fc) + d1*C0*(Fm+Fc) ] / (d0*C1 + d1*C0)    
*/

}


static inline double WENO5_Average(double Fmm, double Fm, double Fc, 
                                   double Fp, double Fpp, 
                                   double dL, double dR)
{
  double dC   = 1.0 - dL - dR;
  double d2QL = Fmm - 2.0*Fm + Fc;
  double d2QC = Fm  - 2.0*Fc + Fp;
  double d2QR = Fc  - 2.0*Fp + Fpp;
  double b0, b1, b2;
  double a0, a1, a2;
  double sum_a, t5;
  
  a1   = Fmm - 4.0*Fm + 3.0*Fc; 
  b0   = thirteen_12*d2QL*d2QL + 0.25*a1*a1;

  a1   = Fm - Fp; 
  b1   = thirteen_12*d2QC*d2QC + 0.25*a1*a1;

  a1   = 3.0*Fc - 4.0*Fp + Fpp; 
  b2   = thirteen_12*d2QR*d2QR + 0.25*a1*a1;

// "Classical" WENO-5 Version 
  // a0 = dL/(1.e-16 + b0);
  // a1 = dC/(1.e-16 + b1);
  // a2 = dR/(1.e-16 + b2);

// Improved version (WENO-Z)
  t5 = fabs(b0-b2);
  a0 = dL*(1.0 + t5/(b0+1.e-40));
  a1 = dC*(1.0 + t5/(b1+1.e-40));
  a2 = dR*(1.0 + t5/(b2+1.e-40));

  sum_a = a0 + a1 + a2;
  return (a0*d2QL + a1*d2QC + a2*d2QR)/sum_a;
}

static inline double WENO5_SplitAverage(double Fmm, double Fm, double Fc, 
                                        double Fp, double Fpp, 
                                        double dL, double dR)
{
  double dC   = 1.0 - dL - dR;
  double d2QL = Fmm - 2.0*Fm + Fc;
  double d2QC = Fm  - 2.0*Fc + Fp;
  double d2QR = Fc  - 2.0*Fp + Fpp;
  double b0, b1, b2;
  double a0, a1, a2;
  double sum_a, t5;
  
  a1   = Fmm - 4.0*Fm + 3.0*Fc; 
  b0   = thirteen_12*d2QL*d2QL + 0.25*a1*a1;

  a1   = Fm - Fp; 
  b1   = thirteen_12*d2QC*d2QC + 0.25*a1*a1;

  a1   = 3.0*Fc - 4.0*Fp + Fpp; 
  b2   = thirteen_12*d2QR*d2QR + 0.25*a1*a1;

  double theta = 3.0;
  double dL_p = 0.5*(dL+theta*fabs(dL)), dL_m = dL_p - dL;
  double dC_p = 0.5*(dC+theta*fabs(dC)), dC_m = dC_p - dC;
  double dR_p = 0.5*(dR+theta*fabs(dR)), dR_m = dR_p - dR;

// "Classical" WENO-5 Version 
  // a0 = dL_p/(1.e-16 + b0);
  // a1 = dC_p/(1.e-16 + b1);
  // a2 = dR_p/(1.e-16 + b2);

// Improved version (WENO-Z)
  t5 = fabs(b0-b2);
  a0 = dL_p*(1.0 + t5/(b0+1.e-40));
  a1 = dC_p*(1.0 + t5/(b1+1.e-40));
  a2 = dR_p*(1.0 + t5/(b2+1.e-40));

  sum_a = a0 + a1 + a2;
  double d2Q_p = (a0*d2QL + a1*d2QC + a2*d2QR)/sum_a;

// "Classical" WENO-5 Version 
  // a0 = dL_m/(1.e-16 + b0);
  // a1 = dC_m/(1.e-16 + b1);
  // a2 = dR_m/(1.e-16 + b2);

// Improved version (WENO-Z)
  t5 = fabs(b0-b2);
  a0 = dL_m*(1.0 + t5/(b0+1.e-40));
  a1 = dC_m*(1.0 + t5/(b1+1.e-40));
  a2 = dR_m*(1.0 + t5/(b2+1.e-40));

  sum_a = a0 + a1 + a2;
  double d2Q_m = (a0*d2QL + a1*d2QC + a2*d2QR)/sum_a;

  double sigp = dL_p + dC_p + dR_p;
  double sigm = dL_m + dC_m + dR_m;

  return sigp*d2Q_p - sigm*d2Q_m;
}

////////////////////////////////////////////////////////////////////////
/* De-average functions with 4 indices                                */
////////////////////////////////////////////////////////////////////////

// Function prototype
template <int dir>
inline double HO_DeAverage4D(Ary4D &Q, int, int, int, int) ;

// Explicit template instantiation for Ary4D
template double HO_DeAverage4D<IDIR>(Ary4D &Q, int i, int j, int k, int nv);
template double HO_DeAverage4D<JDIR>(Ary4D &Q, int i, int j, int k, int nv);
template double HO_DeAverage4D<KDIR>(Ary4D &Q, int i, int j, int k, int nv);

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template <int dir>
inline double HO_DeAverage4D(Ary4D &Q, int i, int j, int k, int nv) 
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int ip1 = i, jp1 = j, kp1 = k;
  int im1 = i, jm1 = j, km1 = k;
  int ip2 = i, jp2 = j, kp2 = k;
  int im2 = i, jm2 = j, km2 = k;
  #if HO_DEAVERAGE == 3
  const double dL = 1.0/6.0;
  const double dR = 1.0/6.0;
  #elif HO_DEAVERAGE == 4
  const double dL = 1.0/3.0;
  const double dR = 2.0/3.0; 
  #elif HO_DEAVERAGE == 5
  static const double dL = -9.0/80.0;
  static const double dR = -9.0/80.0;
  #endif
  double d2Q;

  #if HO_DEAVERAGE == 0
  return 0.0;
  #endif

  if constexpr (dir == IDIR) {
    ip1 = i+1; im1 = i-1;
    ip2 = i+2; im2 = i-2;
  } else if constexpr (dir == JDIR) {
    jp1 = j+1; jm1 = j-1;
    jp2 = j+2; jm2 = j-2;
  } else if constexpr (dir == KDIR) {
    kp1 = k+1; km1 = k-1;
    kp2 = k+2; km2 = k-2;
  }

#if HO_DEAVERAGE == 1
  d2Q = Q(ip1,jp1,kp1,nv) - 2.0*Q(i,j,k,nv) + Q(im1,jm1,km1,nv);
//  d2Q = Q(im1,jm1,km1,nv) - 2.0*Q(i,j,k,nv) + Q(ip1,jp1,kp1,nv);
#elif HO_DEAVERAGE == 3
  d2Q = WENO5_Average (Q(im2,jm2,km2,nv), Q(im1,jm1,km1,nv), Q(i,j,k,nv), 
                       Q(ip1,jp1,kp1,nv), Q(ip2,jp2,kp2,nv), dL, dR);
#elif HO_DEAVERAGE == 4  
  d2Q = WENO3_Average(Q(im1, jm1, km1,nv), Q(i, j, k,nv), Q(ip1, jp1, kp1,nv), 
                      dL, dR, 6.0);
#elif HO_DEAVERAGE == 5
  d2Q = WENO5_SplitAverage (Q(im2,jm2,km2,nv), Q(im1,jm1,km1,nv), Q(i,j,k,nv), 
                            Q(ip1,jp1,kp1,nv), Q(ip2,jp2,kp2,nv), dL, dR);
#endif

  return d2Q;
}

////////////////////////////////////////////////////////////////////////
/* De-average functions with 3 indices                                */
////////////////////////////////////////////////////////////////////////

// Function prototype
template <int dir>
inline double HO_DeAverage3D(Ary3D &Q, int, int, int k); 

// Explicit template instantiation for Ary3D
template double HO_DeAverage3D<IDIR>(Ary3D &Q, int i, int j, int k);
template double HO_DeAverage3D<JDIR>(Ary3D &Q, int i, int j, int k);
template double HO_DeAverage3D<KDIR>(Ary3D &Q, int i, int j, int k);

template <int dir>
inline double HO_DeAverage3D(Ary3D &Q, int i, int j, int k) 
{
  int ip1 = i, jp1 = j, kp1 = k;
  int im1 = i, jm1 = j, km1 = k;
  int ip2 = i, jp2 = j, kp2 = k;
  int im2 = i, jm2 = j, km2 = k;

  #if HO_DEAVERAGE == 3
  static const double dL = 1.0/6.0;
  static const double dR = 1.0/6.0;
  #elif HO_DEAVERAGE == 4
  const double dL = 1.0/3.0;
  const double dR = 2.0/3.0; 
  #elif HO_DEAVERAGE == 5
  static const double dL = -9.0/80.0;
  static const double dR = -9.0/80.0;
  #endif
  double d2Q;

  #if HO_DEAVERAGE == 0
  return 0.0;
  #endif

  if constexpr (dir == IDIR) {
    ip1 = i+1; im1 = i-1;
    ip2 = i+2; im2 = i-2;
  } else if constexpr (dir == JDIR) {
    jp1 = j+1; jm1 = j-1;
    jp2 = j+2; jm2 = j-2;
  } else if constexpr (dir == KDIR) {
    kp1 = k+1; km1 = k-1;
    kp2 = k+2; km2 = k-2;
  }

#if HO_DEAVERAGE == 1
  d2Q = Q(ip1,jp1,kp1) - 2.0*Q(i,j,k) + Q(im1,jm1,km1);
//  d2Q = Q(im1,jm1,km1) - 2.0*Q(i,j,k) + Q(ip1,jp1,kp1);
#elif HO_DEAVERAGE == 3
  d2Q = WENO5_Average (Q(im2,jm2,km2), Q(im1,jm1,km1), Q(i,j,k), 
                       Q(ip1,jp1,kp1), Q(ip2,jp2,kp2), dL, dR);
#elif HO_DEAVERAGE == 4  
  d2Q = WENO3_Average (Q(im1,jm1,km1), Q(i,j, k), Q(ip1, jp1, kp1), 
                       dL, dR, 6.0);
#elif HO_DEAVERAGE == 5
  d2Q = WENO5_SplitAverage (Q(im2,jm2,km2), Q(im1,jm1,km1), Q(i,j,k), 
                            Q(ip1,jp1,kp1), Q(ip2,jp2,kp2), dL, dR);
#endif

  return d2Q;
}

////////////////////////////////////////////////////////////////////////
/* Average functions with 4 indices                                   */
////////////////////////////////////////////////////////////////////////

// Function prototype
template <int dir>
inline double HO_Average4D(Ary4D &Q, int, int, int, int) ;

// Explicit template instantiation for Ary4D
template double HO_Average4D<IDIR>(Ary4D &Q, int i, int j, int k, int nv);
template double HO_Average4D<JDIR>(Ary4D &Q, int i, int j, int k, int nv);
template double HO_Average4D<KDIR>(Ary4D &Q, int i, int j, int k, int nv);

template <int dir>
inline double HO_Average4D(Ary4D &Q, int i, int j, int k, int nv) 
{
  int ip1 = i, jp1 = j, kp1 = k;
  int im1 = i, jm1 = j, km1 = k;
  int ip2 = i, jp2 = j, kp2 = k;
  int im2 = i, jm2 = j, km2 = k;

  #if HO_AVERAGE == 3
  double dL = 1.0/6.0;
  double dR = 1.0/6.0;
  #elif HO_AVERAGE == 4
  const double dL = 1.0/4.0;
  const double dR = 3.0/4.0; 
  #elif HO_AVERAGE == 5
  double dL = -17.0/240.0;
  double dR = -17.0/240.0;
  #endif
  double d2Q;

  #if HO_AVERAGE == 0
  return 0.0;
  #endif

  if constexpr (dir == IDIR) {
    ip1 = i+1; im1 = i-1;
    ip2 = i+2; im2 = i-2;
  } else if constexpr (dir == JDIR) {
    jp1 = j+1; jm1 = j-1;
    jp2 = j+2; jm2 = j-2;
  } else if constexpr (dir == KDIR) {
    kp1 = k+1; km1 = k-1;
    kp2 = k+2; km2 = k-2;
  }

#if HO_AVERAGE == 1
  d2Q = Q(ip1,jp1,kp1,nv) - 2.0*Q(i,j,k,nv) + Q(im1,jm1,km1,nv);
//  d2Q = Q(im1,jm1,km1,nv) - 2.0*Q(i,j,k,nv) + Q(ip1,jp1,kp1,nv);
#elif HO_AVERAGE == 3
  d2Q = WENO5_Average (Q(im2,jm2,km2,nv), Q(im1,jm1,km1,nv), Q(i,j,k,nv), 
                       Q(ip1,jp1,kp1,nv), Q(ip2,jp2,kp2,nv), dL, dR);
#elif HO_AVERAGE == 4  
  d2Q = WENO3_Average(Q(im1, jm1, km1,nv), Q(i, j, k,nv), Q(ip1, jp1, kp1,nv), 
                      dL, dR, 4.0);
#elif HO_AVERAGE == 5
  d2Q = WENO5_SplitAverage (Q(im2,jm2,km2,nv), Q(im1,jm1,km1,nv), Q(i,j,k,nv), 
                            Q(ip1,jp1,kp1,nv), Q(ip2,jp2,kp2,nv), dL, dR);
#endif

  return d2Q;
}

////////////////////////////////////////////////////////////////////////
/* Average functions with 3 indices                                   */
////////////////////////////////////////////////////////////////////////

// Function prototype
template <int dir>
inline double HO_Average3D(Ary3D &Q, int, int, int); 

// Explicit template instantiation for Ary3D
template double HO_Average3D<IDIR>(Ary3D &Q, int i, int j, int k);
template double HO_Average3D<JDIR>(Ary3D &Q, int i, int j, int k);
template double HO_Average3D<KDIR>(Ary3D &Q, int i, int j, int k);

template <int dir>
inline double HO_Average3D(Ary3D &Q, int i, int j, int k) 
{
  int ip1 = i, jp1 = j, kp1 = k;
  int im1 = i, jm1 = j, km1 = k;
  int ip2 = i, jp2 = j, kp2 = k;
  int im2 = i, jm2 = j, km2 = k;
  #if HO_AVERAGE == 3
  static const double dL = 1.0/6.0;
  static const double dR = 1.0/6.0;
  #elif HO_AVERAGE == 4
  const double dL = 1.0/4.0;
  const double dR = 3.0/4.0; 
  #elif HO_AVERAGE == 5
  static const double dL = -17.0/240.0;
  static const double dR = -17.0/240.0;
  #endif
  double d2Q;

  #if HO_AVERAGE == 0
  return 0.0;
  #endif

  if constexpr (dir == IDIR) {
    ip1 = i+1; im1 = i-1;
    ip2 = i+2; im2 = i-2;
  } else if constexpr (dir == JDIR) {
    jp1 = j+1; jm1 = j-1;
    jp2 = j+2; jm2 = j-2;
  } else if constexpr (dir == KDIR) {
    kp1 = k+1; km1 = k-1;
    kp2 = k+2; km2 = k-2;
  }

#if HO_AVERAGE == 1
  d2Q = Q(ip1,jp1,kp1) - 2.0*Q(i,j,k) + Q(im1,jm1,km1);
//  d2Q = Q(im1,jm1,km1) - 2.0*Q(i,j,k) + Q(ip1,jp1,kp1);
#elif HO_AVERAGE == 3
  d2Q = WENO5_Average (Q(im2,jm2,km2), Q(im1,jm1,km1), Q(i,j,k), 
                       Q(ip1,jp1,kp1), Q(ip2,jp2,kp2), dL, dR);
#elif HO_AVERAGE == 4  
  d2Q = WENO3_Average(Q(im1, jm1, km1), Q(i, j, k), Q(ip1, jp1, kp1), 
                      dL, dR, 4.0);
#elif HO_AVERAGE == 5
  d2Q = WENO5_SplitAverage (Q(im2,jm2,km2), Q(im1,jm1,km1), Q(i,j,k), 
                            Q(ip1,jp1,kp1), Q(ip2,jp2,kp2), dL, dR);
#endif

  return d2Q;
}

