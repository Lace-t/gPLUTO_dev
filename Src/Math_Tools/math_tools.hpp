#ifndef MATH_TOOLS
#define MATH_TOOLS
////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Math tools header file.

  This module provides a number of standard numerical routines 
  to achieve simple basic tasks such as

  - Numerical quadrature: GaussQuadrature();
  - LU decomposition functions: LUDecompose(), LUBackSubst();
  - Scalar and multi-dim. root finder methods: Brent(), Ridder(), Brodyen();
  - Vector transformation, VectorCartesianComponents();
  - Ordinary differential equation solvers: ODESolve() with RK4Solve()
    or CK45Solve();
  - Table handling: LocateIndex(), Table2DInterpolate(), GenerateTable2D(),
    InverseLookupTable2D().         
  

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    March 09, 2026
*/
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/*     \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
       Function prototyping                                           */                                                               
////////////////////////////////////////////////////////////////////////

// --------------------------------------------------------
// Functions in math_lu_decomp.c
// --------------------------------------------------------

void GaussianElimination (double A[][3], double *, double *, int);
int  LUDecompose (double **, int, int *, double *);
void LUBackSubst (double **, int, int *, double b[]);
void SwapRows (double A[][3], double *, int, int, int);
void TridiagonalSolve(double *, double *, double *, double *,  double *, int);

// --------------------------------------------------------
// Functions in math_func.c
// --------------------------------------------------------

double BesselJ0(double);
double BesselJ1(double);
double BesselIO(double);
double BesselI1(double);
double BesselKO(double);
double BesselK1(double);
double BesselKn(int, double);
double Dawson(double);

// --------------------------------------------------------
// Functions in math_misc.c
// --------------------------------------------------------

void   QuickSort(int *, int, int);
void   SortArray (double *, int);
void   VectorCartesianComponents(double *, double, double, double);

// /* -------------------------------------------------------
//     Functions contained in math_quadrature.c
//    -------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////
inline void GaussianWeights (double *xg, double *wg, int Ng)
/*!
  Compute Gaussian abscissa and weights.
  
  \param [out]   *xg   an array of abscissas
  \param [out]   *wg   an array of weights
  \param [in]    Ng    the number of desired Gaussian points
*/
////////////////////////////////////////////////////////////////////////
{
  double one_third   =  1.0/3.0; 
  double three_fifth =  3.0/5.0, six_fifth     = 6.0/5.0;
  double two_seventh =  2.0/7.0, three_seventh = 3.0/7.0;
  double ten_seventh = 10.0/7.0;

  // ------------------------------------------------------
  // 0. Setting gaussian points and weights
  // ------------------------------------------------------

  if (Ng == 1) {
    xg[0] = 0.0;
    wg[0] = 2.0;
  } else if (Ng == 2) {     
    xg[0] = -sqrt(one_third);   xg[1] = sqrt(one_third);         
    wg[0] =  1.0;               wg[1] = 1.0;
  } else if (Ng == 3) { 
    xg[0] = -sqrt(three_fifth); xg[1] = 0.0;      xg[2] = sqrt(three_fifth);          
    wg[0] =  5.0/9.0;           wg[1] = 8.0/9.0;  wg[2] = 5.0/9.0;
  } else if (Ng == 4) {
    xg[0] = -sqrt(three_seventh-sqrt(six_fifth)*two_seventh);
    xg[1] =  sqrt(three_seventh-sqrt(six_fifth)*two_seventh);
    xg[2] = -sqrt(three_seventh+sqrt(six_fifth)*two_seventh);
    xg[3] =  sqrt(three_seventh+sqrt(six_fifth)*two_seventh);

    wg[0] = (18.0+sqrt(30.0))/36.0;
    wg[1] = (18.0+sqrt(30.0))/36.0; 
    wg[2] = (18.0-sqrt(30.0))/36.0;
    wg[3] = (18.0-sqrt(30.0))/36.0;
  } else if (Ng == 5) {
    xg[0] =  0.0;
    xg[1] =  sqrt(5.0-2.0*sqrt(ten_seventh))*one_third;
    xg[2] = -sqrt(5.0-2.0*sqrt(ten_seventh))*one_third;
    xg[3] =  sqrt(5.0+2.0*sqrt(ten_seventh))*one_third;
    xg[4] = -sqrt(5.0+2.0*sqrt(ten_seventh))*one_third;

    wg[0] =  128.0/225.0;
    wg[1] = (322.0+13.0*sqrt(70.0))/900.0;
    wg[2] = (322.0+13.0*sqrt(70.0))/900.0;
    wg[3] = (322.0-13.0*sqrt(70.0))/900.0;
    wg[4] = (322.0-13.0*sqrt(70.0))/900.0;
  } else {
//    QuitPLUTO("! GaussQuadrature(): Ng must be an integer 0 < Ng < 6",1);
  }
};

template <typename T>
double GaussQuadrature(double (*func)(double, T *), T *par, double xb, double xe,
                       int nstep, int order)
/*!
 * Perform numerical quadrature of the function f(x) between 
 * the lower bound xb and upper bound xe by subdividing the interval 
 * into 'nstep' steps.
 * A 3 or 5-point Gaussian quadrature rule is used depending on the 
 * input variable order (=3 or =5)
 * 
  \param [in] *func   a pointer to a generic function func(x, par) 
 *                    (returning double) to be integrated.
 *                     Here par is a generic void pointer that may (or not)
 *                     be used to pass parameters.
  \param [in] *par    a pointer to a generic void data type (to be passed 
 *                     to func)  
  \param [in] xb      the lower interval bound
  \param [in] xe      the upper interval bound
  \param [in] nstep   the number of sub-intervals into which the 
 *                     original interval [xb,xe] has to be divided
  \param [in] order   the number of Gaussian points (from 1 to 5)    
*/
////////////////////////////////////////////////////////////////////////
{
  int    i, n;
  double w[8], z[8], x;
  double I, Isub;
  double xb0, xe0, dx;
  
  GaussianWeights (z, w, order);
  
  xb0 = xb; xe0 = xe; // save original interval endpoints 
  dx = (xe - xb)/(double)nstep; // sub-interval length 
  I  = 0.0;
  for (i = 0; i < nstep; i++){
    xb = xb0 + i*dx;
    xe = xb + dx;
    
    Isub = 0.0;  // intgrate sub-interval 
    for (n = 0; n < order; n++){
      x     = 0.5*(xe - xb)*z[n] + (xe + xb)*0.5;
      Isub += w[n]*func(x, par);
    }    
    Isub *= 0.5*(xe - xb);
    I    += Isub;
  }
  
  return I;
};

// double GaussQuadrature(double (*func)(double, auto*), auto *, double, double, 
//                        int, int);

  // ------------------------------------------------
  // Functions contained in math_random.c
  // ------------------------------------------------

#define PRNG_DEFAULT   0  // Standard drand48() 
#define PRNG_ECUYER    1  // prng of L'Ecuyer   
#define PRNG_MT        2  // Mersenne Twister pseudorandom number 

#ifndef PRNG
  #define PRNG  PRNG_DEFAULT
#endif

#pragma acc routine seq
static inline uint64_t hashMix_new(uint64_t z) {
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    z ^= (z >> 31);
    return z;
}

#pragma acc routine seq
static inline double GridRandomNumber_new(int64_t ix, int64_t iy, int64_t iz, int g_rank) {
    uint64_t seed = (uint64_t)ix * 73856093ULL 
                  ^ (uint64_t)iy * 19349663ULL
                  ^ (uint64_t)iz * 83492791ULL
                  ^ (uint64_t)g_rank * 2654435761ULL;

    seed = hashMix_new(seed);
    return (double)(seed & 0x1FFFFFFFFFFFFFULL) / (double)0x1FFFFFFFFFFFFFULL; 
}


#pragma acc routine seq
double   GridRandomNumber(int64_t, int64_t, int64_t);
uint64_t hashMix(uint64_t);
double   GaussianRandomNumber(double, double);
void     RandomSeed (long int, long int);
// #pragma acc routine seq
double   RandomNumber (double, double);

#pragma omp declare target
#pragma acc routine seq
double InitGaussAverage(intList *nv_list, 
                        double xb, double xe, 
                        double yb, double ye, 
                        double zb, double ze, double *uav,
                        int Ng0, int Ng1, int Ng2, RunConfig &run_config);
#pragma omp end declare target

#pragma acc routine seq
double CellVolume(double, double, 
                  double, double, 
                  double, double,
                  int, int, int, RunConfig &);
#pragma acc routine seq
double CellBarycenter(double, double, 
                      double, double, 
                      double, double,
                      int, int, int, int, RunConfig &);
#pragma acc routine seq
double My_fmod(double, double);

double   PowerLawRandomNumber(double, double, double);
unsigned int SeedGenerator(long);
#if PRNG == PRNG_ECUYER
double   NR_ran2(long int *);
#elif PRNG == PRNG_MT
void               init_genrand64(unsigned long long);
unsigned long long genrand64_int64(void);
double             genrand64_real1(void);
#endif


// \endcond 

#endif
