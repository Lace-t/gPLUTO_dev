#if MAPPED_GRID == YES
// Mapped Grid macro

#define MG_EXPONENTIAL  1
#define MG_SINUSOIDAL   2
#define MG_STRETCHED    3
#define MG_LOG          4 
#define MG_SUS          5
#define MG_UP           6
#define MG_US           7
#define MG_UNIFORM      9  
#define MG_USERDEF      99

// Function prototype
double MappedGrid (Grid *, int, double);

double MappedGrid_Exp(Grid *, int, double);
double MappedGrid_SUS(Grid *, int, double, double, double, double, double);
double MappedGrid_US       (Grid *, int, double, double, double);
double MappedGrid_UP       (Grid *, int, double, double, double);
double MappedGrid_Sin      (Grid *, int, double);
double MappedGrid_Sinh     (Grid *, int, double, double);
double MappedGrid_Stretched(Grid *, int, double, double);
double MappedGrid_Uniform  (Grid *, int, double);

double MappedGrid_StretchRatio(double);
#endif

////////////////////////////////////////////////////////////////////////
/*! 
  Define metric coefficients and their inverses. 
  Default values are already set to 1. 
  Here we only update coordinate-dependent scale factors.
*/
////////////////////////////////////////////////////////////////////////
#if GEOMETRY != CARTESIAN
inline void MetricSet (Metric &metric, const double r, const double s)
{
  double r_inv = (fabs(r) > 1.e-18 ? 1.0/r:0.0);
  #if GEOMETRY == CYLINDRICAL
  metric.h2     = fabs(r);
  metric.h2_inv = fabs(r_inv);
  #elif GEOMETRY == SPHERICAL
  double s_inv = (fabs(s) > 1.e-18 ? 1.0/s:0.0);
  metric.h2     = fabs(r);
  metric.h3     = fabs(r*s);
  metric.h2_inv = fabs(r_inv);
  metric.h3_inv = fabs(r_inv*s_inv);
  #endif
}
#endif
