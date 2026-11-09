void   RotateBoundary (const Ary4D &,const Ary3D (&)[6], RBox *, int, Grid *);
double RotateGet_ta(void);
double RotateGet_tb(void);
#if ACC_INIT
#pragma acc routine seq
#endif
void   RotateSet(int, int, int, int, RunConfig &);
#if ACC_INIT
#pragma acc routine seq
#endif
void   RotateVector(double *, int);
