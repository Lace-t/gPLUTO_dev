////////////////////////////////////////////////////////////////////////
/* gPLUTO function prototypes                                         */
////////////////////////////////////////////////////////////////////////

uint32_t AdvanceStep      (Data *, Grid *, timeStep *);
uint32_t AdvanceStepEuler (Data *, Grid *, timeStep *);

void    AllocateSweepArrays (Sweep *, Grid *);

void     Analysis (Data *, Grid *);
char    *Array1D (int, size_t);
char   **Array2D (int, int, size_t);
char  ***Array3D (int, int, int, size_t);
char ****Array4D (int, int, int, int, size_t);

template<int dir> void ArrayReconstruct(const Data *, const Ary3D &, int, 
                                        int, Grid *, RBox *, const int );

// #pragma acc routine vector
// #pragma omp declare target
// void ArrayReconstruct(double *, const Ary3D &, uint32_t *, int, int, int, int, double *, double *, int, Grid *);
// #pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
void   BackgroundField (State &, RunConfig &);
#pragma omp end declare target

unsigned int BitPosition (uint32_t);
void   BinaryPrint (uint32_t);

#pragma omp declare target
#pragma acc routine seq
double BodyForcePotential (const State&, const RunConfig &);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
void   BodyForceVector (const State&, const RunConfig &, double *);
#pragma omp end declare target

void  Boundary (Data *, int, Grid *);
void  BoundaryConservative (Data *, Grid *);
void  ClockBeg (string);
void  ClockEnd (string);
int   CheckData (Ary4D &, Ary4D &, Ary3D (&)[6], int, RBox *, const char *, Grid *);
void  CheckFluxSymmetry(const Data *, RBox *, Grid *);

#pragma omp declare target
#pragma acc routine seq
uint32_t ConsToPrim (State &, RunConfig &);
#pragma omp end declare target

void  ConsToPrim3D (Data *, int, Grid *, RBox *);
float *Convert_dbl2flt (double *, double, int, Grid *);
float ***Convert_dbl2flt_AL (double ***, double, int, Grid *);

template<int dir>
void DataToSweep (Data *, Grid *, RBox *);

string ErrorShow (uint32_t);

void   FillOutputBuffer_AL (int, Output *, RBox *);
void   FillOutputBuffer (int, Output *, RBox *);
void   Fallback (Data *, Grid *, timeStep *);
void   Fallback_Trigger (Ary3Duint32_t&, Ary3Duint32_t&, RBox *);
int    FileClose (FILE *, int) ;
int    FileDelete (char *);
FILE * FileOpen (char *, int, const char *);
void   FileReadData_AL (void *, size_t, int, FILE *, int, int, Grid *);
void   FileReadData (void *, size_t, int, int, Grid *, char *, FileOffset);
void   FileWriteData(void *, size_t, FileDataType, int, Grid *, char *, FileOffset);
void   FileWriteData_AL (void *, size_t, int, FILE *, int, Grid *);

void   FileWriteHeader (char *, char *, int);
void   FileWriteArray (void *, long int, long int, size_t, char *);
void   FlagBoundary(Ary3Duint32_t &, const int, Grid *);

void   FlagShock (Data *, Grid *, RBox *);
string FlagShow (uint32_t);
void   FlipSign (int, int, int *);
void   FreeArray1D (void *);
void   FreeArray2D (void **);
void   FreeArray3D (void ***);
void   FreeArray4D (void ****);

#pragma omp declare target
#pragma acc routine seq
void  GaussianElimination (double (*)[3], double *, double *, int);
#pragma omp end declare target

void  GetOutputFrequency (Output *, const char *);
void  GetCGSUnits (double *, RunConfig &);
void  GetNeighbourRanks (const Grid *, int **, int ***, int ***);
int   GetNghost (void);
void  GnuplotSetting (Runtime *, Grid *);

char  *IndentString ();
void   IndexCycle (int, int *, int *, int *);
#if ACC_INIT
#pragma omp declare target
#pragma acc routine seq
#endif
void   Init (double *, double, double, double, RunConfig &);
#if ACC_INIT
#pragma omp end declare target
#endif
void   InitDomain (Data *, Grid*);
void   Initialize (Data *, timeStep *, cmdLine *, Runtime *, RunConfig *, Grid *);
void   InputDataClose (int);
void   InputDataGridSize (int, int *);
double InputDataInterpolate (int, double, double, double);
int    InputDataOpen (const char *, const char *, const char *, long int, int);
void   InputDataReadSlice (int, int);int  IsLittleEndian (void) ;
template<int dir> void InternalBoundaryReflect (Data *, Grid *, RBox *);

void   mynvtxstart_ (const char *, int);
void   mynvtxstop_ ();
void   measure_acc_mem_usage ();

// #pragma omp declare target
// #pragma acc routine seq
// double MP5_States (double *, int);
// #pragma omp end declare target

// #pragma omp declare target
// #pragma acc routine seq
// double WENOZ_States (double *, int);
// #pragma omp end declare target

// #pragma omp declare target
// #pragma acc routine seq
// double WENO3_States (double *, int, double);
// #pragma omp end declare target

// #pragma omp declare target
// #pragma acc routine seq
// double LimO3_States (double *, int, double);
// #pragma omp end declare target

void  LinkUserVar (Ary3D *, const char *, int);
void  LogFileClose (void);
void  LogFileFlush (void);
void  LogFileOpen (char *, const char *);

#pragma omp declare target
#pragma acc routine seq
int   LUDecompose (double *, int, int *, double *);
#pragma acc routine seq
void  LUBackSubst (double *, int, int *, double b[]);
#pragma acc routine seq
double MeanMolecularWeight(const double *);
#pragma omp end declare target

double NextTimeStep (timeStep *, RunConfig *, Runtime *, Grid *);

void  OutputBin (double ****, Grid *);
void  OutputLogPost (double *, double *, double);

void  OutputGrid (Grid *, Runtime *);

void  print (const char *fmt, ...);
void  printLog (const char *fmt, ...);
void  ParabolicFlux_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir);
void  ParabolicFlux (Data *data, timeStep *Dts, Grid *grid, RBox *box);
void  ParallelExchange (const Ary4D &, const Ary3D (&)[6], Grid *);

template <int flag_op> 
void  ParallelExchangeFlag(Ary3Duint32_t &, const Grid *);

void  ParallelExchangeDouble (Ary3D &, Grid *);
void  PrintColumnLegend (const char *legend[], int, FILE *);
void  print_mem_usage();
int   ParamFileRead    (char *);
char *ParamFileGet     (const char *, int );
int   ParamExist       (const char *);
int   ParamFileHasBoth (const char *, const char *);
void  ParseCmdLineArgs (int argc, char *argv[], cmdLine *);
void  Predict (Data *, Grid *);
//void  PrimToCons3D (Ary4D &, Ary4D &, Grid *, RunConfig *, RBox *);
void  PrimToCons3D (const Data *, int, Grid *, RBox *);

void  QuitPLUTO(string, int);



#if GEOMETRY == CYLINDRICAL || (GEOMETRY == SPHERICAL && DIMENSIONS == 3)
void SingularAxisBoundary   (const Ary4D &, const Ary3D (&)[6], RBox *, int, Grid *);
#endif
void OutflowBoundary    (const Ary4D &, const Ary3D (&)[6], const Ary3D &, 
                         RBox *, int);
void PeriodicBoundary   (const Ary4D &, const Ary3D (&)[6], RBox *, int, Grid *);
void ReflectiveBoundary (const Ary4D &, const Ary3D (&)[6], int *, RBox *, int);
void FlipSign (int, int, int *);




#pragma omp declare target
#pragma acc routine seq
void  RBoxDefine (int, int, int, int, int, int, RBox *);
#pragma omp end declare target

void  RBoxDefine (int *, int *, RBox *);
void  RBoxGetIndexes (RBox *, int *, int *);
long int RBoxGetSize (RBox *);
void  RBoxOffset (RBox *, int *);
template<int dir>
void  RBoxPermuteIndices (RBox *, int &, int &, int &, int &, int &, int &);
void  RBoxShow(RBox *);

void  ReadOutputBuffer_AL (int, Output *, RBox *);
void  ReadOutputBuffer (int, Output *, RBox *);
void  ReadHDF5 (Output *output, Grid *grid);
template<int dir> void Reconstruct     (Data *, Grid *, RBox *);
template<int dir> void ReconstructFlat (Data *, Grid *, RBox *);

void  RestartFromFile (Data *, Runtime *, timeStep *, int, int, Grid *);
void  RestartDump (Runtime *, double *, timeStep *);
void  RestartGet (Data *, Runtime *, timeStep *, int, int, int);

void  RightHandSide (Data *, timeStep *, const Grid *const, double, RBox *, const int);
void  RightHandSideSource (Data *, Grid *, double, RBox *box, const int);

#if RING_AVERAGE > 1
void RingAverageCons(Data *, Ary3D (&)[6], Grid *);
template<int dir> void RingAverageReconstruct (Data *, Ary3D (&)[6], Grid *, RBox *box);
template<int dir> void CT_RingAverageReconstruct (Data *, Ary3D (&)[6], Grid *, RBox *box);
void RingAverageReconstructNew (double *v, int beg, int end, Grid *grid);
void RingAverageSize(Grid *);
#endif

Runtime *RuntimeGet (void);
int   RuntimeSetup (Runtime *, cmdLine *, char *);
void  RuntimeSet (Runtime *runtime);

RiemannSolver *SetSolver (const char *solver);
#if RADIATION
RiemannSolver *RadSetSolver (const char *solver);
#endif
void  SetDefaultVarNames (Output *output);
void  SetJetDomain (Data const*, int, int, Grid *);  
#if PHYSICS == GRMHD
void  SetGeometry (Grid *, RunConfig &);
void  SetGrid (Runtime *, int *, Grid *, int, RunConfig &);
#else
void  SetGeometry (Grid *);
void  SetGrid (Runtime *, int *, Grid *, int);
#endif
void  SetOutput (Data *, Runtime *, Grid *);
int   SetOutputVar (const char *, int, int);
void  ShowConfig (char *, RunConfig &);
void  ShowData (Ary4D &, Ary4D &, Ary3D (&)[6], RBox *, const char *);
void  ShowMatrix (double **, int, double);
void  ShowState (double *vc, int);
void  ShowVector (double *v, int);
void  ShowUnits (RunConfig &);
int   StringArrayIndex (char **, int , char *);
void  Sources (Data *, Ary4D &, double, Grid *, RBox *);
void  SplitSource (Data *, double , timeStep *, Grid *);
void  Startup (Data *, Grid *);
void  SwapEndian (void *, const int);

#pragma omp declare target
#pragma acc routine seq
void  SwapRows (double (*)[3], double *, int, int, int);
#pragma omp end declare target


void  UnsetJetDomain (int, Grid *);  
void  UpdateStage (Data *, Ary4D &, Ary3D (&)[6], Grid *, timeStep *, RBox *, double);
void  UpdateStage (Data *, Ary4D &, Ary3D (&)[6], Grid *, timeStep *, RBox *, double,
                   bool pcp_flux_possible);
void  UserDefBoundary (const Data *, RBox *, int, Grid *);
void  UserDefOutput   (const Data *, Grid *);
void  UserDefOutput   (const Data *, Grid *);
void  VectorTransform(double *, double, double, double, int);

#if PHYSICS == GRMHD
void  VectorCurl(Ary3D *, Ary3D *, Grid *, RBox *, RunConfig &);
#else
void  VectorCurl(Ary3D *, Ary3D *, Grid *, RBox *);
#endif

void  Where (const string, Data *, uint32_t, RBox *, Grid *);
void  WriteData (Data *, timeStep *, Output *, Grid *);
void  WriteHDF5 (Output *, timeStep *, Grid *grid);
void  WriteTabArray (Output *, char *, Grid *);
void  WriteVTK_Header (char *, timeStep *, Grid *);
void  WriteVTK_Vector (FILE *, double ****, double, char *, Grid *);
void  WriteVTK_Scalar (char *, double *, double, char *, Grid *);
void  WriteVTK_Scalar_AL (char *, double ***, double, char *, Grid *);

#pragma omp declare target
#pragma acc routine seq
size_t idx4 (int, int, int, int, Grid *);
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
double* elem4 (double*, int, int, int, int, Grid *);
#pragma omp end declare target

#define ARRAY_1D(nx,type)          (type    *)Array1D(nx,sizeof(type))
#define ARRAY_2D(nx,ny,type)       (type   **)Array2D(nx,ny,sizeof(type))
#define ARRAY_3D(nx,ny,nz,type)    (type  ***)Array3D(nx,ny,nz,sizeof(type))
#define ARRAY_4D(nx,ny,nz,nv,type) (type ****)Array4D(nx,ny,nz,nv,sizeof(type))
#define ARRAY_BOX(i0,i1, j0,j1, k0,k1,type)  \
        (type ***)ArrayBox(i0,i1,j0,j1,k0,k1,sizeof(type))



//void  PowerLawCooling (Ary4D &, double, timeStep &, Grid &);

