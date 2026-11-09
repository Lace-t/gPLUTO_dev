#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        2
#define  PARTICLES                      PARTICLES_LP
#define  USER_DEF_PARAMETERS            9

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   EIGHT_WAVES

// Particles Options

#define  PARTICLES_LP_SPECTRA           YES
#define  PARTICLES_LP_SHOCK_SPECTRA     YES
#define  PARTICLES_LP_MREC_SPECTRA      NO

// User-defined parameters (labels)

#define  GAMMA_UP                       0
#define  NU                             1
#define  BETA                           2
#define  RADX                           3
#define  RADY                           4
#define  RADZ                           5
#define  SIGMA                          6
#define  TEMP                           7
#define  INJECTION_DT                   8

// [Beg] user-defined constants (do not change this line)

#define  INTERNAL_BOUNDARY              YES
#define  LIMITER                        UMIST_LIM
#define  GFORCE_OMEGA                   0.5
#define  SHOCK_FLATTENING               MULTID
#define  WARNING_MESSAGES               3

// [End] user-defined constants (do not change this line)
