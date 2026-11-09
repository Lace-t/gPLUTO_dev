#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        2
#define  PARTICLES                      PARTICLES_LP
#define  USER_DEF_PARAMETERS            11

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   DIV_CLEANING

/* -- Particles Options -- */

#define  PARTICLES_LP_SPECTRA           YES
#define  PARTICLES_LP_SHOCK_SPECTRA     YES
#define  PARTICLES_LP_MREC_SPECTRA      NO

/* -- user-defined parameters (labels) -- */

#define  GAMMA_UP                       0
#define  NU                             1
#define  BETA                           2
#define  RADX                           3
#define  RADY                           4
#define  GAMMA_DOWN                     5
#define  RADZ                           6
#define  THETA                          7
#define  MACH                           8
#define  SIGMA                          9
#define  TEMP                           10

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  DEBUG_MIN_MAX                  NO
#define  SHOCK_FLATTENING               MULTID
#define  LIMITER                        VANLEER_LIM

/* [End] user-defined constants (do not change this line) */
