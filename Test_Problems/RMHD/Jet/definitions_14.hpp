#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

// Physics dependent declarations

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   DIV_CLEANING

// User-defined parameters (labels)

#define  MACH                           0
#define  LORENTZ                        1
#define  RHOJ                           2
#define  RHOA                           3
#define  SIGMA_POL                      4
#define  SIGMA_TOR                      5

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  COVARIANT_MX2                  NO
#define  COVARIANT_BX2                  YES
#define  GLM_ALPHA                      0.2
#define  LIMITER                        VANLEER_LIM
#define  WARNING_MESSAGES               3

// [End] user-defined constants (do not change this line)
