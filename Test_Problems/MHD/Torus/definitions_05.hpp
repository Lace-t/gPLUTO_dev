#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               YES
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  RMIN                           0
#define  RMAX                           1
#define  RHO_CUT                        2
#define  BETA                           3
#define  ETA                            4
#define  SCALE_HEIGHT                   5

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CHECK_DIVB_CONDITION           NO
#define  CT_EMF_AVERAGE                 UCT_HLL
#define  LIMITER                        VANLEER_LIM
#define  USE_DIPOLE                     YES
#define  ACC_INIT                       NO

// [End] user-defined constants (do not change this line)
