#define  PHYSICS                        ResRMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  DIVE_CONTROL                   CONSTRAINED_TRANSPORT

// User-defined parameters (labels)

#define  RHO_IN                         0
#define  PRS_IN                         1
#define  RHO_OUT                        2
#define  PRS_OUT                        3
#define  BMAG                           4
#define  RADIUS                         5
#define  ETA                            6

// [Beg] user-defined constants (do not change this line)

#define  CT_EMF_AVERAGE                 CT_MAXWELL
#define  CT_ENERGY_CORRECTION           YES
#define  LIMITER                        MINMOD_LIM

// [End] user-defined constants (do not change this line)
