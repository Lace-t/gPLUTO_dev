#define  PHYSICS                        RHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

// Physics dependent declarations

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO

// User-defined parameters (labels)

#define  BETA                           0
#define  RHO_IN                         1
#define  RHO_OUT                        2
#define  PRESS_IN                       3

// [Beg] user-defined constants (do not change this line)

#define  ARTIFICIAL_VISC                0.1
#define  CHAR_LIMITING                  YES
#define  EPS_PSHOCK_FLATTEN             2.5
#define  INCLUDE_JDIR                   NO
#define  LIMITER                        MC_LIM
#define  SHOCK_FLATTENING               MULTID

// [End] user-defined constants (do not change this line)
