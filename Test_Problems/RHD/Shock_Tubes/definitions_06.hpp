#define  PHYSICS                        RHD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            9

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 ALWAYS
#define  RADIATION                      NO

// User-defined parameters (labels)

#define  DN_L                           0
#define  VX_L                           1
#define  VY_L                           2
#define  PR_L                           3
#define  DN_R                           4
#define  VX_R                           5
#define  VY_R                           6
#define  PR_R                           7
#define  GAMMA_EOS                      8

// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        MC_LIM
#define  SHOCK_FLATTENING               NO

// [End] user-defined constants (do not change this line)
