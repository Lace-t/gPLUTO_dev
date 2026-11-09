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
#define  USER_DEF_PARAMETERS            4

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  DIVE_CONTROL                   NO

// User-defined parameters (labels)

#define  ETA                            0
#define  P0                             1
#define  EPS                            2
#define  CRG0                           3

// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        VANLEER_LIM
#define  CT_EMF_AVERAGE                 CT_MAXWELL

// [End] user-defined constants (do not change this line)
