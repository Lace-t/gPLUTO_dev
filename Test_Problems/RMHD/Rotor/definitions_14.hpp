#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

// User-defined parameters (labels)


// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        UMIST_LIM
#define  CT_EMF_AVERAGE                 UCT_HLL
#define  CT_ENERGY_CORRECTION           NO
#define  WARNING_MESSAGES               2

// [End] user-defined constants (do not change this line)
