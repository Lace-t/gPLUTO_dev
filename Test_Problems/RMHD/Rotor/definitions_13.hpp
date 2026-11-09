#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   EIGHT_WAVES

// User-defined parameters (labels)


// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        OSPRE_LIM
#define  EIGHT_WAVES_REDUCED            YES
#define  WARNING_MESSAGES               2

// [End] user-defined constants (do not change this line)
