#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

// Physics dependent declarations

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  ETA                            0
#define  MACH                           1
#define  VPHI                           2
#define  BMAG                           3

// [Beg] user-defined constants (do not change this line)

#define  INCLUDE_JDIR                   NO
#define  LIMITER                        MC_LIM

// [End] user-defined constants (do not change this line)
