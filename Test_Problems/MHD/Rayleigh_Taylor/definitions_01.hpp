#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            2

// Physics dependent declarations

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  ETA                            0
#define  CHI                            1

// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        MC_LIM
#define  USE_RANDOM_PERTURBATION        NO
#define  GRAV                           (-0.1)

// [End] user-defined constants (do not change this line)
