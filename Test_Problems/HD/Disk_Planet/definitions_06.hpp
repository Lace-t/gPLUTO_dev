#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

// Physics dependent declarations

#define  DUST_FLUID                     NO
#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  Mstar                          0
#define  Mdisk                          1
#define  Mplanet                        2
#define  Viscosity                      3

// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        VANLEER_LIM
#define  ISO_SOUNDSPEED2                ((CONST_PI*0.1)*(CONST_PI*0.1))

// [End] user-defined constants (do not change this line)
