#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            1

// Physics dependent declarations

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  SCRH                           0

// [Beg] user-defined constants (do not change this line)

#define  ARTIFICIAL_VISC                0.1
#define  INTERNAL_BOUNDARY              YES
#define  INTERNAL_BOUNDARY_REFLECT      YES
#define  SHOCK_FLATTENING               MULTID

// [End] user-defined constants (do not change this line)
