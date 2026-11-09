#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            10

// Physics dependent declarations

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  RHO_LEFT                       0
#define  VX1_LEFT                       1
#define  VX2_LEFT                       2
#define  VX3_LEFT                       3
#define  PRS_LEFT                       4
#define  RHO_RIGHT                      5
#define  VX1_RIGHT                      6
#define  VX2_RIGHT                      7
#define  VX3_RIGHT                      8
#define  PRS_RIGHT                      9

// [Beg] user-defined constants (do not change this line)

#define  GNUPLOT_HEADER                 YES

// [End] user-defined constants (do not change this line)
