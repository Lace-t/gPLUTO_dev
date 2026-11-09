#define  PHYSICS                        MHD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            17

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  RHO_LEFT                       0
#define  VX1_LEFT                       1
#define  VX2_LEFT                       2
#define  VX3_LEFT                       3
#define  BX2_LEFT                       4
#define  BX3_LEFT                       5
#define  PRS_LEFT                       6
#define  RHO_RIGHT                      7
#define  VX1_RIGHT                      8
#define  VX2_RIGHT                      9
#define  VX3_RIGHT                      10
#define  BX2_RIGHT                      11
#define  BX3_RIGHT                      12
#define  PRS_RIGHT                      13
#define  BX1_CONST                      14
#define  DIVIDE_BY_4PI                  15
#define  GAMMA_EOS                      16

// [Beg] user-defined constants (do not change this line)

#define  LIMITER                        MC_LIM

// [End] user-defined constants (do not change this line)
