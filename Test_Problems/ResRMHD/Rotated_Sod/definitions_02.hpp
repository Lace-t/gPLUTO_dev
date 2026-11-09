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
#define  USER_DEF_PARAMETERS            21

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  DIVE_CONTROL                   DIV_CLEANING

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
#define  GAMMA_EOS                      15
#define  INT_JX                         16
#define  INT_JY                         17
#define  INT_KX                         18
#define  INT_KZ                         19
#define  ETA                            20

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CT_EMF_AVERAGE                 CT_MAXWELL
#define  DIVIDE_BY_4PI                  NO
#define  LIMITER                        VANLEER_LIM

// [End] user-defined constants (do not change this line)
