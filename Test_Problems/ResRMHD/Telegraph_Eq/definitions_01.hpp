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
#define  USER_DEF_PARAMETERS            3

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  DIVE_CONTROL                   NO

// User-defined parameters (labels)

#define  ETA                            0
#define  THETA_YZ                       1
#define  APPLY_ROTATION                 2

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CT_EMF_AVERAGE                 CT_MAXWELL
#define  GLM_COMPUTE_DIVE               NO
#define  LIMITER                        VANLEER_LIM
#define  ACC_INIT                       NO

// [End] user-defined constants (do not change this line)
