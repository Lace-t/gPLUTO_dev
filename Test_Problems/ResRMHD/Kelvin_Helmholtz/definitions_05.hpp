#define  PHYSICS                        ResRMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENOZ
#define  TIME_STEPPING                  RK4
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  DIVE_CONTROL                   NO

// User-defined parameters (labels)

#define  SIGMA_POL                      0
#define  SIGMA_TOR                      1
#define  ETA                            2

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CT_EMF_AVERAGE                 CT_MAXWELL
#define  LIMITER                        MC_LIM
#define  HO_LAP_LIMITER                 HO_RDER_LIM
#define  HO_LAP_LIMITER_EPS             0.3
#define  HO_ORDER_REDUCTION             LINEAR
#define  PARALLEL_COMM_ASYNC            NO

// [End] user-defined constants (do not change this line)
