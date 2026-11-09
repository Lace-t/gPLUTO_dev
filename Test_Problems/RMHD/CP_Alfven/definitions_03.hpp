#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 MP5
#define  TIME_STEPPING                  RK4
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

// User-defined parameters (labels)


// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CT_EMF_AVERAGE                 UCT_GFORCE
#define  RMHD_REDUCED_ENERGY            YES
#define  SHOW_TIME_STEPS                YES
#define  ACC_INIT                       NO

// [End] user-defined constants (do not change this line)
