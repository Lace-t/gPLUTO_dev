#define  PHYSICS                        ResRMHD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENOZ
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 ALWAYS
#define  DIVB_CONTROL                   NO
#define  DIVE_CONTROL                   NO

// User-defined parameters (labels)


// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CT_EMF_AVERAGE                 CT_MAXWELL
#define  LIMITER                        MC_LIM
#define  ACC_INIT                       NO

// [End] user-defined constants (do not change this line)
