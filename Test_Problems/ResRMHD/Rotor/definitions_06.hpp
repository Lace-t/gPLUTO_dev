#define  PHYSICS                        ResRMHD
#define  DIMENSIONS                     3
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

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  DIVE_CONTROL                   NO

// User-defined parameters (labels)

#define  ETA                            0

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CT_EMF_AVERAGE                 CT_MAXWELL
#define  LIMITER                        MC_LIM
#define  RECONSTRUCT_4VEL               NO
#define  ROTATE                         0

// [End] user-defined constants (do not change this line)
