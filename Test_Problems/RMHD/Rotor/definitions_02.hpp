#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
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

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  LIMITER                        VANLEER_LIM
#define  RECONSTRUCT_4VEL               YES

// [End] user-defined constants (do not change this line)
