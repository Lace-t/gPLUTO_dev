#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CYLINDRICAL
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
#define  BACKGROUND_FIELD               YES
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  VEL_0                          0

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  LIMITER                        MC_LIM

// [End] user-defined constants (do not change this line)
