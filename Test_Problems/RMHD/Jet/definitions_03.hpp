#define  PHYSICS                        RMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

// Physics dependent declarations

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

// User-defined parameters (labels)

#define  MACH                           0
#define  LORENTZ                        1
#define  RHOJ                           2
#define  RHOA                           3
#define  SIGMA_POL                      4
#define  SIGMA_TOR                      5

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CT_EMF_AVERAGE                 CT_CONTACT
#define  INCLUDE_JDIR                   NO
#define  LIMITER                        VANLEER_LIM

// [End] user-defined constants (do not change this line)
