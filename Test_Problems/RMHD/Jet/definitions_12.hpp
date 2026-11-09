#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

// Physics dependent declarations

#define  EOS                            IDEAL
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

#define  ACC_INIT                       NO
#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CHECK_DIVB_CONDITION           NO
#define  COVARIANT_MX2                  NO
#define  CT_EMF_AVERAGE                 CT_CONTACT
#define  CT_EMF_BOUNDARY                YES
#define  CT_ENERGY_CORRECTION           YES
#define  LIMITER                        OSPRE_LIM
#define  RECONSTRUCT_4VEL               YES
#define  SHOCK_FLATTENING               NO
#define  WARNING_MESSAGES               3

// [End] user-defined constants (do not change this line)
