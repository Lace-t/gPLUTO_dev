#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  P_IN                           0
#define  P_OUT                          1
#define  BMAG                           2
#define  THETA                          3
#define  PHI                            4
#define  RADIUS                         5
#define  GAMMA_EOS                      6

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  CT_EMF_AVERAGE                 UCT_HLL
#define  CT_ENERGY_CORRECTION           YES
#define  LIMITER                        VANLEER_LIM
#define  NGHOST_USR                     3
#define  RING_AVERAGE                   8
#define  RING_AVERAGE_FOURIER           NO
#define  RING_AVERAGE_REC               PARABOLIC
#define  RING_AVERAGE_REC_PRIM          YES
#define  RING_AVERAGE_CT_REC            FLAT
#define  XOFFSET                        YES
#define  WARNING_MESSAGES               3

// [End] user-defined constants (do not change this line)
