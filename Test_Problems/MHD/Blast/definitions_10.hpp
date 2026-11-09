#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENOZ
#define  TIME_STEPPING                  RK4
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

#define  CHAR_LIMITING                  YES
#define  CT_EMF_AVERAGE                 UCT_HLL
#define  CT_ENERGY_CORRECTION           YES
#define  LIMITER                        MINMOD_LIM
#define  HO_LAP_LIMITER                 HO_RDER_LIM
#define  HO_LAP_LIMITER_EPS             0.1
#define  HO_ORDER_REDUCTION             WENO3

// [End] user-defined constants (do not change this line)
