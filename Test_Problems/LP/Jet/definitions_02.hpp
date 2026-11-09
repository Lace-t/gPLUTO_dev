#define  PHYSICS                        RMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        2
#define  PARTICLES                      PARTICLES_LP
#define  USER_DEF_PARAMETERS            13

// Physics dependent declarations

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

// Particles Options

#define  PARTICLES_LP_SPECTRA           YES
#define  PARTICLES_LP_SHOCK_SPECTRA     YES
#define  PARTICLES_LP_MREC_SPECTRA      NO

// User-defined parameters (labels)

#define  JET_GAMMA                      0
#define  JET_ETA                        1
#define  JET_PRS                        2
#define  JET_RAD                        3
#define  JET_SIGMAB                     4
#define  JET_ANG                        5
#define  JET_TRC                        6
#define  JET_TON                        7
#define  JET_TOFF                       8
#define  T_RAMP_UP                      9
#define  T_RAMP_DN                      10
#define  INJECT_FREQ                    11
#define  HALO_N0                        12

// [Beg] user-defined constants (do not change this line)

#define  SHOCK_FLATTENING               MULTID
#define  EPS_PSHOCK_FLATTEN             5
#define  CT_EMF_AVERAGE                 CT_CONTACT
#define  CT_ENERGY_CORRECTION           YES
#define  INTERNAL_BOUNDARY              YES
#define  DEBUG_MIN_MAX                  YES
#define  PARTICLES_LP_SHOCK_THRESHOLD   2
#define  PARTICLES_LP_SHOCK_GRADP       2
#define  DIV_JET                        NO
#define  PERTURB_INLET                  NO
#define  UNIT_LENGTH                    (CONST_pc)
#define  WARNING_MESSAGES               2

// [End] user-defined constants (do not change this line)
