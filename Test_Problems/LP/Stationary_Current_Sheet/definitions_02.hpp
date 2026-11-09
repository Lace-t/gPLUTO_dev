#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      PARTICLES_LP
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

// Particles Options

#define  PARTICLES_LP_SPECTRA           YES
#define  PARTICLES_LP_SHOCK_SPECTRA     NO
#define  PARTICLES_LP_MREC_SPECTRA      YES

// User-defined parameters (labels)

#define  V0                             0
#define  SHEAR_WIDTH                    1
#define  ETA                            2
#define  SIGMA                          3
#define  RHO_0                          4
#define  B_GUIDE                        5
#define  B0_PARAM                       6

// [Beg] user-defined constants (do not change this line)

#define  WARNING_MESSAGES               NO
#define  CT_EMF_AVERAGE                 CT_CONTACT
#define  REC_UPDATE_DEBUG               0
#define  ACC_INIT                       NO
#define  PARTICLES_DT_CONTROL           YES
#define  PARTICLES_NCOLORS              6
#define  ANALYSYS                       YES

// [End] user-defined constants (do not change this line)
