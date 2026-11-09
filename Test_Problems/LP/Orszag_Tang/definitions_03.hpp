#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENOZ
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      PARTICLES_LP
#define  USER_DEF_PARAMETERS            0

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

/* -- Particles Options -- */

#define  PARTICLES_LP_SPECTRA           YES
#define  PARTICLES_LP_SHOCK_SPECTRA                    YES
#define  PARTICLES_LP_MREC_SPECTRA                     NO

// User-defined parameters (labels)


// [Beg] user-defined constants (do not change this line)

#define  CT_EMF_AVERAGE                 UCT_HLLD
#define  CT_ENERGY_CORRECTION           NO
#define  MULTIPLE_LOG_FILES             NO
#define  DATAMINMAX                     NO

// [End] user-defined constants (do not change this line)
