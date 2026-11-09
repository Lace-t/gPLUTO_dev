#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        2
#define  PARTICLES                      PARTICLES_LP
#define  USER_DEF_PARAMETERS            9

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

// Particles Options

#define  PARTICLES_LP_SPECTRA           YES
#define  PARTICLES_LP_SHOCK_SPECTRA     YES
#define  PARTICLES_LP_MREC_SPECTRA      NO

// User-defined parameters (labels)

#define  GAMMA_UP                       0
#define  NU                             1
#define  BETA                           2
#define  RADX                           3
#define  RADY                           4
#define  RADZ                           5
#define  SIGMA                          6
#define  TEMP                           7
#define  INJECTION_DT                   8

// [Beg] user-defined constants (do not change this line)

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CT_EMF_AVERAGE                 UCT_HLL
#define  CT_EMF_BOUNDARY                YES
#define  INTERNAL_BOUNDARY              YES
#define  INTERNAL_BOUNDARY_REFLECT      YES
#define  LIMITER                        VANLEER_LIM
#define  SHOCK_FLATTENING               MULTID
#define  WARNING_MESSAGES               3

// [End] user-defined constants (do not change this line)
