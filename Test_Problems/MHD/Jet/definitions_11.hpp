#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        SNEq
#define  RECONSTRUCTION                 LimO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  ETA                            0
#define  JET_VEL                        1
#define  SIGMA_Z                        2
#define  SIGMA_PHI                      3
#define  PERT_AMPLITUDE                 4
#define  PERT_PERIOD                    5

// [Beg] user-defined constants (do not change this line)

#define  INCLUDE_JDIR                   NO
#define  COOLING_TIME_STEPPING          SELECTIVE
#define  ROS34_MAXSTEPS                 40
#define  WARNING_MESSAGES               2

// [End] user-defined constants (do not change this line)
