#define  PHYSICS                        RHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       SPHERICAL
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

// Physics dependent declarations

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO

// User-defined parameters (labels)

#define  ENRG0                          0
#define  DNST0                          1
#define  MASS                           2

// [Beg] user-defined constants (do not change this line)

#define  SHOCK_FLATTENING               NO
#define  RING_AVERAGE                   8
#define  LIMITER                        VAN_LEER

// [End] user-defined constants (do not change this line)
