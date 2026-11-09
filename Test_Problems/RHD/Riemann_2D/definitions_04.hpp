#define  PHYSICS                        RHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LimO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            1

// Physics dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO

// User-defined parameters (labels)

#define  SCRH                           0

// [Beg] user-defined constants (do not change this line)

#define  ARTIFICIAL_VISC                0.1
#define  FALLBACK_MODE                  FALLBACK_RECOMPUTE
#define  FALLBACK_TTL_MAX               4
#define  FALLBACK_BUFFER                3
#define  FALLBACK_SUBSTEPS              2
#define  CHAR_LIMITING                  YES
#define  LIMITER                        VANLEER_LIM
#define  NGHOST_USR                     3

// [End] user-defined constants (do not change this line)
