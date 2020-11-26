#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     VECTOR
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          YES
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             RK_LEGENDRE
#define  VISCOSITY                      RK_LEGENDRE
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  THETA                          0
#define  CK                             1
#define  SIGMA                          2
#define  GAMMA                          3
#define  PERT                           4
#define  MPOLY                          5
#define  SEED                           6

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES

/* [End] user-defined constants (do not change this line) */
