#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             SUPER_TIME_STEPPING
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  THETA                          0
#define  MPOLY                          1
#define  AMP                            2
#define  SEED                           3
#define  CK                             4
#define  SIGMA                          5

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   1.0
#define  UNIT_LENGTH                    1.0
#define  UNIT_VELOCITY                  1.0

// New flag for boundary definition inside the run
#define  FLAG_BOUNDARY                  128 

/* [End] user-defined constants (do not change this line) */
