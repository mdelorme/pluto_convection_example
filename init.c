/* ///////////////////////////////////////////////////////////////////// */
/*! 
  Cartesian slab for the Whole-Sun convection benchmark.
  This setup is inspired by Cattaneo et al 1991.

  The domain is initialized using a polytropic model 
  T(z)   = (1.0 + theta*z)
  rho(z) = (1.0 + theta*z)**m
  P(z)   = (1.0 + theta*z)**(m+1)

  and closed using a ideal gas equation of state : P = rho * T

  with :
   . z the depth
   . m the polytropic index
   . theta the temperature gradient between the top and the bottom of the box

  in this setup, gravity is expressed as g = theta*(m+1) and applied along positive z. 
  It is imposed as a constant acceleration over the whole domain. The value is derived
  from hydrostatic equilibrium at initialization.

  Instability is triggered by imposing a small random perturbation on pressures at startup.

  Vertical boundary conditions :
   . Top of the domain    (z=0) -> Ttop=1.0
   . Bottom of the domain (z=1) -> delta T=theta
   . In both cases : 
     * Horizontal velocity gradients are null
     * Vertical velocity is null
     * Pressure is obtained using hydrostatic equilibrium
     * Density is recovered using the EOS with the pressure and temperature.
  Horizontal boundary conditions are taken periodic

  The domain is taken to be [0;4]x[0;4]x[0;1]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <assert.h>

int first = 1; 
const double Ttop = 1.0;

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3) {

  // Run parameters
  int    seed  = g_inputParam[SEED];
  double mpoly = g_inputParam[MPOLY];
  double theta = g_inputParam[THETA];
  double amp   = g_inputParam[AMP];

  // One-time init per process
  if (first) {
    g_gamma = 5.0/3.0;
    srand(seed*prank);
    first = 0;
  }

  // z is the vertical direction
  double z = x3;

  // Random perturbation
  double pert = (((float)rand() / RAND_MAX) - 0.5) * amp;

  // We avoid applying the perturbation at the boundaries
  if (z < 0.1 || z > 0.9)
    pert = 0.0;

  // ICs
  v[RHO] = pow(1.0 + theta*z, mpoly);
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = pow(1.0 + theta*z, mpoly+1.0) * (1.0 + pert);
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid) {}
void Analysis (const Data *d, Grid *grid) {}

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {
  /**
   * Important notes on boundary conditions definitions:
   *
   * PLUTO is a finite-volume code which means that traditionnaly we use
   * values in ghost-cells outside of the domain to impose conditions.
   * Now due to the restrictions of the run, and the complexity of flux
   * calculations, we impose values in the boundary to make sure that
   * the diffusive terms are computed correctly (thermal conduction and
   * viscosity). The actual hydro flux calculation is overwritten in 
   * update_stage.c to ensure that the mass flux is 0 through the boundaries
   *
   * To do this, we add a new flag in PLUTO : FLAG_BOUNDARY that allows us
   * to flag the specific cells that are at the boundary. This flag is
   * set in this function along with the values used for the diffusive solvers.
   * 
   * For hydro step, once we have the flag set, we detect it in update_stage.c
   * and make sure the mass flux is null through the boundary.
   *
   * NOTE : The way the pressure is defined at the top of domain can lead
   *        to negative densities when using high stratification (theta=20)
   *        with only few points (less than 128). Resolution must be pushed
   *        for runs with theta=20, or pressure must be defined in another 
   *        way.
   **/

  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  // Paramètres du run
  double theta = g_inputParam[THETA];
  double m     = g_inputParam[MPOLY];
  double gval  = theta * (m+1);
  
  if (side == X3_BEG) {  // z=0 -> Top of the domain
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        // Reference index for the cells velocities (mirror of the boundary)
        int kref = 2*KBEG - (k+1);
        
        // Reference state is the first cell in the domain
        double rho_ref = d->Vc[RHO][KBEG][j][i];
        double P_ref   = d->Vc[PRS][KBEG][j][i];
        double T_ref   = P_ref / rho_ref;

	// We compute the pressure using hydrostatic equilibrium
	// and extrapolate temperature from the value at the boundary (1.0)
        double dz      = x3[KBEG] - x3[k];
        double P       = P_ref - dz * rho_ref*gval;
        double T       = 1.0 + x3[k]*theta;

	// Density is given by the EOS
        double rho     = P / T;
          
        // Finally we fill in the ghost value and flag the boundary
        d->Vc[PRS][k][j][i] = P;
        d->Vc[RHO][k][j][i] = rho;
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][kref][j][i];
        d->Vc[VX2][k][j][i] =  d->Vc[VX2][kref][j][i];
        d->Vc[VX3][k][j][i] = -d->Vc[VX3][kref][j][i] * d->Vc[RHO][kref][j][i] / rho; // We make sure rhovz is 0 at boundary
        d->flag[KBEG-1][j][i] |= FLAG_BOUNDARY;
      }
    }
  }

  if (side == X3_END){  // z=1 -> Bottom of the domain
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        int kref = 2*KEND - (k-1);
          
        // Reference state is the first cell inside the domain
        double rho_ref = d->Vc[RHO][KEND][j][i];
        double P_ref   = d->Vc[PRS][KEND][j][i];
        double T_ref   = P_ref / rho_ref;
        double dz  = x3[k] - x3[KEND];

        // Temperature is imposed using the gradient
	// Pressure is recomputed using hydrostatic equilibrium
	// Density is given by the EOS
        double T   = T_ref + dz * theta;
        double P   = P_ref + dz * rho_ref*gval;
        double rho = P / T;

        // We fill in the ghost and flag the boundary
        d->Vc[PRS][k][j][i] = P;
        d->Vc[RHO][k][j][i] = rho;
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][kref][j][i];
        d->Vc[VX2][k][j][i] =  d->Vc[VX2][kref][j][i];
        d->Vc[VX3][k][j][i] = -d->Vc[VX3][kref][j][i] * d->Vc[RHO][kref][j][i] / rho; // We make sure rhovz is 0 at boundary

        d->flag[KEND][j][i] |= FLAG_BOUNDARY;
      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3) {

  // Gravitational acceleration is theta * (m+1)
  // We set it to 0 outside of the domain
  double sgn = 1.0;
  double z = x3;
  if (z > 1.0 || z < 0.0)
    sgn = 0.0;
  double gval = g_inputParam[THETA] * (g_inputParam[MPOLY] + 1.0) * sgn;
  
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = gval;  
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
{
  return 0.0; 
}
#endif
