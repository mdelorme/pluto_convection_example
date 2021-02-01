# Example Pluto run for the Whole Sun convection benchmark 
This repository holds an example run for the Whole Sun convection benchmark. This example is meant to be compiled and run with the [PLUTO](http://plutocode.ph.unito.it/) code. This example corresponds to runs #3 of the benchmark with theta=10, sigma=0.1. 

Important note : This example is compatible with version 4.3 of PLUTO and has not yet been adapted for 4.4 !

## Files description
The files in the repository are the following : 
 * `init.c`: the initial conditions, the boundary conditions and the body-force definitions
 * `definitions.h`: definitions for the run, including user parameters
 * `pluto.ini`: the input configuration file for the run
 * `tc_kappa.c`: definition of the Thermal conduction parameter (kappa)
 * `visc_nu.c`: definition of the viscosity parameter (nu)

## How to compile and run
Provided you are working with a UNIX-compliant system, that PLUTO is installed and the `$PLUTO_DIR` variable set up, the run is setup using the pluto script :
```bash
python $PLUTO_DIR/setup.py
```
The file `definitions.h` should already have every parameter setup for the run, so press enter until reaching the makefile configuration. Pick a makefile adapted to your environment and exit the configuration tool. The code is then built using make 

```bash
make
```

and finally run, either in serial mode:

```bash
./pluto
```

or in parallel using mpi (if you picked a mpi-compilation makefile in the configuration tool):

```bash
mpirun -np 8 ./pluto
```

## Remarks
The configuration file `pluto.ini` defines the upper and lower boundaries as outflows. This is because PLUTO is a finite-volume code where boundary conditions are imposed using ghost cells. It is very difficult to get an impenetrable wall with a specific value for the temperature using this method because the value at the actual boundary is found after solving the Riemann problem and will actually depend on the Riemann solver user, a discussion of this problem can be found in (Freytag 2012). To solve this problem, instead of using the ghost layers, we setup an internal boundary and define the values of the first and last layer of cells along the z (vertical) direction. This can be seen in `init.c` at lines 156-189 
