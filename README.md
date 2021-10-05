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
 * `update_stage.c`: the update algorithm. The file is modified to force boundary mass flux to zero.
 * `plot_run.py`: the plotting and info extraction script (see below)

## How to compile and run
Provided you are working with a UNIX-compliant system, that PLUTO is installed and the `$PLUTO_DIR` variable is set up, the run is setup using the pluto script :
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
Since the last version of this repo (feb. 2021) substantial changes have been made. In particular in the treatment of boundary conditions. Before the update, we used to define the boundary conditions in the first cell of the domain (at top and bottom boundaries) to avoid mass loss due to the lack of control in the results of the Riemann solver. Since we have adopted a different strategy where the flux is set manually in `update_stage.c`. The values in the ghost layers are still set for the diffusive kernels (thermal conduction and viscosity) to take place. This allows us finer control over what is happening at the boundary and hence perfect mass conservation.

## Plotting and data extraction
We have also provided in the repo a python script `plot_run.py` which extracts all the information necessary for the plotting of the runs. The script undergoes the following steps :

1- Rendering the simulation. The results are stored as a series of png files in the folder `render`. 
2- Extracting the time evolution of the simulation. This will result in the creation of the `pluto_time.csv` file in the format required for the benchmark as well as the creation of a `time_evolution.png` image with the evolution in time of the mass as well as the kinetic, internal and total energies.
3- Extrating the time evolution of temperature profiles. This will result in the creation of a `temperatures.png` file where each line corresponds to the temperature profile of a snapshot. By default we draw one line every 50 snapshots (hence 20 lines for 1000 snapshots). All the lines should be starting at `z=0` from the top temperature (`Ttop=1.0`) and should end at `z=1` with the same gradient.
4- Finally, the averaged fluxes and profiles for the analysis. These are horizontally averaged and time averaged between the times t=0.895 and t=0.905. This step writes a `pluto_prof.csv` file with the fluxes formatted correctly for the benchmark analysis.

Each step can be deactivated by using the following command line arguments after calling the python script : `--no-render` (step 1); `--no-time-evolution` (step 2); `--no-temperatures` (step 3) and `--no-profiles` (step 4).
