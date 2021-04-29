# HMOranges
M2Oranges builds off [HMOranges](https://github.com/TheOfficialDarthVader/HMOranges) by adding a the Abramzon-Sirignano droplet model (M2). The implementation is as per HMOranges, with extra variables and a changes to the kernel to include the new droplet model.

## Project Structure
The project is organised into relevant directories. The contents of each of these is described below.

#### analysis/
This directory contains all types of the analysis and a Python script for post-simulation data processing.
The Python script `Stats.py` which is used to collate and graph the data from simulations.
The subdirectory `stokes_numbers` is inherited from HMOranges but has not been updated to use with M2 yet.

#### kernels/
This directory contains all of the OpenCL kernels for the project.
`iterate_particle_m2.cl` is an outdated version of the M2 kernel as it is now merged into `iterate_particle.cl`.
`iterate_particle_tgv.cl` is the kernel specific for `analysis/tgvmag/tgv_sim.c`.
When running the project this directory must be in the same directory as the working directory so that the program can access the kernel files.
Note, the kernel utilities file is in `util`, not `kernels`.

####  sims/
This directory contains all simulation programs as well as the main `simRunner` files.
`simRunner` is accessed by the simulations and handles all of the simulation running code.
`simRunner_tgv` has the same function as `simRunner` but is specified for `analysis/tgvmag/tgv_sim.c` only.
A variety of simulations are available and are all CMake targets.

#### structures/
This directory contains all of the particle data structures used in the host code.
Matching device structures are found in `util/kernelUtils.cl`.

#### tests/
This directory contains all of the unit tests for the project.
These are accessed from `test/run_tests/` in simulations or can be run separately with more debug messages with the `run_tests` target.

#### util/
This directory contains all of the utility functions used in the project.
Most of these are found in `util/clUtils/` which contains all of the OpenCL utility functions that make calling OpenCL functions much simpler.
This directory also contains `kernelUtils.cl` which contains all of the utility functions and data structures for the OpenCL kernels.

#### verification/
This directory contains all of the verification test cases.
Actual simulations can be run using their targets and the resulting data can be graphed against the analytic solutions with the Python script in each directory.
