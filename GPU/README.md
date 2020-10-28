Tinker-HP: High Performance Massively Parallel Evolution of Tinker
==================================================================


<H2><B>WARNING ! This GPU version is NOT an official release of Tinker-HP.
      This is the version that is fighting against the COVID-19 since the beginning of March 2020. </B></H2>
          

# Getting started with Tinker-HP
   - Installation Guide

## Installation Guide
   -  [Prerequisites](Prerequisites.md)
   -  [Build Tinker-HP (GPU version)](build.md)

## Run Tinker-HP (CPU/GPU)
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature we are looking for is available on the GPU.

### GPU available features
   - dynamic analyze and minimize program
   - Integrators (respa, respa1, baoab, baoabrespa1, verlet)
   - Amoeba forcefield
   - PCG and DC-DIIS solver for polarisation (DC-DIIS is not adapted to device ! use PCG instead)
   - Bussi Thermostat for NVT simulations  (it is default)
   - Montecarlo and Berendsen barostat for NPT simulations (default is Berendsen)
   - aMD and GaMD Simulations
   - Octahedron and Orthogonal Box Shape (to be used only on single MPI process)
