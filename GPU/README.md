Tinker-HP: High Performance Multi-GPUs Massively Parallel Evolution of Tinker
==================================================================


<H2><B>WARNING ! </b></h2>   <b>This GPU version is NOT an official release of Tinker-HP.</b>

      But this is the version that is fighting against the COVID-19 since the beginning of March 2020.
      We decided to make it freely available to the scientific community. The code has been extensively tested on 2080, V100 and A100 NVIDIA cards. It will be part of the major Tinker-HP 1.3 2021 release.

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
