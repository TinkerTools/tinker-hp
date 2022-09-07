Tinker-HP: High Performance Multi-GPUs Massively Parallel Evolution of Tinker
==================================================================


<b>This phase-advance GPU version (1.2 ++) is not (yet) an official release of Tinker-HP but is made freely available in link with the COVID-19 HPC community effort.</b>

This work will be part of a larger 2021 Tinker-HP 1.3 official release.
In addition to GitHub, a GPUs container (quick install!) is available thanks to NVIDIA on the NVIDIA NGC's website: https://ngc.nvidia.com/catalog/containers/hpc:tinkerhp

# Getting started with Tinker-HP
   - Installation Guide

## Installation Guide
   -  [Prerequisites](Prerequisites.md)
   -  [Build Tinker-HP (GPU version)](build.md)

## Run Tinker-HP (CPU/GPU)
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature you are looking for is available on the GPU version. The present version is optimized to accelerate simulations using the AMOEBA polarizable force field. Some minimal non-polarizable capabilities are present (enhanced support will be available in 2021). The code has been extensively tested on 1080, 2080, 3090, P100, V100 and A100 NVIDIA GPU cards and support multi-GPUs computations. It will be part of the major Tinker-HP 1.3 2022 release but this present version will continue to evolve. 

### GPU available features
   - dynamic analyze minimize and bar programs
   - Integrators (RESPA, RESPA1, BAOAB-RESPA, BAOAB-RESPA1, VERLET)
   - Amoeba polarizable force field, classical force fields (AMBER/CHARMM/OPLS)
   - New implementation of PCG and DC-DIIS solver for polarization (DC-DIIS is not adapted to the device, use PCG instead)
   - Bussi Thermostat for NVT simulations  (it is default)
   - Montecarlo and Berendsen barostat for NPT simulations (default is Berendsen)
   - Accelerate Molecular Dynamics : aMD and GaMD Simulations
   - Steered Molecular Dynamics (SMD)
   - Orthogonal and Octahedron PBC Box Shapes (the latest to be used only on a single MPI process for now)
   - Plumed support available
   - Colvars support available
   - Lambda dynamic support available
   - Neural Network Potentials (ANI-2X, etc...) support available in beta version (Checkout the Deep-HP branch)
   - More to come
