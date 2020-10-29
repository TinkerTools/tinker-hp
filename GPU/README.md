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
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature we are looking for is available on the GPU. The present version is optimized to accelerate simulations using the AMOEBA polarizable force field. Some minimal non-polarizable capabilities are present (enhanced support will be available in 2021).

### GPU available features
   - dynamic analyze and minimize program
   - Integrators (respa, respa1, baoab, baoabrespa1, verlet)
   - Amoeba polarizable force field, classical force fields (AMBER/CHARMM/OPLS)
   - New implementation of PCG and DC-DIIS solver for polarization (DC-DIIS is not adapted to device ! use PCG instead)
   - Bussi Thermostat for NVT simulations  (it is default)
   - Montecarlo and Berendsen barostat for NPT simulations (default is Berendsen)
   - aMD and GaMD Simulations
   - Octahedron and Orthogonal Box Shape (to be used only on single MPI process)
   
   If you use the code please cite :
   
   Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems.
O. Adjoua,  L. Lagardère, L.-H. Jolly, Arnaud Durocher, Z. Wang, T. Very, I. Dupays, T. Jaffrelot Inizan, F. Célerse, P. Ren, J. Ponder, J-P. Piquemal, 2020, preprint to come
   
   and 
   
Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access) DOI: 10.1039/C7SC04531J
