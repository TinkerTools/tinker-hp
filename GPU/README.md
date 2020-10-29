Tinker-HP: High Performance Multi-GPUs Massively Parallel Evolution of Tinker
==================================================================


<H2><B>WARNING ! </b></h2>   <b>This phase-advance GPU version is NOT an official release of Tinker-HP.</b>

But this is the version that is currently fighting against the COVID-19 outbreak since the beginning of March 2020.
We decided to make it freely available to the scientific community.

This work will be part of a larger 2021 Tinker-HP 1.3 official release.

# Getting started with Tinker-HP
   - Installation Guide

## Installation Guide
   -  [Prerequisites](Prerequisites.md)
   -  [Build Tinker-HP (GPU version)](build.md)

## Run Tinker-HP (CPU/GPU)
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature we are looking for is available on the GPU. The present version is optimized to accelerate simulations using the AMOEBA polarizable force field. Some minimal non-polarizable capabilities are present (enhanced support will be available in 2021). The code has been extensively tested on 1080, 2080, V100 and A100 NVIDIA GPU cards and support multi-GPUs computations. It will be part of the major Tinker-HP 1.3 2021 release but this present version will continue to evolve. Feedbacks are welcomed : 

### GPU available features
   - dynamic analyze and minimize program
   - Integrators (RESPA, RESPA1, BAOAB, BAOAB-RESPA1, VERLET)
   - Amoeba polarizable force field, classical force fields (AMBER/CHARMM/OPLS)
   - New implementation of PCG and DC-DIIS solver for polarization (DC-DIIS is not adapted to device ! use PCG instead)
   - Bussi Thermostat for NVT simulations  (it is default)
   - Montecarlo and Berendsen barostat for NPT simulations (default is Berendsen)
   - Accelerate Molecular Dynamics : aMD and GaMD Simulations
   - Steered Molecuar Dynamcis (SMD)
   - Orthogonal and Octahedron PBC Box Shapes (the latest to be used only on a single MPI process)
   - Plumed support (to appear)
   - More to come
   
   If you use the code please cite :
   
   Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems.
O. Adjoua,  L. Lagardère, L.-H. Jolly, Arnaud Durocher, Z. Wang, T. Very, I. Dupays, T. Jaffrelot Inizan, F. Célerse, P. Ren, J. Ponder, J-P. Piquemal, 2020, preprint to come
   
   and 
   
Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access) DOI: 10.1039/C7SC04531J

License : Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories and supercomputer centers through the global Tinker license (https://dasher.wustl.edu/tinker/downloads/license.pdf). Non-academic entities (e.g., companies, for profit organizations) should contact the managing universities (see license).

Funding : this work has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 810367), project EMC2 (see preprint for full acknowledgments)
