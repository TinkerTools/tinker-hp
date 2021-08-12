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
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature you are looking for is available on the GPU version. The present version is optimized to accelerate simulations using the AMOEBA polarizable force field. Some minimal non-polarizable capabilities are present (enhanced support will be available in 2021). The code has been extensively tested on 1080, 2080, 3090, P100, V100 and A100 NVIDIA GPU cards and support multi-GPUs computations. It will be part of the major Tinker-HP 1.3 2021 release but this present version will continue to evolve. 

### GPU available features
   - dynamic analyze minimize and bar programs
   - Integrators (RESPA, RESPA1, BAOAB, BAOAB-RESPA1, VERLET)
   - Amoeba polarizable force field, classical force fields (AMBER/CHARMM/OPLS)
   - New implementation of PCG and DC-DIIS solver for polarization (DC-DIIS is not adapted to the device, use PCG instead)
   - Bussi Thermostat for NVT simulations  (it is default)
   - Montecarlo and Berendsen barostat for NPT simulations (default is Berendsen)
   - Accelerate Molecular Dynamics : aMD and GaMD Simulations
   - Steered Molecular Dynamics (SMD)
   - Orthogonal and Octahedron PBC Box Shapes (the latest to be used only on a single MPI process for now)
   - Plumed support available (updated : 03/2021)
   - More to come

   For more detailed informations on how to use the application, see section V to VIII of the readme of the CPU 1.2 version: https://github.com/TinkerTools/tinker-hp/blob/master/v1.2/Readme_v1.2.pdf
   Beware that not all the features available on CPU are available on GPU (see above, for example the TCG solver is only available on CPU as is the Langevin Piston barostat).
   
   <B>If you use the code please cite :</B>
   
   Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems.
O. Adjoua,  L. Lagardère, L.-H. Jolly, Arnaud Durocher, Z. Wang, T. Very, I. Dupays, T. Jaffrelot Inizan, F. Célerse, P. Ren, J. Ponder, J-P. Piquemal, J. Chem. Theory. Comput., 2021, 17 (4), 2034–2053 (Open Access) https://doi.org/10.1021/acs.jctc.0c01164
   
   and 
   
Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access) https://doi.org/10.1039/C7SC04531J

<B>License :</B> 

Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories and supercomputer centers through the global Tinker license (https://dasher.wustl.edu/tinker/downloads/license.pdf). Non-academic entities (e.g., companies, for profit organizations) should contact the managing universities (see license).

<B>Tinkertools :</B> Tinker-HP is part of the Tinker distribution and uses the same tools as Tinker. These tools can be found here : https://github.com/TinkerTools/tinker if you use the Tinkertools please cite :

Tinker 8: Software Tools for Molecular Design. J. A. Rackers, Z. Wang, C. Lu, M. L. Maury, L. Lagardère, M. J. Schnieders, J.-P. Piquemal, P. Ren, J. W. Ponder, J. Chem. Theory. Comput., 2018, 14 (10), 5273–5289 DOI: http://dx.doi.org/10.1021/acs.jctc.8b00529 PMC free text : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6335969/

<B>Support :</B>

We provide support to registered users only (http://tinker-hp.ip2ct.upmc.fr/?Download-instructions).

Email: TinkerHP_Support@ip2ct.upmc.fr

<B>Funding :</B> 
- this work has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 810367), project EMC2 (see preprint for full acknowledgments)

- we thank GENCI, NVIDIA and HPE as well as the engineering team of the IDRIS Supercomputer center (CNRS/GENCI, France). 
