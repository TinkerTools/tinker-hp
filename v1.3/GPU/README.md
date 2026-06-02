[![DOI:10.1021/acs.jctc.0c01164](https://zenodo.org/badge/DOI/10.1021/acs.jctc.0c01164.svg)](https://doi.org/10.1021/acs.jctc.0c01164)

Tinker-HP: High Performance Multi-GPUs Massively Parallel Evolution of Tinker
=============================================================================


<b>This phase-advance GPU version (1.2 ++) is not (yet) an official release of Tinker-HP but is made freely available in link with the COVID-19 HPC community effort.</b>

This work will be part of a larger 2023 Tinker-HP 1.3 official release.
In addition to GitHub, a [GPU container](https://ngc.nvidia.com/catalog/containers/hpc:tinkerhp) (quick install!) is available thanks to NVIDIA on the NVIDIA NGC's website: 

# Getting started with Tinker-HP


## Installation Guide and tutorials
   -  [Prerequisites](Prerequisites.md)  
   Some setup may prove to be unstable. Please, find out [here](Prerequisites.md)
   -  [Build Tinker-HP (GPU version)](build.md)
   -  [Deep-HP](Deep-HP.md)
   -  [Quantum-HP](Quantum-HP_tutorial.md)


## Run Tinker-HP (CPU/GPU)
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature you are looking for is available on the GPU version. The present version is optimized to accelerate simulations using the AMOEBA polarizable force field. Some minimal non-polarizable capabilities are present (enhanced support will be available in 2021). The code has been extensively tested on 1080, 2080, 3090, P100, V100 and A100 NVIDIA GPU cards and support multi-GPUs computations. It will be part of the major Tinker-HP 1.3 2022 release but this present version will continue to evolve. 


### GPU available features
   - **Applications :** *dynamic analyze minimize bar pimd radial* programs
   - **Integrators :** *(RESPA, RESPA1, BAOAB, BAOAB-RESPA, BAOAB-RESPA1, VERLET, BEEMAN, BBK)*
   - **Force field :** polarizable ones *(Amoeba/Amoeba+)*, classical ones *(AMBER/CHARMM/OPLS)*
   - New implementation of **PCG**(default) and **DC-DIIS** solver for polarization (DC-DIIS is not adapted to the device, use PCG instead!)
   - **Thermostats :** *BUSSI* (default), *Andersen* for NVT simulations
   - **Barostats :** *Montecarlo* and *Berendsen*(default) for NPT simulations
   - Accelerate Molecular Dynamics : *aMD* and *GaMD* Simulations
   - Steered Molecular Dynamics *(SMD)*
   - **PBC Box Shapes :** *Orthogonal* and *Octahedron* (**ONLY** to be used with a single process)
   - **Plumed** support available
   - **Colvars** support available
   - **Neural Network** Interface via *Deep-HP*
   - **Lambda dynamics** support available for *free energy* calculations
   - **Nuclear Quantum Effects** via Ring Polymer Molecular Dynamics (RPMD) and adaptive Quantum Thermal Bath (adQTB)
   -  **More to come**

## Contact
For any bug report, please refere to :
   *  TinkerHP_Support@ip2ct.upmc.fr

## Please Cite
```tex
@article{2021tinkerhpGPU,
  title={Tinker-HP: Accelerating molecular dynamics simulations of large complex systems with advanced point dipole polarizable force fields using GPUs and multi-GPU systems},
  author={Adjoua, Olivier and Lagard{\`e}re, Louis and Jolly, Luc-Henri and Durocher, Arnaud and Very, Thibaut and Dupays, Isabelle and Wang, Zhi and Inizan, Th{\'e}o Jaffrelot and C{\'e}lerse, Fr{\'e}d{\'e}ric and Ren, Pengyu and others},
  journal={Journal of chemical theory and computation},
  volume={17},
  number={4},
  pages={2034--2053},
  year={2021},
  publisher={ACS Publications}
}
```