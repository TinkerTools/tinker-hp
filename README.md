# Tinker-HP: High-Performance Massively Parallel Evolution of Tinker on CPUs & GPUs {#mainpage}

## News {#news}
---
* **Update 01/2026:** Further speedups for the FeNNix-Bio1 foundation machine learning model via Multiple Time Steps and Distillation (DMTS). Check the [DMTS paper](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.5c03720) (J. Phys. Chem. Lett. 2026, DOI: 10.1021/acs.jpclett.5c03720)
* **Update 05/2025:** Integration of the FeNNix-Bio1 foundation machine learning model for molecular dynamics simulations. Check the [FeNNix-Bio1 paper](https://doi.org/10.26434/chemrxiv-2025-f1hgn-v4) (ChemRxiv)
* **Update 07/2024:** Integration of the Lambda-ABF method for alchemical free energy simulations. Check the [Lambda-ABF paper](https://doi.org/10.1021/acs.jctc.3c01249) (J. Chem. Theory Comput. 2024, 20, 11, 4481–4498)
* **Update 08/2023:** Integration of Quantum-HP for Nuclear Quantum Effects (NQE) through RPMD and (ad)QTB methods. Check the [Quantum-HP paper](https://doi.org/10.1021/acs.jctc.2c01233) (J. Chem. Theory Comput. 2023, 19, 5, 1432–1445)
* **Update 02/2023:** Support for neural networks potentials (ANI-2X, DeepMD etc...) is available. Check the [Deep-HP module paper](https://doi.org/10.1039/D2SC04815A) (Chem. Sci., 2023,14, 5438-5452)
* **Update 02/2023:** Support for the **AMOEBA+** potential now available.
* **Update 10/2022:** **New website for Tinker-HP**, check it out! https://tinker-hp.org
* **Update 02/2021:** **PLUMED** Support for version 1.2 GPUs
* **Update 11/2021:** **PLUMED** Support for version 1.2 (CPUs)
* **Update 24/2020:** All versions have been pushed to GitHub.

---

## Versions {#versions}

* **Current Github version:** [1.3 (CPUs)](v1.3/CPU) + [1.3 (multi)-GPUs](v1.3/GPU), [1.1v (enhanced AVX512 vectorized CPUs version)](v1.1v/)

All releases of the Tinker-HP code are now being performed on Github. For news, benchmarks, and additional tutorials, please visit the [Tinker-HP website](https://tinker-hp.org/) and follow us on [Twitter](https://twitter.com/TINKERtoolsMD).

In addition to GitHub, a GPUs container (quick install!) is available thanks to NVIDIA on the [NVIDIA NGC's website.](https://ngc.nvidia.com/catalog/containers/hpc:tinkerhp)

---

## Description {#description}

**Tinker-HP** is a **CPUs and GPUs** based, multi-precision, **MPI** massively parallel package dedicated to long **polarizable molecular dynamics** simulations and to polarizable **QM/MM**. Tinker-HP is an evolution of the popular Tinker package that conserves its simplicity of use but brings new capabilities allowing performing very long molecular dynamics simulations on modern supercomputers that use thousands of cores. 

The Tinker-HP approach offers various strategies using domain decomposition techniques for periodic boundary conditions in the framework of the *(n)log(n) Smooth Particle Mesh Ewald*. 



Tinker-HP proposes a high-performance scalable computing environment for polarizable **(AMOEBA, AMOEBA+, HIPPO...)** and classical **(Amber, Charmm, OPLS...) force fields** giving access to large systems up to **millions of atoms**. It can be used on supercomputers as well as on lab clusters. Tinker-HP supports **Intel** (**AVX-512** enhanced version) and AMD CPUs platforms as well as **NVIDIA GPUs** *(GTX-10xx, RTX-20xx, 30xx, 40xx, P100, V100, A100)*. 

---

## Documentation

Various pages documenting the capabilities of the package as well as guidelines to use it are given in the [doxygen documentation](html/pages.html), you will find:
  - [an overview of the suite detailing the programs and their cabapibilites](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_binaries_doxygen.html)
  - [prerequisites and a guide to compile the CPU version](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_Build_CPU.html)
  - [prerequisites for the GPU version](html/md_Prerequisites_doxygen.html) and a [guide to build it](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_Build_GPU.html)
  - [detailed information about the dynamic program to run molecular dynamics](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_dynamic_doxygen.html)
  - [a list of the potential energy functions available](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_potential_doxygen.html)
  - [information about I/O](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_io_doxygen.html)
  - [general information about free energy calculations](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_free_energy.html)
  - [detailed information about the Lambda-abf method for free energy simulation](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_lambda-abf_doxygen.html)
  - [a guide to use the Deep-HP interface to use Machine Learning Potentials](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_Deep-HP.html)
  - [a guide to use Quantum-HP to include Nuclear Quantum Effects in molecular dynamics](https://htmlpreview.github.io/?https://github.com/TinkerTools/tinker-hp/blob/master/html/md_Quantum-HP_doxygen.html)

A complete description of the sources following doxygen standard can also be found.

## Licence {#licence}

Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories, and supercomputer centers through the global [Tinker license](https://dasher.wustl.edu/tinker/downloads/license.pdf).  
Non-academic entities (e.g., companies, for-profit organizations) should contact the managing universities (see [license](license-Tinker.pdf)).

---

## Please Cite {#citation}

* If you use **Tinker-HP**, please cite:  
  [Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems with Advanced Polarizable Force Fields. L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access)](https://doi.org/10.1039/C7SC04531J)

* If you use the **GPUs version**, please also cite:  
  [Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems Olivier Adjoua, Louis Lagardère, Luc-Henri Jolly, Arnaud Durocher, Thibaut Very, Isabelle Dupays, Zhi Wang, Théo Jaffrelot Inizan, Frédéric Célerse, Pengyu Ren, Jay W. Ponder, Jean-Philip Piquemal, J. Chem. Theory. Comput., 2021, 17 (4), 2034–2053 (Open Access)](https://doi.org/10.1021/acs.jctc.0c01164)

* For the **AVX512 vectorized version** dedicated to Intel's CPUs (Skylake, CascadeLake etc...), please also cite:  
  [Raising the Performance of the Tinker-HP Molecular Modeling Package [Article v1.0]. L. H. Jolly, A. Duran, L. Lagardère, J. W. Ponder, P. Y. Ren, J.-P. Piquemal, LiveCoMS, 2019, 1 (2), 10409  (Open Access)](https://doi.org/10.33011/livecoms.1.2.10409)

Tinker-HP is part of the Tinker distribution and uses the same tools as Tinker. These tools can be found [here](https://github.com/TinkerTools/tinker).

* If you use the Tinkertools please cite:  
  [Tinker 8: Software Tools for Molecular Design. J. A. Rackers, Z. Wang, C. Lu, M. L. Maury, L. Lagardère, M. J. Schnieders, J.-P. Piquemal, P. Ren, J. W. Ponder,  J. Chem. Theory. Comput., 2018, 14 (10), 5273–5289](http://dx.doi.org/10.1021/acs.jctc.8b00529)  
  PMC free text: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6335969/

---

## Contact {#contact}

We provide support to users:

**Email:** TinkerHP_Support@ip2ct.upmc.fr
