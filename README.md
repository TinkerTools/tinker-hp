# Tinker-HP: High-Performance Massively Parallel Evolution of Tinker on CPUs & GPUs
-----------------------------------------------------------------------------------------------------------------------------------------------

<H2><B>News</B></H2>

------------------------------
<B>-Update 10/2022: New website for Tinker-HP, check it out!!! https://tinker-hp.org

<B>-Update 09/2022: support for neural networks potentials (ANI-2X, DeepMD etc...) is available in beta version (see Deep-HP branch on Github). Check the Deep-HP module paper: https://arxiv.org/abs/2207.14276</B>

<B>-Update 03/02/2021 : PLUMED Support for version 1.2 GPUs !</B>

<B>-Update 02/11/2021 : PLUMED Support for version 1.2 (CPUs) !</B>

<B>-Update 11/24/2020 : all versions have been pushed to GitHub</B>
 
 ------------------------------
 
 <H2><B>Versions</B></H2>

Current Github version: 1.1v (enhanced AVX512 vectorized CPUs version), 1.2 (CPUs) + 1.2 (multi)-GPUs

Current Development version: 1.3 (CPUs + multi-GPUs)


All releases of the Tinker-HP code are now being performed on Github. For news, benchmarks and additional tutorials, please visit the Tinker-HP website
http://tinker-hp.org/  
 
<B>In addition to GitHub, a GPUs container (quick install!) is available thanks to NVIDIA on the NVIDIA NGC's website: https://ngc.nvidia.com/catalog/containers/hpc:tinkerhp </B>
 
 ------------------------------

Tinker-HP is a CPUs and GPUs based, multi-precision, MPI massively parallel package dedicated to long polarizable molecular dynamics simulations and to polarizable QM/MM. Tinker-HP is an evolution of the popular Tinker package that conserves it simplicity of use but brings new 
capabilities allowing performing very long molecular dynamics simulations on modern supercomputers that use thousands of cores. 
The Tinker-HP approach offers various strategies using domain decomposition techniques for periodic boundary conditions in the 
framework of the (n)log(n) Smooth Particle Mesh Ewald. Tinker-HP proposes a high performance scalable computing environment for 
polarizable (AMOEBA, Amberpol...) and classical (Amber, Charmm, OPLS...) force fields giving access to large systems up to millions of atoms. It can be used on supercomputers as well as on lab clusters. Tinker-HP supports Intel (AVX5212 enhanced version) and AMD CPUs platforms as well as NVIDIA GPUs (1080, 2080, 3090, P100, V100, A100). 

i) Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories and supercomputer centers through the global Tinker license (https://dasher.wustl.edu/tinker/downloads/license.pdf).
Non-academic entities (e.g., companies, for profit organizations) should contact the managing universities (see license).

ii) <B>If you use Tinker-HP, please cite:</B>

- Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems 
with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, 
M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access)
https://doi.org/10.1039/C7SC04531J

- if you are using the GPUs version, please also cite:
Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems
Olivier Adjoua, Louis Lagardère, Luc-Henri Jolly, Arnaud Durocher, Thibaut Very, Isabelle Dupays, Zhi Wang, Théo Jaffrelot Inizan, Frédéric Célerse, Pengyu Ren, Jay W. Ponder, Jean-Philip Piquemal, J. Chem. Theory. Comput., 2021, 17 (4), 2034–2053(Open Access) https://doi.org/10.1021/acs.jctc.0c01164

- if you are using the AVX512 vectorized version dedicated to Intel's CPUs (Skylake, CascadeLake etc...), please also cite:
Raising the Performance of the Tinker-HP Molecular Modeling Package [Article v1.0].
L. H. Jolly, A. Duran, L. Lagardère, J. W. Ponder, P. Y. Ren, J.-P. Piquemal, LiveCoMS, 2019, 1 (2), 10409  (Open Access)
 https://doi.org/10.33011/livecoms.1.2.10409
 
iii) <b>Tinkertools</b>
 
Tinker-HP is part of the Tinker distribution and uses the same tools as Tinker. These tools can be found here : https://github.com/TinkerTools/tinker

If you use the Tinkertools please cite :

Tinker 8: Software Tools for Molecular Design.
J. A. Rackers, Z. Wang, C. Lu, M. L. Maury, L. Lagardère, M. J. Schnieders, J.-P. Piquemal, P. Ren, J. W. Ponder,  J. Chem. Theory. Comput., 2018, 14 (10), 5273–5289 DOI: http://dx.doi.org/10.1021/acs.jctc.8b00529 ; PMC free text : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6335969/

iv) <B>Support:</B>

We provide support to users:

Email: TinkerHP_Support@ip2ct.upmc.fr


