<!-- [![DOI:10.48550/arXiv.2207.14276](https://zenodo.org/badge/DOI/10.48550/arXiv.2207.14276.svg)](https://doi.org/10.48550/arXiv.2207.14276) [![DOI:10.1039/C7SC04531J](https://zenodo.org/badge/DOI/10.1039/C7SC04531J.svg)](https://doi.org/10.1039/C7SC04531J) [![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.1039/c7sc04531j)](https://juleskreuer.eu/projekte/citation-badge/) ![GitHub followers](https://img.shields.io/github/followers/TinkerTools?style=social) ![Twitter Follow](https://img.shields.io/twitter/follow/TINKERtoolsMD?style=social)  -->



# Deep-HP: Multi-GPUs platform for hybrid Machine Learning Polarizable Potential
Deep-HP is a multi-GPU deep learning potential platform, part of the Tinker-HP molecular dynamics package and aims to couple deep learning with force fields for biological simulations. 



# What is Deep-HP?
Deep-HP aims to democratize the use of deep learning in biological simulations. <br /> 
More precisely, Deep-HP aims at scaling up **deep learning potential code from laptop to hexascale and from quantum chemistry to biophysics**. Deep-HP ultimate goal is the unification within a reactive molecular dynamics many-body interaction potential of the short-range quantum mechanical accuracy and of long-range classical effects, at force field computational cost. Application of Deep-HP span from drug-like molecule to proteins and DNA.


What can I do? Here's a few examples:
* Combine trained machine learning potential with force fields (long-range interactions, polarization effects ...)
* Compute solvation free energies of drug-like molecules.
* Compute binding free energies.
* Conformational sampling with state-of-the-art enhanced sampling techniques (Colvars, Plumed).
* ...
* Everything than before but with machine learning and faster
* For more check-out [TinkerTools](https://tinkertools.org/), [Tinker-HP](https://tinker-hp.org/)


# Dependencies

Deep-HP relies on the FeNNol library to execute ML models. Follow instructions on FeNNol's github repo (https://github.com/thomasple/FeNNol) to install it before compiling Tinker-HP. <br /> 
**Make sure to install `cffi` and `pycuda` and that the correct python environment is loaded when compiling Tinker-HP with Deep-HP**


# Run Deep-HP
Deep-HP has only three main **KEYWORDS**: `MLPOT`, `ML-MODEL` and `MLPOT-CUTOFF`. 

`MLPOT` sets the type of simulation:
* `MLPOT NONE` deactivates the machine learning potential evaluation
* `MLPOT ONLY` activates only the machine learning potential
* `MLPOT` activates the machine learning potential but also the force field. As both are evaluated at each time step don't forget to disable the terms of the forcefield you don't want to use (Example 2)
<!-- * `MLPOT EMBD` actives the machine learning potential on a group of atoms, like a QM/MM embedding. (Example 4, 5) -->

`ML-MODEL` sets the path to the machine learning potential model. This should be a FeNNol model file (typically with the .fnx extension)
  
`MLPOT-CUTOFF` sets the cutoff distance for the machine learning potential. It is a distance in Angstrom. <br />
**The cutoff distance must be the same as the one used to train your model**: this can be checked with the terminal command `fennol_inspect path/to/model.fnx` <br />


# Example
We provide an example located in the `examples` directory of Tinker-HP. This example can be run with the command:
```bash
cd examples
mpirun -n 1 ../bin/dynamic Deep-HP_example 1000 0.5 100 2 300
```
This will run a short MD simulation of a box of water using the ANI2x model stored in the `ml_models` directory.



# Please Cite
```tex
@misc{https://doi.org/10.48550/arxiv.2207.14276,
  doi = {10.48550/ARXIV.2207.14276},
  url = {https://arxiv.org/abs/2207.14276},
  author = {Jaffrelot Inizan, Théo and Plé, Thomas and Adjoua, Olivier and Ren, Pengyu and Gökcan, Hattice and Isayev, Olexandr and Lagardère, Louis and Piquemal, Jean-Philip},
  keywords = {Chemical Physics (physics.chem-ph), FOS: Physical sciences, FOS: Physical sciences},
  title = {Scalable Hybrid Deep Neural Networks/Polarizable Potentials Biomolecular Simulations including long-range effects},
  publisher = {arXiv},
  year = {2022},
  copyright = {Creative Commons Attribution 4.0 International}
}

@article{ple2024fennol,
    author = {Plé, Thomas and Adjoua, Olivier and Lagardère, Louis and Piquemal, Jean-Philip},
    title = {FeNNol: An efficient and flexible library for building force-field-enhanced neural network potentials},
    journal = {The Journal of Chemical Physics},
    volume = {161},
    number = {4},
    pages = {042502},
    year = {2024},
    month = {07},
    doi = {10.1063/5.0217688},
    url = {https://doi.org/10.1063/5.0217688},
}
```
