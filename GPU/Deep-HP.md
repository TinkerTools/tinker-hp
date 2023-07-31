[![DOI:10.48550/arXiv.2207.14276](https://zenodo.org/badge/DOI/10.48550/arXiv.2207.14276.svg)](https://doi.org/10.48550/arXiv.2207.14276) [![DOI:10.1039/C7SC04531J](https://zenodo.org/badge/DOI/10.1039/C7SC04531J.svg)](https://doi.org/10.1039/C7SC04531J) [![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.1039/c7sc04531j)](https://juleskreuer.eu/projekte/citation-badge/) ![GitHub followers](https://img.shields.io/github/followers/TinkerTools?style=social) ![Twitter Follow](https://img.shields.io/twitter/follow/TINKERtoolsMD?style=social) 



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



# Key features
* **Compatible with TensorFlow and Pytorch**, among the most popular and efficient machine learning libraries that encompass almost all possible deep learning architectures.  
* **Easy combination with TorchANI and DeePMD models**  
* **Performing highly efficient classical molecular dynamics with your deep learning model** thanks to the highly optimized GPU-resident Tinker-HP code. Tinker-HP and Tinker9 are the fastest code for many-body polarizable force fields.
* **And also ... quantum molecular dynamics** with our new extension, coming soon, with path-integral, adaptive quantum thermal bath, ...
* **Combine molecular dynamics with state-of-the-art enhanced sampling techniques** through the recent coupling with PLUMED and Colvars.



# Run Deep-HP
Deep-HP has only two **KEYWORDS**: `MLPOT` and `ML-MODEL`. 

`MLPOT` sets the type of simulation:
* `MLPOT NONE` deactivates the machine learning potential evaluation
* `MLPOT ONLY` activates only the machine learning potential (Example 1, 3, 6)
* `MLPOT` activates the machine learning potential but also the force field. As both are evaluated at each time step don't forget to disable the terms of the forcefield you don't want to use (Example 2)
* `MLPOT EMBEDDING` actives the machine learning potential on a group of atoms, like a QM/MM embedding. (Example 4, 5)

`ML-MODEL` sets the machine learning model. For Tinker-HP native machine learning model, [ANI1X](https://pubs.rsc.org/en/content/articlehtml/2017/sc/c6sc05720a), [ANI1CCX](https://www.nature.com/articles/s41467-019-10827-4), [ANI2X](https://pubs.acs.org/doi/full/10.1021/acs.jctc.0c00121) and [ML_MBD](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.2c00936), **KEYWORDS** are explicit:
* `ML-MODEL ANI2X` uses ANI2X as machine learning potential. Same for ANI1X, ANI1CXX and ML_MBD. (Example 1, 2, 4, 5, 6)

For non native machine learning models, things work similarly. You should put your machine learning potential model inside your simulation directory (with your .xyz and .key files) and write explicity the name, including extension, of your model. <br />
The begining of the **KEYWORDS** depends of whether the model was built with TorchANI or [DeePMD](https://docs.deepmodeling.com/projects/deepmd/en/master/index.html):
* `ML-MODEL ANI_GENERIC my_mlp_model.pt # my_mlp_model.json my_mlp_model.pkl my_mlp_model.yaml`
* `ML-MODEL DEEPMD my_mlp_model.pb` (Example 3)

:warning: Almost all Tinker-HP features are compatible with Deep-HP such as RESPA and RESPA1 integrators **BUT** the splitting used is a bit different from AMOEBA and is **only available when using the embedding scheme** `MLPOT EMBEDDING`. To understand why, have a look in section 2.3.5 of [our paper](https://arxiv.org/pdf/2207.14276.pdf)


## TorchANI format
If your `model` was trained with TorchANI, you don't have to do anything in particular but make sure that :warning: **your predicted energies is in Hartree and atomic coordinates is in Angstrom** :warning: (default units in TorchANI). Tinker-HP will convert it to kcal/mol. <br /> 
For more informations have a look at [TorchANI](https://aiqm.github.io/torchani/examples/energy_force.html) or directly in our source code extension located in your tinkerml environment `/home/user/.../anaconda3/envs/tinkerml/lib/python3.9/site-packages/torchani`. <br />
Our TorchANI extension has also functions that convert your `model` in `pkl`, `yaml` and `json` formats. These formats are more compact and readable than TorchANI's `jit` format. <br />
But more importantly, when you are saving a `model` in `jit` format you are saving the whole code and you will not be able to use Tinker-HP's full capabilities for example its highly efficient neighbor list and thus use multi-GPU, Particle Mesh Ewald,... <br />
We recommend to save your model in `pkl`, `yaml` or `json` formats. Here is a python code that explain how to do it, starting from this [Example](https://aiqm.github.io/torchani/examples/nnp_training.html) of TorchANI:

```python
# ...
# once you have trained your model:
model = torchani.nn.Sequential(aev_computer, nn).to(device)
# ...
_, predicted_energies = model((species, coordinates))
# ...
# you can save it with:
model.to_json("my_mlp_model.json")
model.to_pickle("my_mlp_model.pkl")
model.to_yaml("my_mlp_model.yaml")
```

## DeePMD
For DeePMD it is similar but we don't provide other formats than the original `pb`. <br />
:warning: **In DeePMD the predicted energies are in eV and atomic coordinates are in Angstrom** (default units in DeePMD). Tinker-HP will convert it to kcal/mol.


# Example
We provide 6 examples that encompass the basics of Deep-HP inputs with which you can cover a wide range of applications. They are located in `/home/user/.../tinker-hp/GPU/examples/`. Some toy machine learning potential models are located in `/home/user/.../tinker-hp/GPU/ml_models/`.

* **Example 1:** <br />
:green_circle: *Objective*: Perform machine learning potential simulation - on the full system. <br />
:large_blue_circle: Simulation parameter: NPT with montecarlo barostat, velocity-verlet integrator and ANI2X potential. <br />
<ins>Command</ins>: :arrow_right: `mpirun -np 1 ../bin/dynamic_<ext> Deep-HP_example1 1000 0.2 100 4 300 1` <br />

* **Example 2:** <br />
:green_circle: *Objective*: Perform hybrid machine learning potential/MM simulation - on the full system.<br />
:large_blue_circle: Simulation parameter: NPT with montecarlo barostat, velocity-verlet integrator and hybrid ANI2X/AMOEBA VdW energy.<br />
<ins>Command</ins>: :arrow_right: `mpirun -np 1 ../bin/dynamic_<ext> Deep-HP_example2 1000 0.2 100 4 300 1`<br />

* **Example 3:**<br />
:green_circle: *Objective*: Perform DeePMD machine learning potential simulation - on the full system.<br />
:large_blue_circle: Simulation parameter: NPT with montecarlo barostat, velocity-verlet integrator and a simple trained DeePMD potential. <br />
<ins>Command</ins>: :arrow_right: `mpirun -np 1 ../bin/dynamic_<ext> Deep-HP_example3 1000 0.2 100 4 300 1`<br />

* **Example 4:**<br />
:green_circle: *Objective*: Perform hybrid machine learning potential/MM simulation - on a ligand of the SAMPL4 challenge.<br />
:large_blue_circle: Simulation parameter: NPT with montecarlo barostat, RESPA integrator with 0.2fs inner time-step/ 2fs outer time-step and ANI2X potential applied only to ligand-ligand interactions (atoms 1 to 24), ligand-water and water-water interactions use AMOEBA.<br />
:yellow_circle: Note: You can try with the `integrator respa1` with the following settings `dshort 0.0002` (0.2fs), `dinter 0.001` (1fs) and use 3fs as outer time-step. It is safer for the dynamics stability to also enable `heavy-hydrogen`.
<ins>Command</ins>: :arrow_right: `mpirun -np 1 ../bin/dynamic_<ext> Deep-HP_example4 1000 2.0 100 4 300 1`<br />

* **Example 5:**<br />
:green_circle: *Objective*: Perform hybrid machine learning potential/MM simulation - on a host-guest complex of the SAMPL4 challenge.<br />
:large_blue_circle: Simulation parameter: NPT with montecarlo barostat, RESPA integrator with 0.2fs inner time-step/ 2fs outer time-step and ANI2X potential applied only to ligand-ligand interactions (atoms 1 to 24), the rest of the interactions use AMOEBA.<br />
:yellow_circle: Note: You can try with the `integrator respa1` with the following setting `dshort 0.0002` (0.2fs), `dinter 0.001` (1fs) and use 3fs as outer time-step. It is safer for the dynamics stability to also enable `heavy-hydrogen`.
<ins>Command</ins>: :arrow_right: `mpirun -np 1 ../bin/dynamic_<ext> Deep-HP_example5 1000 2.0 100 4 300 1`<br />

* **Example 6:**<br />
:green_circle: *Objective:* Perform machine learning potential simulation - on the full system (100000 atoms) with 2 GPUs.<br />
:large_blue_circle: Simulation parameter: NPT with montecarlo barostat, velocity-verlet integrator and ANI2X potential.<br />
<ins>Command</ins>: :arrow_right: `mpirun -np 2 ../bin/dynamic_<ext> Deep-HP_example5 1000 0.2 100 4 300 1`<br />

# Contact
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.<br />

If you want to add your favorite machine learning potential code inside Deep-HP, Contact us:
* jean-philip.piquemal@sorbonne-universite.fr


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
```
