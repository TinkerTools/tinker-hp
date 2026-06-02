# Free energy calculations {#free-energy}

Tinker-HP supports various free-energy techniques, mainly through the support of the external COLVARS and PLUMED libraries, if these are linked as described in the documentation concerning building of the package.

Regarding PLUMED, the line:
  - plumed input.dat output.dat, has to be added to the keyfile, with input.dat the PLUMED input file and output.dat the output file.

Regarding COLVARS, a COLVARS configuration file with the same prefix as the keyfile and the extension .colvars has to be added. All the necessary information about coupling of Tinker-HP and Colvars can be found in the dedicated [manual](http://colvars.github.io/master/colvars-refman-tinkerhp.pdf).

Tinker-HP also handles alchemical free energy simulations where a list of 'alchemical' atoms has to be provided to the keyfile with the line:
- LIGAND 

'Fixed lambda' trajectories can be produced by providing information about atoms included in the alchemical region and providing alchemical lambda values for van der Waals and electrostatics (with polarization for AMOEBA-like force fields). The two lambda values can be modified with the keywords:
- ele-lambda
- vdw-lambda

One can also control a global lambda parameters of which the two specific lambdas above are a function. By default vdw-lambda grows linearly wrt lambda from 0 to 0.5 and similarly ele-lambda grows linearly from 0.5 to 1. 

The value of lambda where vdw-lambda is 1 can be controlled by the keyword:
- BOUND-VDW-LAMBDA , default 0.5

The value of lambda where ele-lambda is 0 can be controlled by the keyword:
- BOUND-ELE-LAMBDA , default 0.5

Furthermore, softcore potential for van der Waals interactions are implemented and the associated parameters can be modified with the following keywords:

For Lennard-Jones interactions:
- VDW-SC-ALPHA , default 0.7
- VDW-SC-EXP, default  5
- VDW-SC-K , default 6
- VDW-SC-S , default 2
- VDW-SC-T , default 1

For Halgren buffered 14-7 potential:
- VDW-SC-ALPHA , default 0.7
- VDW-SC-EXP, default 5

By default, scaled van der Waals interactions are "decoupled" meaning that intramolecular interactions are kept intact, but one can also progressively scale down intramolecular van der Waals interactions, by using the keyword:
- VDW-ANNIHILATE

More advanded and efficient alchemical free-energy simulations can be run by leveraging the Lambda-ABF methodology as described in the [dedicated page](md_lambda-abf_doxygen.html).

Note that both PLUMED and COLVARS based simulations are compatible with all multi-timestep integrators, are these external libraries are called at the outer timestep in this case.

