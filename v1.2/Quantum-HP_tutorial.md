[![DOI:10.1021/acs.jctc.2c01233](https://zenodo.org/badge/DOI/10.1021/acs.jctc.2c01233.svg)](https://doi.org/10.1021/acs.jctc.2c01233) 



# Quantum-HP: Fast Nuclear Quantum Effects for Large-Scale MD Simulations

Quantum-HP is an extension of the Tinker-HP molecular dynamics (MD) package that allows for the explicit inclusion of nuclear quantum effects (NQEs) in large-scale simulations. Quantum-HP leverages both Ring-Polymer Molecular Dynamics (RPMD) and adaptive Quantum Thermal Bath (adQTB) simulation methods to efficiently incorporate NQEs in MD simulations. This extension opens up new possibilities for studying the quantum nature of important interactions in biological matter and other condensed-phase systems.

## What is Quantum-HP?
Quantum-HP is a multi-CPU and multi-GPU massively parallel platform designed to incorporate nuclear quantum effects in Tinker-HP MD simulations. It utilizes two primary simulation strategies:

1. **Ring-Polymer Molecular Dynamics (RPMD)**: Based on Feynman's path integrals formalism, this approach provides exact structural properties of a quantum system by simulating an effective classical system composed of multiple replicas of the original system. RPMD is computationally expensive but highly accurate.<br/>
Further reading: https://doi.org/10.1146/annurev-physchem-040412-110122

2. **Adaptive Quantum Thermal Bath (adQTB)**: adQTB imposes the quantum distribution of energy on a classical system using a generalized Langevin thermostat. This strategy offers computationally affordable (close to classical MD cost) and reasonably accurate (though approximate) treatment of NQEs.<br/>
Further reading: https://doi.org/10.1021/acs.jctc.8b01164

## Key Features
- Explicit inclusion of nuclear quantum effects (NQEs) in large-scale MD simulations.
- Utilizes both Ring-Polymer Molecular Dynamics (RPMD) and adaptive Quantum Thermal Bath (adQTB) simulation strategies.
- Efficient implementation with parallelization for multi-CPU and multi-GPU systems.
- Compatible with BAOAB-RESPA multi-timestep integrator 
- Multiple ring-polymer contraction schemes for efficient RPMD simulations.
- Compatible with all Tinker-HP's force fields, including the AMOEBA polarizable force field.
- Compatible with Tinker-HP's Deep-HP machine learning potentials module, enabling combined quantum and machine learning simulations.
- Allows for alchemical free energy estimation in conjunction with Tinker-HP.
- Capable of simulating systems with a significant number of atoms, making it suitable for studying large condensed-phase systems.

## Running Quantum-HP
To use Quantum-HP, follow these steps:

1. **Installation**: Quantum-HP is provided by default with Tinker-HP so you just need to install a recent version of Tinker-HP (later than 2023/08)

2. **Setting up the Input**: Prepare the input files for your MD simulation using Tinker-HP's standard input format. Include the necessary parameters to specify the use of Quantum-HP for nuclear quantum effects (see list of keywords and examples below).

3. **Selecting Simulation Strategy**: Choose between RPMD and adQTB based on the level of accuracy required and the computational resources available. RPMD provides exact results but is computationally more expensive than adQTB.

4. **Parallelization**: Quantum-HP is optimized for parallel computing. Make sure to take advantage of multi-CPU and (multi-)GPU systems to accelerate your simulations.

5. **Running the Simulation**:
    - **adQTB:** Running (ad)QTB MD is very similar to running classical MD. Simply run the `dynamic` executable with the usual command line arguments.
    The QTB or adQTB is activated by setting the `THERMOSTAT` keyword to `adQTB` in the .key file. Available integrators for adQTB are `BAOAB`,`BAOABRESPA` and `BAOABRESPA1`.  
    - **RPMD:** path integrals simulations can be performed using the `pimd` executable. Command line arguments for `pimd` are the same as for `dynamic`. Available integrators for path-integrals simulations are `BAOAB` and `BAOABRESPA`.

## Examples
The following examples show the basic usage of the Quantum-HP module. We will start from a RPMD simulation of a water box and check the convergence with increasing number of replicas. Then, we will run a standard QTB simulation and compare the results. We will then analyze the zero-point energy leakage (ZPEL) that is the main pitfall of the standard QTB methodology. Finally, we will perform an adQTB simulation to see how it fixes the ZPEL by adjusting the parameters of the thermal bath.<br/>
We will use the fixed-charge q-SPC/Fw (https://doi.org/10.1063/1.2386157) force field as our water model in order to perform quick simulations. For polarizable AMOEBA simulations, you can use the parameter files `qwater22-QTB.prm` and `qwater22-PI.prm` provided in the `params` directory. These are the quantum-parameterized water models from https://doi.org/10.1021/acs.jpcb.2c04454

- **Example 1: RPMD simulation of a water box**
  
    In this example, we will run a standard RPMD simulation of a small water box with the q-SPC/Fw force field. The input files are provided in the `examples/NQE/liquid_water/RPMD` folder. To run the simulation, execute the following command (with the correct path to the `pimd` executable):
    ```
    mpirun -np 1 /path/to/pimd watersmall 10000 2 1 2 300
    ```
    This command will run a 10000 step simulation at 300K with 32 replicas (if the simulation takes too long, do not hesitate to parallelize by editing the `-np` argument in the previous command). Notice that the estimated temperature is much larger than 300K (~900-1000K). This temperature is computed from the quantum kinetic energy estimator (centroid-virial) and thus includes zero-point energy effects ! We can check that the simulation is indeed running at 300K by checking the average centroid temperature. <br/>
    To check the convergence of the RPMD simulation, you can vary the number of replicas with the `NBEADS` keyword in `watersmall.key` from 1 replica (classical MD) to as many as your machine can handle ! Note that 32 replicas are typically used for simulations containing light atoms (Hydrogen) at 300K. This number should increase at lower temperature and convergence should always be checked when simulating a new system.<br/>
    For RPMD simulations, Tinker-HP saves a trajectory file for each replica (by default) here named `watersmall_beads*.arc`. These files can be used independently to compute configurational observables (for example radial distribution functions) just like classical trajectories. Observables can then be averaged over the replicas.<br/>
    <br/>
    The provided script `examples/NQE/liquid_water/compute_rdf.py` can be used to compute the radial distribution functions (RDFs) of water from the trajectory files (python3 and numpy are required). To compute the RDFs, run the following commands from the `examples/NQE/liquid_water/RPMD` directory:
    ```
    cat watersmall_beads*.arc > watersmall.arc
    python3 ../compute_rdf.py watersmall.arc
    ```
    which will generate the file `gr.dat` containing the O-O, O-H and H-H RDFs.
    The files `gr_32beads.ref` and `gr_classical.ref` provide reference results for the radial distribution function of water computed from the 32-beads RPMD simulation and a classical MD simulation, respectively.


- **Example 2: standard QTB simulation**
    
    We will now use the Quantum Thermal Bath to accelerate the simulation. The input files are provided in the `examples/NQE/liquid_water/QTB` folder. To run the simulation, execute the following command (with the correct path to the `dynamic` executable):
    ```
    mpirun -np 1 /path/to/dynamic watersmall 100000 2 1 2 300
    ```
    This should run much faster than the 32-beads RPMD simulation. The kinetic energy should be very similar to the path integrals one and potential energy should be slightly overestimated. One can then compute radial distribution functions (RDF) from the trajectory file `watersmall.arc` and compare them to the reference results in `gr_QTB.ref` and to the results from the previous example. If you visualize the RDFs, you will see that the peaks are broader in QTB than in RPMD and classical MD. This is a typical indication of zero-point energy leakage effects. <br/>
    We can analyze more in details this effect using the `QTB_spectra_*.out` file that has been generated by the dynamics (or the provided reference `QTB_spectra_*.out.ref`). Columns 2 and 3 of this file contain the velocity-velocity and velocity-random force correlation spectra estimated during the dynamics that are used to quantify the ZPEL (see https://doi.org/10.1021/acs.jctc.8b01164, section III). Their difference (column 4) is a measure of where (in the frequency domain) is energy missing or overflowing. In this example, one can see a dip at ~3600 cm^{-1} indicating that energy is missing in the O-H stretching mode. Indeed, this mode contains a large amount of zero-point energy that "leaks" into the lower-frequency modes during the QTB simulation. <br/>
    In the next example, we will see that the adaptive QTB (adQTB) fixes the ZPEL.

- **Example 3: adQTB simulation**
    We will now run an adQTB simulation of the same water box. The input files are provided in the `examples/NQE/liquid_water/adQTB` folder. To run the simulation, execute the following command (with the correct path to the `dynamic` executable):
    ```
    mpirun -np 1 /path/to/dynamic watersmall 100000 2 1 2 300
    ```
    During the simulation, Tinker-HP should periodically print "Adapting gammar for adQTB" (with a period given by the `TSEG` keyword), which indicates that the parameters of the thermal bath are iteratively adapted to fix the ZPEL (i.e. minimizing the difference between the two spectra mentionned in the previous section). The file `gamma_restart.out` contains the values of the current parameters for each atom type, as a function of frequency. The base value should be 20 ps^{-1} and you should see a decrease at low frequencies and an increase at higher frequencies. <br/>
    When the parameters converge and start fluctuating around a fixed value, this means that the adaptation procedure compensated the ZPEL (reference parameters are provided in `gamma_restart.out.ref`). You can check that the difference between the two spectra is now very small (column 4 of `QTB_spectra_*.out`) compared to the previous QTB calculation. <br/>
    You can now compute the RDF from the trajectory file `watersmall.arc` and compare it to the reference results in `gr_adQTB.ref` and to the results in the previous examples. You should see that the peaks are now sharper than in the QTB simulation and very close to the RPMD results. <br/>
    You can experiment with the adaptation procedure (with the keyword `ADQTB_OPTIMIZER`) and its associated parameters in `watersmall.key` to modify the speed at which the parameters are optimized (and the associated level of noise induced in `gamma_restart.out`).

- **Example 4: adQTB simulation of solvated Benzene**
    In adQTB a set of adjustable parameters (`gamma_restart.out`) is associated to each atom type. The corresponding spectra `QTB_spectra_*.out` are then averaged over all equivalent degrees of freedom. Spectra for types that have lots of representative atoms are then estimated much quicker than for types with a few representative atoms. This situation is typical of a small molecule solvated in a box of water. In this case, one can adjust the parameters for the solvent much more agressively than that of the solute which require a slower tuning.<br/>
    Tinker-HP thus provides keywords for specifying the adaptation scheme for specific atom types, as can be seen in `examples/NQE/solvated_benzene/adQTB/benzenebox.key`. 

# List of Keywords
This section provides a list and a brief explanation of the keywords specific to Quantum-HP.

- **(ad)QTB:**
   - `FRICTION` (real): the friction coefficient of the Langevin thermostat   (default: 1., unit: ps^{-1}). Typical values for (ad)QTB are 10-20 ps^{-1}   (default value is set for classical MD).
   - `OMEGACUT` (real): the cutoff frequency for computing spectra and beyond which   the QTB power spectrum is zero (default: 15000, unit:cm^{-1})
   - `TSEG` (real): the duration of a trajectory segment on which to compute spectra   (default: 1., unit: ps)
   - `QTB_VERBOSE` (logical): print to standard output when refreshing the QTB noise and updating the adQTB parameters
   - `QTB_BATCH_SIZE` (int): number of atoms on which to generate the colored noise   simultaneously (only for the GPU version). Allows to reduce memory requirements   for large systems but might be slower when using small batch sizes. Defaults to   all the local atoms (i.e. no batches).
   - `CORR_FACT_QTB` (real): set the value of the kinetic energy correction. If not present, the correction is automatically estimated along the dynamics.
   - `NO_CORR_POT` (logical): disables the potential energy correction.
   - `REGISTER_SPECTRA` (logical): save spectra for the zero-point energy leakage   diagnostic (Cvv, Cvf, Delta_FDT). Automatically ON when performing adQTB   simulation.
   - `SKIPSEG` (int): skip the first `SKIPSEG` segments before starting the   adaptation (default: 3)
   - `STARTSAVESPEC` (int): segment from which spectra start to be averaged. Allows   to save spectra only after an equilibration/adaptation period (default: 25)
   - `GAMMA_HISTORY` (logical): saves the evolution of the gamma coefficients at   each segment (warning: may produce large files for long dynamics).
   - `NOQTB` (logical): use the QTB noise generation with a white noise kernel to   perform classical MD. Mostly for debug or testing purposes.
   - `ADQTB_OPTIMIZER` (string): the adaptation method to use (possible values:   `SIMPLE`,`RATIO`)
   - `A_GAMMA` (real): coefficient that controls the adaptation speed in the   `SIMPLE` adaptation mode (default: 0.5). 
   - `ADQTB_TAU_AVG` (real): coefficient that controls the adaptation speed in the   `RATIO` adaptation mode (default:10*`TSEG`, unit: ps)
   - `ADQTB_TAU_ADAPT` (real): 
       - for the `SIMPLE` method: typical time of window averaging of Delta_FDT.   Allows in principle to improve the adaptation, similarly to stochastic   gradient descent with momentum  (default:0, unit: ps)  
       - for the `RATIO` method: typical time of window averaging of the gamma   coefficients. Smoothens gamma but slows down the adaptation (default:0, unit:   ps) 
   - `ADQTB_SMOOTH` (real): smoothens the adaptation by window-averaging the spectra   with an exponentially decaying window of size of `ADQTB_SMOOTH` cm^{-1} (default:   0, unit: cm^{-1})
   - `ADQTB_MODIFIER` (int,string,real,real): specify adaptation method (`RATIO` or   `SIMPLE`), associated parameter (`ADQTB_TAU_AVG` or `A_GAMMA`) and `ADQTB_SMOOTH` for a specific type. The atom type is specified by the first integer argument. The second argument is the adaptation method. The third and fourth arguments are the associated parameter (`ADQTB_TAU_AVG` or `A_GAMMA`) and `ADQTB_SMOOTH` (the default parameter is obtained by setting to a negative value).

- **RPMD:**
    - `NBEADS` (int): number of beads to use in the simulation (default: 1)
    - `NBEADS_CTR` (int): number of beads of the contracted ring polymer (default:    0=not used)
    - `CENTROID_LONGRANGE` (logical): compute long-range forces on the centroid only. 
    - `POLAR_CENTROID` (logical): compute the full polarization interactions on the   centroid only. Only used if `CENTROID_LONGRANGE` is activated. If not present   while `CENTROID_LONGRANGE` is present, long-range polarization is calculated on  the centroid and short-range polarization on the beads (or contracted beads if   `NBEADS_CTR`>1).
    - `LAMBDA_TRPMD`: lambda value for the TRPMD friction coefficient. See https://doi.org/10.1063/1.4883861 section II-C. Note that there is a factor 2 between the definition in Quantum-HP and the reference paper: `LAMBDA_TRPMD`=1 corresponds to lambda=0.5 in the article. The default behavior is the optimally thermostatted TRPMD (`LAMBDA_TRPMD`=1), standard RPMD can be obtained with `LAMBDA_TRPMD`=0.
    - `CAY_CORR` (logical): activates the Cayleigh correction to use the BCOCB    integrator of https://doi.org/10.1063/1.5134810. Allows in principle to use   larger time steps when large number of beads are required.
    - `ADD_BUFFER_PI` (logical): add a buffer to the neighbor lists that accounts for the delocalization of the ring polymer. Might be useful at low temperatures when delocalization is larger than default buffers.
    - `PI_START_ISOBARIC`: time at which to activate the barostat. Allows to    thermalize the ring polymer in NVT (to have unbiased pressure estimator) before   performing the NPT simulation (default: 0, unit: ps)
    - `PITIMER`: display timers specific to the path-integral simulation.
    - `SAVE_BEADS` (string): controls which beads trajectories are saved. Possible    values: `ALL` (default), `RANDOM` and `ONE`.
         - `ALL`: save all the beads trajectories
         - `RANDOM`: save the trajectory of a random bead (sampled at each frame)
         - `ONE`: save the trajectory of the first bead only

## (ad)QTB-specific Output files:
- `QTB_spectra_*.out`: spectra for the zero-point energy leakage diagnostic(mCvv, Cvf/gamma, Delta_FDT, gamma) for each atom type.
- `corr_pot.dat`: potential energy correction (used at the beginning of thedynamics if present or generated automatically if not present).
- `corr_fact_qtb_restart.out`: kinetic energy correction (used at the beginningof the dynamics if present or generated automatically if `CORR_FACT_QTB` is not set).
- `gamma_restart.out`: gamma coefficients for all the atom types (used at thebeginning of the dynamics to restart from previous simulation).

# Extra Features
### Path-Integral Bennett Acceptance Ratio (PI-BAR)

PI-BAR free energy estimation can be performed using the `pibar` executable. As for the standard `bar` executable, two modes are available from the command line: 
 1. the first mode creates the .bar file by computing the potential in the two thermodynamical states for each frame in both trajectories. The potential is averaged over the beads to be able to compute the PI-BAR estimator. If the .key file defines a contracted ring polymer, the average potential is also computed on the contracted ring-polymer and saved in a .bar.ctr file.
 2. the second mode computes the free energy from the .bar file. If a .bar.ctr file is present, the user is prompted to enter if they want to use it or do the estimation from the full ring polymer.

### On-the-fly calculation of infrared spectra

Infrared spectra can be computed from the total dipole moment (defined by the multipoles or charges). To activate the calculation of IR spectra, add the `IR_SPECTRA` keyword to the .key file.  
The IR spectra are computed in segments as defined for QTB dynamics and are thus influenced by the `OMEGACUT` and `TSEG` keywords.
Frequencies and IR spectra are saved in cm^{-1}. 
In the case of multipoles, the dipole moment can be approximated using charges only with the keyword `IR_LINEAR_DIPOLE`.  
The spectra can be deconvoluted to remove the influence of the friction coefficient in a Langevin simulation with the keyword `IR_DECONVOLUTION` (https://doi.org/10.1063/1.4990536).

# Contact
For any bug report, technical question or feature request, contact us at:
* thomas.ple@sorbonne-universite.fr
* louis.lagardere@sorbonne-universite.fr
* jean-philip.piquemal@sorbonne-universite.fr


# Please Cite
```tex
@article{https://doi.org/10.1021/acs.jctc.2c01233,
  title={Routine molecular dynamics simulations including nuclear quantum effects: from force fields to machine learning potentials},
  author={Pl{\'e}, Thomas and Mauger, Nastasia and Adjoua, Olivier and Inizan, Th{\'e}o Jaffrelot and Lagard{\`e}re, Louis and Huppert, Simon and Piquemal, Jean-Philip},
  journal={Journal of Chemical Theory and Computation},
  volume={19},
  number={5},
  pages={1432--1445},
  year={2023},
  publisher={ACS Publications}
}
```
