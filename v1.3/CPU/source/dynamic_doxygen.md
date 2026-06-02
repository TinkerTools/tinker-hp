# More information about the dynamic program 

The `dynamic` executable allows to run molecular dynamics simulation. As
for Tinker-8, the command line used to run the MD should give first
the number of MD steps to make, then the size of each time step (in
femtoseconds), then the time between each writing of geometry (in
picoseconds), then the statistical ensemble to sample : 1 is NVE, 2 is
NVT, 4 is NPT. For NVT and NPT simulations, this number should be
followed by the temperature (in Kelvin) of the simulation, and for NPT
simulation by the pressure (in Atmosphere) of the simulation.

For example the command lines:

- `mpirun -np 16 ../bin/dynamic dhfr2 1000 1 1 1`
  will give you as an output 1000 MD steps in NVE for the dhfr2 system,
  with a 1 *fs* time step and a 1 *ps* frequency output.

- `mpirun -np 16 ./dynamic dhfr2 1000 1 1 2 300`
  will give you as an output 1000 MD steps in NVT at 300K for the dhfr2
  system, with a 1 *fs* time step and a 1 *ps* frequency output.

- `mpirun -np 16 ./dynamic dhfr2 1000 1 1 4 300 1`
  will give you as an output 1000 MD steps in NPT at 300K and 1atm for
  the dhfr2 system, with a 1 *fs* time step and a 1 *ps* frequency
  output.

Various keywords to be included in the simulation keyfile are specific to the 'dynamic' program such as the description of integrators, thermostats, barostats  or free-energies. Here is the list of these keywords with a description of their impact on the simulation.

**Integrators**:
Similar to Tinker-8, default integrator is 'Beeman'. Various integrators are available including some multiple-timestep ones.
  - BEEMAN: beeman integrator (variation of velocity verlet integrator)
  - VERLET: velocity verlet integrator
  - BAOAB: Langevin dynamics with BAOAB discretization, already includes the Langevin theromstat
  - RESPA: Respa-like two-level multiple timestep integrator with a velocity verlet inner loop. The inner loop evaluates the bonded potential energy terms and the outer one the non-bonded ones.
  - BAOABRESPA: similar as RESPA with a BAOAB inner loop.
  - RESPA1: Respa1-like three-level multiple timestep integrator with a velocity verlet inner loop. The inner loop evaluates the bonded potential energy terms, the intermediate one the short-range non-bonded ones and the outer one the long-range non-bonded ones.
  - BAOABRESPA1: same as RESPA1 with a BAOAB inner loop.

Multiple timestep integrators can have their shorter (and intermediate for three-level integrators) timesteps controled with the following keywords:
  - DSHORT: shorter timestep in picoseconds, default 0.00025 (0.25 fs)
  - DINTER: intermediate timestep in picoseconds, default 0.002 (2 fs)

All of these integrators are associated with a stability limit (maximum timestep) that can be increased by distributing the mass of heavy atoms on the hydrogens their are bound to, without compromising accuracy of computed observables, this is done with the keyword:
  - HEAVY-HYDROGEN

**Thermostats**:
Aside from Langevin dynamics, temperature controlled can be made through the 'thermostat' keyword with the following arguments:
  - BUSSI: Bussi thermostat, default.
  - BERENDSEN: Berendsen thermostat.
  - TAU-TEMPERATURE: characteristic time of thermostat in 1/ps, default  0.2
  - ANDERSEN: Andersen thermostat.
Regarding Langevin dynamics, the friction can be controllod with the keyword:
  - FRICTION: friction of Langevin dynamics in 1/ps, default 1.

**Barostats**:
Pressure can be controlled can be made through the 'barostat' keyword with the following arguments:
  - BERENDSEN: Berendsen barostat, default.
  - TAU-PRESSURE: characteristic time of barostat in ps, default 2.
  - MONTECARLO: Monte-carlo barostat
  - LANGEVIN: Langevin barostat, extended dynamics of the volume, requires a mass and a friction than can be controlled with the following keywords:
   - FRICTIONPISTON: friction of Langevin barostat in 1/ps, default 20.
   - MASSPISTON: mass of the piston for Langevin barostat, default 1e5
Semi-isotropic pressure can be imposed by freezing some axis of the simulation box with the following keyword:
  - FREEZE-A-AXIS
  - FREEZE-B-AXIS
  - FREEZE-C-AXIS


**CONSTRAINED DYNAMIC**:
Constrained MD can be run thanks to the RATTLE algorithm in combination with the non-Langevin integrators, the constrains can can be controlled with the following keywords:
  - RATTLE-EPS: convergence criterion of RATTLE, default 1e-6
  - RATTLE WATER: fix all distances and angles of water molecules to their equilibrium
  - RATTLE BONDS: fix all bonds to their equilibrium
  - RATTLE ANGLES: fix all angles to their equilibrium
  - RATTLE DIATOMIC: takes two connected atom indexes in argument, fix the corresponding distance to its equilibrium
  - RATTLE TRIATOMIC: takes three connected atom indexes in argument, fix the corresponding bonds ang angles to their equilibrium
  - RATTLE-DISTANCE: takes two atom indexes and a distance in argument, fix the corresponding distance to this value

**NEIGHBOR-LISTS**:
The non-bonded potential energy terms evaluation relies on neighbor lists built with a cell-list algorithm every 'x' timestep with a skin of a few Angstroms. This be modified with the keywords:
  - NLUPDATE: number of timesteps between each neighbor list update, default 20 for a 2fs integration.
  - LIST-BUFFER: skin (margin) taken to build neighbor lists, default 2 Angstroms.

- **parallelism**:
  - dd-cutoff: additional cutoff for communication between domain for 3d spatial decompostion, can be larger than the largest non-bonded cutoff
As the reciprocal space interactions are known to have a less efficient parallel scaling
(because of FFTs) it is possible to specify a lower number of <span class="smallcaps">mpi</span> processes that will be dedicated to
these computations for both the computation of electrostatic and
polarization interactions. This can be done by using the keyword
**pme-procs** `x`, corresponding to `x`
<span class="smallcaps">mpi</span> processes dedicated to reciprocal
space computations.

To find the ideal `x` value of **pme-procs**, a good starting point,
when a large number of cores is used, is usually to dedicate about
\f$\frac{1}{4}\f$ of the total cores to reciprocal space computations. But
this parameter depends greatly on the machine used and on the setup of
your simulation so it should be adjusted manually by comparing the
timings obtained with different values for **pme-procs**. During a
dynamic, detailed timings are written when the **verbose** keyword is in
the \*.key file. In the future, the parameter **pme-procs** will be
adapted heuristically by the program as it is done in popular MD
packages.

**OTHERS**:
  - REMOVE-INERTIA: number of timesteps between removal of center of mass translation and rotational kinetic energy, default 0 (no such removal)
  - COMPRESS: isothermal compressibility of medium in 1/Atm
  - VOLUME-MOVE: maximum volume move for monte-carlo barostat
  - VIRNUM: compute numerical virial (with finite differences)
  - IR_SPECTRA: compute infrared spectra
  - IR_LINEAR_DIPOLE: include only permanent dipole for IR spectra computation (no induced dipoles)
  - IR_DECONVOLUTION: use deconvolution procedure to compute IR spectra with Langevin dynamics
  - NITER_DECONV: number of iteration for the deconvolution procedure
