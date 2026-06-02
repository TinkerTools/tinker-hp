# List of Potential Energy functions {#list_potential}

Various bonded and non-bonded energy terms are implemented with analytical gradients in Tinker-HP. Parameter files associated to well established force fields such as charmm22, ambers99 or amoeba can be found in the parameters directory of the repository.

Information about coupling of the GPU version of Tinker-HP with machine-learning libraries in order to run simulations leveraging machine learning potentials can be found [here](md_Deep-HP.html).

- **Bond Stretching**:
  - description: 
    - standard bond-stretching wrt to an equilibrium distance, function  of distances between bonded atoms. 
    - can be HARMONIC as Taylor expansion of Morse potential through the 4th power of the bond length deviation
    - can be MORSE or Morse potential expanded to the 4th order (MORSE4)
  - keywords: 
    - BONDTERM only to restrict the potential to bond-stretching, none to remove all bond-stretching interactions
    - BONDTYP: type of bond-stretching interactions, can be **HARMONIC**, **MORSE** or **MORSE4**
    - BONDUNIT: conversion factor of bond stretch energy to kcal/mole
    - BOND-CUBIC: cubic coefficient in bond stretch potential
    - BOND-QUARTIC: quartic coefficient in bond stretch potential
  - files: kbond.f90, ebond.f90, ebond1.f90, ebond3.f90

- **Angle Bending**:
  - description:
    - standard angle-bending wrt to an equilibrium angle, function of angles between bonded atoms.
    - can be in-plane,  or not. For in-plane, can be harmonic, linear, fourier or gaussian
  - keywords:
    - ANGLETERM only to restrict the potential to angle-bending, none to remove all angle-bending interactions
    - ANGLEUNIT: conversion factor of angle bending energy to kcal/mole
    - ANGLE-CUBIC: cubic coefficient in angle bending potential
    - ANGLE-QUARTIC: quartic coefficient in angle bending potential
    - ANGLE-PENTIC: pentic coefficient in angle bending potential
    - ANGLE-SEXTIC: sextic coefficient in angle bending potential
  - files: kangle.f90, eangle.f90, eangle1.f90, eangle3.f90

- **Stretch-Bending**:
  - description:
    - given 3 bonded atoms, coupling between the harmonic energy of the distances wrt to their equilibrium and the harmonic energy of the angle wrt to its equilibrium
  - keywords:
    - STRBNDTERM only to restrict the potential to strech-bend coupling, none to remove all the associated interactions
    - STRBNDUNIT: conversion factor of stretch bending energy to kcal/mole
  - files: kstrbnd.f90, estrbnd.f90, estrbnd1.f90, estrbnd3.f90


- **Urey-Bradley**:
  - description:
    - standard Urey-Bradley term, potential on distances between outer atoms of bonded angles
    - can be repulsive (exponential), quartic or harmonic with taylor expansion up to quartic
  - keywords
    - UREYTERM only to restrict the potential to urey-bradley, none to remove all urey-bradley interactions
    - UREYUNIT: conversion factor of Urey-Bradley energy to kcal/mole
    - UREY-CUBIC: cubic coefficient in Urey-Bradley potential
    - UREY-QUARTIC: quartic coefficient in Urey-Bradley potential
  - files: kurey.f90, eurey.f90, eurey1.f90, eurey3.f90

- **Angle-Angle**:
  - description:
    - given two angles, coupling between harmonic energy of angles wrt to their equilibrium 
  - keywords:
    - ANGANGTERM: only to restrict the potential to angle-angle coupling, none to remove all the associated interactions
    - ANGANGUNIT: conversion factor of angle-angle energy to kcal/mole 
  - files: kangang.f90, eangang.f90, eangang3.f90, eangang1.f90

- **Out-of-plane bending**: 
  - description: 
    - out-of-plane bend at trigonal centers via Wilson-Decius-Cross or Allinger angles
  - keywords:
    - OPBENDTERM: only to restrict the potential to angle-angle coupling, none to remove all the associated interactions
    - OPBENDTYPE: type of out-of-plane bending, can be **W-D-C** or **ALLINGER**
    - OPBENDUNIT: conversion of out-of-plane bending energy to kcal/mole
    - OPBEND-CUBIC: cubic coefficient in out-of-plane bending potential
    - OPBEND-QUARTIC: quartic coefficient in out-of-plane bending potential
    - OPBEND-PENTIC: pentic coefficient in out-of-plane bending potential
    - OPBEND-SEXTIC: sextic coefficient in out-of-plane bending potential

  - files: kopbend.f90, eopbend.f90, eopbend1.f90, eopbend3.f90

- **Out-of-plane distance**:
  - description: potential at trigonal center via the central atom height
  - keywords: 
    - OPDISTTERM only to restrict the potential to out-of-plane distance coupling, none to remove all the associated interactions
    - OPDISTUNIT: conversion of out-of-plane distance potential to kcal/mole
    - OPDIST-CUBIC: cubic coefficient in out-of-plane-distance potential
    - OPDIST-QUARTIC: quartic coefficient in out-of-plane-distance potential
    - OPDIST-PENTIC: pentic coefficient in out-of-plane-distance potential
    - OPDIST-SEXTIC: sextic coefficient in out-of-plane-distance potential

- **Improper Dihedral**:
  - description: standard harmonic potential on improper dihedrals
  - keywords: 
    - IMPROPTERM: only to restrict the potential to improper dihedral term, none to remove all the associated interactions
    - IMPROPUNIT: conversion factor of improper dihedrals energy to kcal/mole
  - files: kimprop.f90, eimprop.f90, eimprop1.f90, eimprop3.f90

- **Improper Torsion**:
  - description: periodic improper torsion potential
  - keywords:
    -IMPTORTERM: only to restrict the potential to improper torsion term, none to remove all the associated interactions

    -IMPTORUNIT: conversion factor of improper torsions energy to kcal/mole
  - files: kimptor.f90, eimptor.f90, eimptor1.f90, eimptor3.f90

- **Torsions**:
  - standard periodic torsion potential
  - keywords:
    - TORSIONTERM: only to restrict the potential to torsion term, none to remove all the associated interactions
    - TORSIONUNIT: conversion factor of torsion energy to kcal/mole
  - files: torsions.f90, ktors.f90, etors.f90, etors1.f90, etors3.f90

- **Pi-Torsions**: 
  - description: Pi-Orbital torsion potential
  - keywords: 
    - PITORSTERM: only to restrict the potential to Pi-torsion term, none to remove all the associated interactions
    - PITORSUNIT: conversion factor of Pi-torsion energy to kcal/mole
  - files: kpitors.f90, epitors.f90, epitors1.f90, epitors3.f90

- **Strech-Torsion**: 
  - description: given 4 atoms involved in a torsion, Stretch-Torsion coupling between the first distance and the torsion
  - keywords:
    - STRTORTERM: only to restrict the potential to Stretch-torsion term, none to remove all the associated interactions
    - STRTORUNIT: conversion factor of stretch-torsion energy to kcal/mole
  - files: kstrtor.f90, estrtor.f90, estrtor1.f90, estrtor3.f90

- **Angle-Torsion**:
  - description: given 4 atoms involved in a torsion, Angle-Torsion coupling between the first angle and the torsion
  - keywords:
    -ANGTORTERM: only to restrict the potential to angle-torsion term, none to remove all the associated interactions
    -ANGTORUNIT: conversion factor of angle-torsion energy to kcal/mole
  - files: kangtor.f90, eangtor.f90, eangtor1.f90, eangtor3.f90

- **Torsion-Torsion**:
  - description: coupling between two torsions
  - keywords:
    - TORTORTERM: only to restrict the potential to torsion-torsion term, none to remove all the associated interactions
    - TORTORUNIT: conversion factor of torsion-torsion energy to kcal/mole
  - files: ktortor.f90, etortor.f90, etortor1.f90, etortor3.f90

- **Lennard-Jones**:
  - description: standard 6-12 Lennard-Jones van der Waals energy
  - keywords:
    - VDWTERM: only to restrict the potential to van der Waals term, none to remove all the associated interactions
    - VDWTYP: type of van der Waals interactions, can be **LENNARD-JONES** or **BUFFERED-14-7** (Halgren)
    - RADIUSTYPE: type of radius for LJ interactions, can be **SIGMA** or **R-MIN** (default)
    - RADIUSSIZE: size of radius for LJ interactions, can be **DIAMTER** or **RADIUS** (default)
    - RADIUSRULE: combination rule for radius, can be **ARITHMETIC** (default), **GEOMETRIC** or **CUBIC-MEAN**
    - EPSILONRULE: combination rule for epsilon, can be **ARITHMETIC**, **GEOMETRIC** (default), **HARMONIC**, **HHG** or **WH**
    - VDW-12-SCALE: scaling factor for van der Waals interactions between 1-2 connected atoms (default 0)
    - VDW-13-SCALE: scaling factor for van der Waals interactions between 1-3 connected atoms (default 0)
    - VDW-14-SCALE: scaling factor for van der Waals interactions between 1-4 connected atoms (default 1)
    - VDW-15-SCALE: scaling factor for van der Waals interactions between 1-5 connected atoms (default 1)
    - VDW-CORRECTION: apply long range van der Waals correction
    - VDW-TAPER: proportion of the cutoff after which to apply tapering for van der Waals interactions (default 0.9)
    - VDW-CUTOFF: cutoff for van der Waals interactions (default 9 Angstroms)
    - VDWSHORT-CUTOFF: cutoff for short range  van der Waals interactions (default 7 Angstroms)
 - files: kvdw.f90, elj.f90, elj1.f90, elj3.f90

- **Halgren**:
  - description: buffered 7-14 Halgren van der Waals energy
  - keywords:
    - VDWTERM: only to restrict the potential to van der Waals term, none to remove all the associated interactions
    - GAMMA-HALGREN: gamma parameter for Halgren potential (default 0.12, AMOEBA force field value)
    - DELTA-HALGREN: delta parameter for Halgren potential (defulat 0.07, AMOEBA force field value)
    - VDW-12-SCALE: same as L-J
    - VDW-13-SCALE: same as L-J
    - VDW-14-SCALE: same as L-J
    - VDW-15-SCALE: same as L-J
    - VDW-CORRECTION: same as L-J
    - VDW-TAPER:  same as L-J
    - VDW-CUTOFF: same as L-J
    - VDWSHORT-CUTOFF: same as L-J
 - files: kvdw.f90, ehal.f90, ehal1.f90, ehal3.f90

- **Repulsion**: 
  - description: 'classical' Pauli repulsion energy based on multipolar expansion
  - keywords: 
    - REPULSIONTERM: only to restrict the potential to repulsion term, none to remove all the associated interactions
    - REP-12-SCALE: scaling factor for repulsion interactions between 1-2 connected atoms
    - REP-13-SCALE: scaling factor for repulsion interactions between 1-3 connected atoms
    - REP-14-SCALE: scaling factor for repulsion interactions between 1-4 connected atoms
    - REP-15-SCALE: scaling factor for repulsion interactions between 1-5 connected atoms
    - REPULS-TAPER: proportion of the cutoff after which to apply tapering for repulsion interactions (default 0.9)
    - REPULS-CUTOFF: cutoff for repulsion interactions (default 9 Angstroms)
    - REPULSSHORT-CUTOFF: cutoff for short range repulsion interactions (default 7 Angstroms)
  - files: krepel.f90, erepel.f90, erepel1.f90, erepel3.f90

- **Dispersion**:
  - description: damped dispersion energy
  - keywords:
    - DISPTERM: only to restrict the potential to dispersion term, none to remove all the associated interactions
    - DISP-TAPER: proportion of the cutoff after which to apply tapering for dispersion interactions (default 0.9)
    - DISP-CUTOFF: cutoff for dispersion interactions (default 9 Angstroms)
    - DISPSHORT-CUTOFF: cutoff for short-range dispersion interactions (default 7 Angstroms)
  - files: kdisp.f90, edisp.f90, edisp1.f90, edisp3.f90

- **Charge**:
  - description: charge-charge permanent electrostatics, only SPME available
  - keywords:
    - CHARGETERM: only to restrict the potential to charge term, none to remove all the associated interactions
    - ELECTRIC: energy factor in kcal/mol for current force field
    - DIELECTRIC: dielectric constant for electrostatic interactions
    - CHG-12-SCALE: scaling factor for charge interactions between 1-2 connected atoms
    - CHG-13-SCALE: scaling factor for charge interactions between 1-3 connected atoms
    - CHG-14-SCALE: scaling factor for charge interactions between 1-4 connected atoms
    - CHG-15-SCALE: scaling factor for charge interactions between 1-5 connected atoms
    - CHG-CUTOFF: cutoff for charge interactions (default 7 Angstroms with PME)
    - CHGSHORT-CUTOFF: cutoff for short-range charge interactions (default 5 Angstroms)

- **Multipole**:
  - description: multipole-multipole (up to quadrupoles) permanent electrostatics, only SPME available
  - keywords:
    - MULTIPOLETERM: only to restrict the potential to multipole term, none to remove all the associated interactions
    - MPOLETERM: same as 'MULTIPOLETERM'
    - ELECTRIC: energy factor in kcal/mol for current force field
    - DIELECTRIC: dielectric constant for electrostatic interactions
    - PENETRATION: type of charge-penetration, can be 'NONE', 'GORDON1' or 'GORDON2'
    - MPOLE-12-SCALE: scaling factor for multipole interactions between 1-2 connected atoms
    - MPOLE-13-SCALE: scaling factor for multipole interactions between 1-3 connected atoms
    - MPOLE-14-SCALE: scaling factor for multipole interactions between 1-4 connected atoms
    - MPOLE-15-SCALE: scaling factor for multipole interactions between 1-5 connected atoms
    - MPOLE-CUTOFF: cutoff for mulipole interactions (with PME, default 7 Angstroms)
    - MPOLESHORT-CUTOFF: cutoff for short-range multipole interactions (default 5 Angstroms)

    
- **Polarization**:
  - description: polarization through induced dipoles
  - keywords:
    - POLARIZETERM: only to restrict the potential to polarization term, none to remove all the associated interactions
    - POLAR-12-SCALE: scaling factor for permanent electric fields of polarization between 1-2 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles')
    - POLAR-13-SCALE: scaling factor for permanent electric fields of polarization between 1-3 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles')
    - POLAR-14-SCALE: scaling factor for permanent electric fields of polarization between 1-4 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles')
    - POLAR-15-SCALE: scaling factor for permanent electric fields of polarization between 1-5 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles')
    - POLAR-12-INTRA: scaling factor for permanent electric fields of polarization between 1-2 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles'), within polarization group
    - POLAR-13-INTRA: scaling factor for permanent electric fields of polarization between 1-3 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles'), within polarization group
    - POLAR-14-INTRA: scaling factor for permanent electric fields of polarization between 1-4 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles'), within polarization group
    - POLAR-15-INTRA: scaling factor for permanent electric fields of polarization between 1-5 connected atoms (scaling factor for 'p' permanent electric field, 'p-induced dipoles'), within polarization group
    - DIRECT-11-SCALE: scaling factor for permanent electric fields of polarization between 1-2 connected atoms (scaling factor for 'd' permanent electric field, 'd-induced dipoles')
    - DIRECT-12-SCALE: scaling factor for permanent electric fields of polarization between 1-3 connected atoms (scaling factor for 'd' permanent electric field, 'd-induced dipoles')
    - DIRECT-13-SCALE: scaling factor for permanent electric fields of polarization between 1-4 connected atoms (scaling factor for 'd' permanent electric field, 'd-induced dipoles')
    - DIRECT-14-SCALE: scaling factor for permanent electric fields of polarization between 1-5 connected atoms (scaling factor for 'd' permanent electric field, 'd-induced dipoles')
    - MUTUAL-11-SCALE: scaling factor for electric fields of current induced dipoles between 1-2 connected atoms
    - MUTUAL-12-SCALE: scaling factor for electric fields of current induced dipoles between 1-3 connected atoms
    - MUTUAL-13-SCALE: scaling factor for electric fields of current induced dipoles between 1-4 connected atoms
    - MUTUAL-14-SCALE: scaling factor for electric fields of current induced dipoles between 1-5 connected atoms
    - INDUCE-12-SCALE: scaling factor for electric fields of current induced dipoles between 1-2 connected atoms, with charge penetration
    - INDUCE-13-SCALE: scaling factor for electric fields of current induced dipoles between 1-3 connected atoms, with charge penetration
    - INDUCE-14-SCALE: scaling factor for electric fields of current induced dipoles between 1-4 connected atoms, with charge penetration
    - INDUCE-15-SCALE: scaling factor for electric fields of current induced dipoles between 1-5 connected atoms, with charge penetration
    - POLAR-ALG: iterative algorithm to solve the polarization equations, 1: PCG (default), 2: Jacobi/DIIS, 3: TCG, 5: DC-Jacobi/DIIS
    - POLAR-ALGSHORT: iterative algorithm to solve the short-range polarization equations, 1: PCG (default), 2: Jacobi/DIIS, 3: TCG, 5: DC-Jacobi/DIIS
    - POLAR-EPS: convergence threshold for polarization equations
    - POLAR-PRT: printing level of the polarization solver
    - TCGORDER: order of TCG (default 2)
    - TCGGUESS: direct guess for TCG (default 0, false)
    - TCGPEEK: peek step for TCG (default 1, true)


- **Charge-Transfer**:
  - description: exponential-based charge transfer between multipole sites
  - keywords:
    - CHGTRNTERM: only to restrict the potential to charge-transfer term, none to remove all the associated interactions
    - CHGTRN-CUTOFF: cutoff for charge-transfer interactions
    - CHGTRNSHORT-CUTOFF: cutoff for short-range charge-transfer interactions

- **Restraints**:
  - description: harmonic restraints on various geometric quantities (positions of atoms, distances between atoms...). Note that a wide variety of harmonic restraints can be designed through the Colvars interface when this library is linked with Tinker-HP.
  - keywords:
    - RESTRAINTERM: only to restrict the potential to geometrical restrain term, none to remove all the associated interactions
    - RESTRAIN-POSITION: positional restraints, takes an atom or a list of atoms (-a b for a list of atoms between index a and b) as an argument. Without further indication the atoms are restrained to their original position with a force constant of 100 kcal/mol/A**2
    - RESTRAIN-BACKBONE: restrain all atoms whose name is 'CA' to their original position with a force constant of 100 kcal/mol/A**2
    - RESTRAIN-DISTANCE: takes to atom indexes as an argument and 2 values, the first is the value of the spring constant in kcal/mol/A**2 and the other is the distances after which a flat bottom harmonic restraint starts
    - RESTRAIN-ANGLE: takes 3 atom indexes defining an angle as an argument and 2 values, same logic as for distances
    - RESTRAIN-TORSION: takes 4 atom indexes defining a dihedral as an argument and 2 values, same logic as for distances
  - files: kgeom.f90, egeom.f90, egeom1.f90, egeom3.f90

- **Long-range electrostatics with Smooth Particle Mesh Ewald**:
  - description: Tinker-HP only handles Periodic Boundary Conditions where electrostatics interactions are handled through Smooth Particle-Mesh Ewald. This concerns charge-charge and multipole-multipole interactions as well as polarization and dispersion with ewald summation. Long-range intercations require 3d-fft computations that are handled in parallel within the 2decomp_fft library.
  - keywords:
    - EWALD: use ewald summation for electrostatics (only option available to deal with electrostatics)
    - EWALD-CUTOFF: cutoff for real space Ewald (default 7 Angstroms)
    - EWALD-ALPHA: alpha parameter for Ewald summation, default corresponds to an error of 10d-8 for real space with the real space cutoff
    - EWALDSHORT-CUTOFF: cutoff for short range real space Ewald (default 5 Angstroms)
    - EWALD-BOUNDARY: type of boundary for Ewald summation, can be 'TINFOIL' (default) or 'VACUUM'
    - PME-GRID: number of points in each axis of the 3d PME grid, chosen by default so that the density of the grid corresponds to a low error of Ewald summation with a 7 Angstroms real space cutoff
    - PME-ORDER: order of the splines for reciprocal part of Ewald summation with SPME, default 5
    - DEWALD: use ewald summation for dispersion
    - DEWALD-CUTOFF: cutoff for real space ewald summation for dispersion
    - DEWALDSHORT-CUTOFF: cutoff for short-range real space ewald summation for dispersion
    - DEWALD-ALPHA: alpha parameter for Ewald summation for dispersion
    - DPME-GRID: same as PME-GRID but for dispersion
    - DPME-ORDER: same as PME-ORDER but for dispersion
  - files: kewald.f90, pmestuff.f90, fft_mpi.f90, 2decomp_fft routines


