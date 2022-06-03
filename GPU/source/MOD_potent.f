c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module potent  --  usage of each potential energy component  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     use_bond       logical flag governing use of bond stretch potential
c     use_angle      logical flag governing use of angle bend potential
c     use_angtor     logical flag governing use of angle-torsion potential
c     use_strbnd     logical flag governing use of stretch-bend potential
c     use_urey       logical flag governing use of Urey-Bradley potential
c     use_angang     logical flag governing use of angle-angle cross term
c     use_opbend     logical flag governing use of out-of-plane bend term
c     use_opdist     logical flag governing use of out-of-plane distance
c     use_improp     logical flag governing use of improper dihedral term
c     use_imptor     logical flag governing use of improper torsion term
c     use_tors       logical flag governing use of torsional potential
c     use_pitors     logical flag governing use of pi-orbital torsion term
c     use_strtor     logical flag governing use of stretch-torsion term
c     use_tortor     logical flag governing use of torsion-torsion term
c     use_vdw        logical flag governing use of vdw der Waals potential
c     use_vdwshort   logical flag governing use of short range vdw potential
c     use_vdwlong    logical flag governing use of long range vdw potential
c     use_charge     logical flag governing use of charge-charge potential
c     use_mpole      logical flag governing use of multipole potential
c     use_mpoleshortreal   logical flag governing use of short range real space multipole potential
c     use_mpolelong  logical flag governing use of long range real space multipole potential
c     use_cshortreal logical flag governing use of short range real space charge potential
c     use_clong      logical flag governing use of long range real space charge potential
c     use_polar      logical flag governing use of polarization term
c     use_polarshortreal   logical flag governing use of short range real space polarization term
c     use_geom       logical flag governing use of geometric restraints
c     use_extra      logical flag governing use of extra potential term
c     use_pmecore    logical flag governing use of separate cores for pme
c     use_emtp       logical flag governing use of emtp formula for electrostatics
c     use_mreal      logical flag governing use of real space multipolar potential
c     use_mrec       logical flag governing use of reciprocal space multipolar potential
c     use_mself      logical flag governing use of self multipolar potential
c     use_creal      logical flag governing use of real space charge potential
c     use_crec       logical flag governing use of reciprocal space charge potential
c     use_cself      logical flag governing use of self charge potential
c     use_preal      logical flag governing use of real space polarization potential
c     use_prec       logical flag governing use of reciprocal space polarization potential
c     use_pself      logical flag governing use of self polarization potential
c     use_dispreal   logical flag governing use of real space polarization potential
c     use_disprec    logical flag governing use of reciprocal space polarization potential
c     use_dispself   logical flag governing use of self polarization potential
c     use_smd_velconst  logical flag governing use of CVSMD
c     use_smd_forconst  logical flag governing use of CFSMD
c     use_repuls        logical flag governing use of Pauli repulsion term
c     use_repulsshort   logical flag governing use of short range vdw potential
c     use_repulslong    logical flag governing use of long range vdw potential
c     use_disp          logical flag governing use of dispersion potential
c     use_dispshort     logical flag governing use of short range dispersion potential
c     use_dispshortreal logical flag governing use of short range real space dispersion potential
c     use_dispslong     logical flag governing use of long range dispersion potential
c     use_disp          logical flag governing use of dispersion potential
c     use_chgtrnshort   logical flag governing use of short range charge transfer term
c     use_chgtrnlong    logical flag governing use of long range charge transfer term
c     use_chgtrn        logical flag governing use of charge transfer term
c     use_chgflx        logical flag governing use of charge flux term
c     use_dewald        logical flag governing use of PME for dispersion
c     use_chgpen        logical flag governing use of charge penetration
c     use_lambdadyn     logical flag governing use of lambda dynamic (with colvar module)
c     use_OSRW          logical flag governing us of Orthogonal Space Random Walk sampling (with colvar module)
c     use_amd_dih       allow for the use of aMD on only dihedrals
c     use_amd_ene       allow for the use of aMD on only potential energy
c     use_amdstep       allow for the use of aMD between first and last step
c
c     fuse_chglj       logical flag enable with Lennard-Jones and point-charge potentials 
c     fuse_bonded      logical flag to compute all bonded interactions in a single device kernel
c     bonded_l         logical flag that follow bonded computation
c     shortnonbonded_l logical flag that follow short range non bonded computation
c     nonbonded_l      logical flag that follow non bonded computation
c     PotentialAll     logical flag enable with all potential (default)
c     PotentialAmoeba* logical flag enable when Amoeba forcefield is being processed
c     PotentialCharmm  logical flag enable with Charmm forcefield
c     PotentialWater*  logical flag enable with Water's main potential terms
c
      module amdpotent
      implicit none
      logical use_amd_dih, use_amd_ene, use_gamd
      logical use_amd_wat1
      end module

      module bondedpotent
      implicit none
      logical use_bond,use_angle,use_strbnd
      logical use_urey,use_angang,use_opbend
      logical use_opdist,use_improp,use_imptor
      logical use_tors,use_pitors,use_strtor,use_angtor
      logical use_tortor
      logical use_geom,use_extra
      end module
      
      module nonbondedpotent
      implicit none
      logical use_vdw
      logical use_vdwshort,use_vdwlong,use_mpolelong
      logical use_charge
      logical use_cshortreal,use_clong
      logical use_creal,use_crec,use_cself
      logical use_mpole
      logical use_mself,use_mpoleshortreal
      logical use_mreal,use_mrec
      logical use_polar,use_solv
      logical use_pself,use_polarshortreal
      logical use_preal,use_prec
      logical use_repuls,use_disp,use_chgtrn,use_chgflx
      logical use_dispreal,use_dispself,use_disprec
      logical use_repulsshort,use_dispshort,use_dispshortreal
      logical use_chgtrnshort
      logical use_repulslong,use_displong,use_chgtrnlong
      logical use_smd_velconst, use_smd_forconst
      logical use_dewald,use_chgpen
!$acc declare create(use_mpole,use_polar)
      end module

      module potent
      use amdpotent
      use bondedpotent
      use nonbondedpotent
      implicit none
      logical use_born,use_pmecore
      logical use_emtp
      logical bonded_l,shortnonBonded_l,nonBonded_l
      logical use_lambdadyn
      logical use_OSRW

      logical fuse_chglj
      logical fuse_bonded

      logical PotentialAll
      logical PotentialAmoeba,PotentialAmoeba18
      logical PotentialAmoeba181,PotentialAmoeba182
      logical PotentialWaterAmoeba,PotentialWaterCharmm
      logical PotentialCharmm
!$acc declare create(use_pmecore)
      end
