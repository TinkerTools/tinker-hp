!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module potent  --  usage of each potential energy component  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     use_bond       logical flag governing use of bond stretch potential
!     use_angle      logical flag governing use of angle bend potential
!     use_angtor     logical flag governing use of angle-torsion potential
!     use_strbnd     logical flag governing use of stretch-bend potential
!     use_urey       logical flag governing use of Urey-Bradley potential
!     use_angang     logical flag governing use of angle-angle cross term
!     use_opbend     logical flag governing use of out-of-plane bend term
!     use_opdist     logical flag governing use of out-of-plane distance
!     use_improp     logical flag governing use of improper dihedral term
!     use_imptor     logical flag governing use of improper torsion term
!     use_tors       logical flag governing use of torsional potential
!     use_pitors     logical flag governing use of pi-orbital torsion term
!     use_strtor     logical flag governing use of stretch-torsion term
!     use_tortor     logical flag governing use of torsion-torsion term
!     use_embd_potoff logical flag governing use of potential in embedding
!     use_embd_bond  logical flag governing use of bond stretch potential in embedding
!     use_embd_angle logical flag governing use of angle bend potential in embedding
!     use_embd_strbnd logical flag governing use of stretch-bend potential in embedding
!     use_embd_urey  logical flag governing use of Urey-Bradley potential in embedding
!     use_embd_angang logical flag governing use of angle-angle cross term in embedding
!     use_embd_opbend logical flag governing use of out-of-plane bend term in embedding
!     use_embd_opdist logical flag governing use of out-of-plane distance in embedding
!     use_embd_improp logical flag governing use of improper dihedral term in embedding
!     use_embd_imptor logical flag governing use of improper torsion term in embedding
!     use_embd_tors  logical flag governing use of torsional potential in embedding
!     use_embd_pitors logical flag governing use of pi-orbital torsion term in embedding
!     use_embd_strtor logical flag governing use of stretch-torsion term in embedding
!     use_embd_tortor logical flag governing use of torsion-torsion term in embedding
!     use_vdw        logical flag governing use of vdw der Waals potential
!     use_vdwshort   logical flag governing use of short range vdw potential
!     use_vdwlong    logical flag governing use of long range vdw potential
!     use_charge     logical flag governing use of charge-charge potential
!     use_mpole      logical flag governing use of multipole potential
!     use_mpoleshortreal   logical flag governing use of short range real space multipole potential
!     use_mpolelong  logical flag governing use of long range real space multipole potential
!     use_cshortreal logical flag governing use of short range real space charge potential
!     use_clong      logical flag governing use of long range real space charge potential
!     use_polar      logical flag governing use of polarization term
!     use_polarshortreal   logical flag governing use of short range real space polarization term
!     use_geom       logical flag governing use of geometric restraints
!     use_extra      logical flag governing use of extra potential term
!     use_pmecore    logical flag governing use of separate cores for pme
!     use_emtp       logical flag governing use of emtp formula for electrostatics
!     use_mreal      logical flag governing use of real space multipolar potential
!     use_mrec       logical flag governing use of reciprocal space multipolar potential
!     use_mself      logical flag governing use of self multipolar potential
!     use_creal      logical flag governing use of real space charge potential
!     use_crec       logical flag governing use of reciprocal space charge potential
!     use_cself      logical flag governing use of self charge potential
!     use_preal      logical flag governing use of real space polarization potential
!     use_prec       logical flag governing use of reciprocal space polarization potential
!     use_pself      logical flag governing use of self polarization potential
!     use_dispreal   logical flag governing use of real space polarization potential
!     use_disprec    logical flag governing use of reciprocal space polarization potential
!     use_dispself   logical flag governing use of self polarization potential
!     use_smd_velconst  logical flag governing use of CVSMD
!     use_smd_forconst  logical flag governing use of CFSMD
!     use_repuls        logical flag governing use of Pauli repulsion term
!     use_repulsshort   logical flag governing use of short range vdw potential
!     use_repulslong    logical flag governing use of long range vdw potential
!     use_disp          logical flag governing use of dispersion potential
!     use_dispshort     logical flag governing use of short range dispersion potential
!     use_dispshortreal logical flag governing use of short range real space dispersion potential
!     use_dispslong     logical flag governing use of long range dispersion potential
!     use_disp          logical flag governing use of dispersion potential
!     use_chgtrnshort   logical flag governing use of short range charge transfer term
!     use_chgtrnlong    logical flag governing use of long range charge transfer term
!     use_chgtrn        logical flag governing use of charge transfer term
!     use_chgflx        logical flag governing use of charge flux term
!     use_dewald        logical flag governing use of PME for dispersion
!     use_chgpen        logical flag governing use of charge penetration
!     use_lambdadyn     logical flag governing use of lambda dynamic (with colvar module)
!     use_OSRW          logical flag governing us of Orthogonal Space Random Walk sampling (with colvar module)
!     use_amd_dih       allow for the use of aMD on only dihedrals
!     use_amd_ene       allow for the use of aMD on only potential energy
!     use_amdstep       allow for the use of aMD between first and last step
!
!     fuse_chglj       logical flag enable with Lennard-Jones and point-charge potentials
!     fuse_bonded      logical flag to compute all bonded interactions in a single device kernel
!     bonded_l         logical flag that follow bonded computation
!     shortnonbonded_l logical flag that follow short range non bonded computation
!     nonbonded_l      logical flag that follow non bonded computation
!     PotentialAll     logical flag enable with all potential (default)
!     PotentialAmoeba* logical flag enable when Amoeba forcefield is being processed
!     PotentialCharmm  logical flag enable with Charmm forcefield
!     PotentialWater*  logical flag enable with Water's main potential terms
!     use_mlpot         flag enable when ML potential is being processed
!
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
   logical use_embd_potoff
   logical use_embd_bond,use_embd_angle,use_embd_strbnd
   logical use_embd_urey,use_embd_angang,use_embd_opbend
   logical use_embd_opdist,use_embd_improp,use_embd_imptor
   logical use_embd_tors,use_embd_pitors,use_embd_strtor
   logical use_embd_angtor,use_embd_tortor
   logical use_embd_geom,use_embd_extra
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
   !logical use_ctransfer,use_dispersion,use_repulsion
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
   logical use_mlpot,use_ml_embedding,use_mlpot_only
   logical use_lambdadyn
   logical use_OSRW

   logical fuse_chglj
   logical fuse_bonded, disable_fuse_bonded

   logical PotentialAll
   logical PotentialAmoeba,PotentialAmoeba18
   logical PotentialAmoeba181,PotentialAmoeba182
   logical PotentialWaterAmoeba,PotentialWaterCharmm
   logical PotentialCharmm
!$acc declare create(use_pmecore)
end
