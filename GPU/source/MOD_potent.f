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
c     use_ctransfer  logical flag governing use of ect formula
c     use_dispersion logical flag governing use of dispersion formula
c     use_repulsion  logical flag governing use of repulsion formula
c     use_mreal      logical flag governing use of real space multipolar potential
c     use_mrec       logical flag governing use of reciprocal space multipolar potential
c     use_mself      logical flag governing use of self multipolar potential
c     use_creal      logical flag governing use of real space charge potential
c     use_crec       logical flag governing use of reciprocal space charge potential
c     use_cself      logical flag governing use of self charge potential
c     use_preal      logical flag governing use of real space polarization potential
c     use_prec       logical flag governing use of reciprocal space polarization potential
c     use_pself      logical flag governing use of self polarization potential
c     use_smd_velconst     logical flag governing use of CVSMD
c     use_smd_forconst     logical flag governing use of CFSMD
c     use_amd_dih       allow for the use of aMD on only dihedrals
c     use_amd_ene       allow for the use of aMD on only potential energy
c     use_amdstep       allow for the use of aMD between first and last step
c
c     PotentialAll    logical flag enable with all potential (default)
c     PotentialAmoeba* logical flag enable when Amoeba forcefield is being processed
c     PotentialCharmm logical flag enable with Charmm forcefield
c     PotentialWater* logical flag enable with Water's main potential terms
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
      logical use_tors,use_pitors,use_strtor
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
!$acc declare create(use_mpole,use_polar)
      end module

      module potent
      use amdpotent
      use bondedpotent
      use nonbondedpotent
      implicit none
      logical use_born,use_pmecore
      logical use_emtp,use_ctransfer,use_dispersion
      logical use_repulsion

      logical use_smd_velconst, use_smd_forconst

      logical PotentialAll
      logical PotentialAmoeba,PotentialAmoeba18
      logical PotentialWaterAmoeba,PotentialWaterCharmm
      logical PotentialCharmm
!$acc declare create(use_pmecore)
      end
