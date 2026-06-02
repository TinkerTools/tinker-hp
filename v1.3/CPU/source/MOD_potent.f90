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
module potent
   implicit none
   logical :: use_bond !<logical flag governing use of bond stretch potential
   logical :: use_angle !<logical flag governing use of angle bend potential
   logical :: use_strbnd !<logical flag governing use of stretch-bend potential
   logical :: use_urey !<logical flag governing use of Urey-Bradley potential
   logical :: use_angang !<logical flag governing use of angle-angle cross term
   logical :: use_opbend !<logical flag governing use of out-of-plane bend term
   logical :: use_opdist !<logical flag governing use of out-of-plane distance
   logical :: use_improp !<logical flag governing use of improper dihedral term
   logical :: use_imptor !<logical flag governing use of improper torsion term
   logical :: use_tors !<logical flag governing use of torsional potential
   logical :: use_pitors !<logical flag governing use of pi-orbital torsion term
   logical :: use_strtor !<logical flag governing use of stretch-torsion term
   logical :: use_angtor !<logical flag governing use of angle-torsion potential
   logical :: use_tortor !<logical flag governing use of torsion-torsion term
   logical :: use_vdw !<logical flag governing use of vdw der Waals potential
   logical :: use_charge !<logical flag governing use of charge-charge potential
   logical :: use_mpole !<logical flag governing use of multipole potential
   logical :: use_polar !<logical flag governing use of polarization term
   logical :: use_solv !<logical flag governing use of implicit solvation
   logical :: use_geom !<logical flag governing use of geometric restraints
   logical :: use_extra !<logical flag governing use of extra potential term
   logical :: use_pmecore !<logical flag governing use of separate cores for pme
   logical :: use_grad !<logical flag being true only if the main is "dynamic"
   logical :: use_mreal !<logical flag governing use of real space multipolar potential
   logical :: use_mrec !<logical flag governing use of reciprocal space multipolar potential
   logical :: use_preal !<logical flag governing use of real space polarization potential
   logical :: use_prec !<logical flag governing use of reciprocal space polarization potential
   logical :: use_creal !<logical flag governing use of real space charge potential
   logical :: use_crec !<logical flag governing use of reciprocal space charge potential
   logical :: use_cself !<logical flag governing use of self charge potential
   logical :: use_polarshortreal !<logical flag governing use of short range real space polarization term
   logical :: use_mself !<logical flag governing use of self multipolar potential
   logical :: use_pself !<logical flag governing use of self polarization potential
   logical :: use_mpoleshortreal !<logical flag governing use of short range real space multipole potential
   logical :: use_cshortreal !<logical flag governing use of short range real space charge potential
   logical :: use_clong !<logical flag governing use of long range real space charge potential
   logical :: use_vdwshort !<logical flag governing use of short range vdw potential
   logical :: use_vdwlong !<logical flag governing use of long range vdw potential
   logical :: use_mpolelong !<logical flag governing use of long range real space multipole potential
   logical :: use_smd_velconst !<logical flag governing use of CVSMD
   logical :: use_smd_forconst !<logical flag governing use of CFSMD
   logical :: use_repuls !<logical flag governing use of Pauli repulsion term
   logical :: use_disp !<logical flag governing use of dispersion potential
   logical :: use_chgtrn !<logical flag governing use of charge transfer term
   logical :: use_chgflx !<logical flag governing use of charge flux term
   logical :: use_dispreal !<logical flag governing use of real space polarization potential
   logical :: use_dispself !<logical flag governing use of self polarization potential
   logical :: use_disprec !<logical flag governing use of reciprocal space polarization potential
   logical :: use_repulsshort !<logical flag governing use of short range vdw potential
   logical :: use_dispshort !<logical flag governing use of short range dispersion potential
   logical :: use_dispshortreal !<logical flag governing use of short range real space dispersion potential
   logical :: use_chgtrnshort !<logical flag governing use of short range charge transfer term
   logical :: use_repulslong !<logical flag governing use of long range vdw potential
   logical :: use_displong !<logical flag governing use of long range dispersion potential
   logical :: use_chgtrnlong !<logical flag governing use of long range charge transfer term
   logical :: use_dewald !<logical flag governing use of PME for dispersion
   logical :: use_chgpen !<logical flag governing use of charge penetration
   logical :: use_lambdadyn !<logical flag governing use of lambda dynamic (with colvar module)
   logical :: use_OSRW !<logical flag governing us of Orthogonal Space Random Walk sampling (with colvar module)

   save
end
