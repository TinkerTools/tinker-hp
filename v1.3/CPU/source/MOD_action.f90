!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module action  -- total number of each energy term computed  ##
!     ##                                                               ##
!     ###################################################################
!
!
module action
   implicit none
   integer :: neb !<number of bond stretch energy terms computed
   integer :: nea !<number of angle bend energy terms computed
   integer :: neba !<number of stretch-bend energy terms computed
   integer :: neub !<number of Urey-Bradley energy terms computed
   integer :: neaa !<number of angle-angle energy terms computed
   integer :: neopb !<number of out-of-plane bend energy terms computed
   integer :: neopd !<number of out-of-plane distance energy terms computed
   integer :: neid !<number of improper dihedral energy terms computed
   integer :: neit !<number of improper torsion energy terms computed
   integer :: net !<number of torsional energy terms computed
   integer :: nept !<number of pi-orbital torsion energy terms computed
   integer :: neat !<number of angle-torsion energy terms computed
   integer :: nebt !<number of stretch-torsion energy terms computed
   integer :: nett !<number of torsion-torsion energy terms computed
   integer :: nev !<number of van der Waals energy terms computed
   integer :: nec !<number of charge-charge energy terms computed
   integer :: ner !<number of Pauli repulsion energy terms computed
   integer :: nedsp !<number of dispersion energy terms computed
   integer :: nect !<number of charge transfer energy terms computed
   integer :: nem !<number of multipole energy terms computed
   integer :: nep !<number of polarization energy terms computed
   integer :: neg !<number of geometric restraint energy terms computed
   integer :: nex !<number of extra energy terms computed
   save
end
