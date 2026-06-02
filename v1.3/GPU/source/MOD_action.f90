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
!     neb     number of bond stretch energy terms computed
!     nea     number of angle bend energy terms computed
!     neba    number of stretch-bend energy terms computed
!     neub    number of Urey-Bradley energy terms computed
!     neaa    number of angle-angle energy terms computed
!     neopb   number of out-of-plane bend energy terms computed
!     neopd   number of out-of-plane distance energy terms computed
!     neid    number of improper dihedral energy terms computed
!     neit    number of improper torsion energy terms computed
!     net     number of torsional energy terms computed
!     nept    number of pi-orbital torsion energy terms computed
!     neat    number of angle-torsion energy terms computed
!     nebt    number of stretch-torsion energy terms computed
!     nett    number of torsion-torsion energy terms computed
!     nev     number of van der Waals energy terms computed
!     ner     number of Pauli repulsion energy terms computed
!     nedsp   number of dispersion energy terms computed
!     nec     number of charge-charge energy terms computed
!     nem     number of multipole energy terms computed
!     nem_    double precision container of number of multipole energy terms computed
!     nect    number of charge transfer energy terms computed
!     nep     number of polarization energy terms computed
!     neg     number of geometric restraint energy terms computed
!     nex     number of extra energy terms computed
!
!
module action
   implicit none
   integer neb,nea,neba,neub
   integer neaa,neopb,neopd
   integer neid,neit,net,nept
   integer neat,nebt,nett,nev,nec
   integer ner,nedsp,nect
   integer nem,nep
   integer neg,nex
   integer nemlpot
   logical :: action_data_ondevice=.FALSE.
   real*8 nem_,nep_,nev_,nec_

contains
   subroutine create_action_data_ondevice
!$acc enter data create(nev,ner,nedsp,nec,nem,nep,nect &
!$acc     ,nem_,nep_,nev_,nec_,nemlpot)
      action_data_ondevice=.TRUE.
   end subroutine
   subroutine delete_action_data_ondevice
!$acc exit data delete(nev,ner,nedsp,nec,nem,nep,nect &
!$acc    ,nem_,nep_,nev_,nec_,nemlpot)
      action_data_ondevice=.FALSE.
   end subroutine
end
