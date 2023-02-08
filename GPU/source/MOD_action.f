c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module action  -- total number of each energy term computed  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     neb     number of bond stretch energy terms computed
c     nea     number of angle bend energy terms computed
c     neba    number of stretch-bend energy terms computed
c     neub    number of Urey-Bradley energy terms computed
c     neaa    number of angle-angle energy terms computed
c     neopb   number of out-of-plane bend energy terms computed
c     neopd   number of out-of-plane distance energy terms computed
c     neid    number of improper dihedral energy terms computed
c     neit    number of improper torsion energy terms computed
c     net     number of torsional energy terms computed
c     nept    number of pi-orbital torsion energy terms computed
c     neat    number of angle-torsion energy terms computed
c     nebt    number of stretch-torsion energy terms computed
c     nett    number of torsion-torsion energy terms computed
c     nev     number of van der Waals energy terms computed
c     ner     number of Pauli repulsion energy terms computed
c     nedsp   number of dispersion energy terms computed
c     nec     number of charge-charge energy terms computed
c     nem     number of multipole energy terms computed
c     nem_    double precision container of number of multipole energy terms computed
c     nect    number of charge transfer energy terms computed
c     nep     number of polarization energy terms computed
c     neg     number of geometric restraint energy terms computed
c     nex     number of extra energy terms computed
c
c
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
!$acc enter data create(nev,ner,nedsp,nec,nem,nep,nect
!$acc&     ,nem_,nep_,nev_,nec_,nemlpot)
      action_data_ondevice=.TRUE.
      end subroutine
      subroutine delete_action_data_ondevice
!$acc exit data delete(nev,ner,nedsp,nec,nem,nep,nect
!$acc&    ,nem_,nep_,nev_,nec_,nemlpot)
      action_data_ondevice=.FALSE.
      end subroutine
      end
