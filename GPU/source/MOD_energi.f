c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module energi  --  individual potential energy components  ##
c     ##                                                             ##
c     #################################################################
c
c
c     esum   total potential energy of the system
c     eb     bond stretch potential energy of the system
c     ea     angle bend potential energy of the system
c     eba    stretch-bend potential energy of the system
c     eub    Urey-Bradley potential energy of the system
c     eaa    angle-angle potential energy of the system
c     eopb   out-of-plane bend potential energy of the system
c     eopd   out-of-plane distance potential energy of the system
c     eid    improper dihedral potential energy of the system
c     eit    improper torsion potential energy of the system
c     et     torsional potential energy of the system
c     ept    pi-orbital torsion potential energy of the system
c     eat    angle-torsion potential energy of the system
c     ebt    stretch-torsion potential energy of the system
c     ett    torsion-torsion potential energy of the system
c     ev     van der Waals potential energy of the system
c     ec     charge-charge potential energy of the system
c     ecrec  reciproqual charge-charge potential energy of the system
c     em     atomic multipole potential energy of the system
c     emrec  reciproqual part of atomic multipole potential energy of the system
c     ep     polarization potential energy of the system
c     eprec  reciproqual part of polarization potential energy of the system
c     eg     geometric restraint potential energy of the system
c     ex     extra term potential energy of the system
c     esave  stored potential energy of the system
c     ensmd  extra term smd potential energy of the system
c     eDaMD  extra term aMD potential for dihedrals of the system
c     ePaMD  extra term aMD potential energy of the system
c     eW1aMD extra term aMD potential for bonds/angles of waters 
c     eW2aMD extra term aMD potential energy of waters
c
c
#include "tinker_precision.h"
#include "tinker_types.h"
      module energi
      implicit none
      real(r_p) esum,eb,ea,eba
      real(r_p) eub,eaa,eopb,eopd
      real(r_p) eid,eit,et,ept
      real(r_p) eat,ebt,ett
      real(r_p) ev,ec,ecrec,em,emrec,ep,eprec
      real(r_p) eg,ex
      real(r_p) esave
      real(r_p) ensmd
      real(r_p) eDaMD,ePaMD,eW1aMD,eW2aMD
      ener_rtyp ev_r,ec_r,em_r,ep_r,eb_r

      interface
      module subroutine create_energi_on_device()
      end subroutine
      module subroutine delete_energi_on_device()
      end subroutine
      module subroutine info_energy(rank)
      integer,intent(in)::rank
      end subroutine
      end interface

      end module

      submodule(energi) subEnergi

      contains
#include "convert.f.inc"

      module subroutine create_energi_on_device()
      implicit none
!$acc enter data create(esum,eb,ea,eba,eub,eaa,eopb,eopd,
!$acc&                  eid,eit,et,ept,ebt,ett,eg,ex,eat,
!$acc&                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd,
!$acc&                  eDaMD,ePaMD,eW1aMD,eW2aMD,
!$acc&                  ev_r,ec_r,em_r,ep_r,eb_r)
      end subroutine
      module subroutine delete_energi_on_device()
      implicit none
!$acc exit data delete(esum,eb,ea,eba,eub,eaa,eopb,eopd,
!$acc&                  eid,eit,et,ept,ebt,ett,eg,ex,eat,
!$acc&                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd,
!$acc&                  eDaMD,ePaMD,eW1aMD,eW2aMD,
!$acc&                  ev_r,ec_r,em_r,ep_r,eb_r)
      end subroutine

      module subroutine info_energy(rank)
      implicit none
      integer,intent(in):: rank

      if (rank.eq.0) then
!$acc wait
 20   format ( 40('-'))
         print 20
 30   format (1x,A,F18.6)
!$acc update host(esum,eb,ea,eba,eub,eaa,eopb,eopd,
!$acc&                  eid,eit,et,ept,ebt,ett,eg,ex,eat,
!$acc&                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd,
!$acc&                  eDaMD,ePaMD,eW1aMD,eW2aMD,
!$acc&                  ev_r,ec_r,em_r,ep_r,eb_r)

         if (eb   /=real(0,r_p)) print 30, 'eb   = ',eb
         if (ea   /=real(0,r_p)) print 30, 'ea   = ',ea 
         if (eba  /=real(0,r_p)) print 30, 'eba  = ',eba
         if (eub  /=real(0,r_p)) print 30, 'eub  = ',eub
         if (eaa  /=real(0,r_p)) print 30, 'eaa  = ',eaa
         if (eid  /=real(0,r_p)) print 30, 'eid  = ',eid
         if (eit  /=real(0,r_p)) print 30, 'eit  = ',eit
         if (et   /=real(0,r_p)) print 30, 'et   = ',et
         if (ept  /=real(0,r_p)) print 30, 'ept  = ',ept
         if (ebt  /=real(0,r_p)) print 30, 'ebt  = ',ebt
         if (ett  /=real(0,r_p)) print 30, 'ett  = ',ett
         if (eat  /=real(0,r_p)) print 30, 'eat  = ',ett
         if (eopb /=real(0,r_p)) print 30, 'eopb = ',eopb
         if (eopd /=real(0,r_p)) print 30, 'eopd = ',eopd
         if (eg   /=real(0,r_p)) print 30, 'eg   = ',eg
         if (ex   /=real(0,r_p)) print 30, 'ex   = ',eg
         if (ec   /=real(0,r_p)) print 30, 'ec   = ',ec
         if (ev   /=real(0,r_p)) print 30, 'ev   = ',ev
         if (em   /=real(0,r_p)) print 30, 'em   = ',em
         if (ep   /=real(0,r_p)) print 30, 'ep   = ',ep
         if (ensmd/=real(0,r_p)) print 30, 'ensmd= ',ensmd
         if (eDaMD/=real(0,r_p)) print 30, 'eDaMD= ',eDaMD
         if (ePaMD/=real(0,r_p)) print 30, 'ePaMD= ',ePaMD
         if(eW1aMD/=real(0,r_p)) print 30,'eW1aMD= ',eW1aMD
         !if (esum/=real(0,r_p)) print 30, 'esum = ',esum
      end if
      end subroutine

      end submodule
