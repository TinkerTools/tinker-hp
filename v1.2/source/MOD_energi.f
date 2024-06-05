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
c     er     Pauli repulsion potential energy of the system
c     edsp   dispersion potential energy of the system
c     ec     charge-charge potential energy of the system
c     em     atomic multipole potential energy of the system
c     ep     polarization potential energy of the system
c     ect    charge transfer potential energy of the system
c     eg     geometric restraint potential energy of the system
c     ex     extra term potential energy of the system
c     esave  stored potential energy of the system
c     ensmd  extra term smd potential energy of the system
c
c
      module energi
      implicit none
      real*8 esum,eb,ea,eba
      real*8 eub,eaa,eopb,eopd
      real*8 eid,eit,et,ept
      real*8 eat,ebt,ett,ev,ec
      real*8 em,ep
      real*8 er,edsp,ect
      real*8 eg,ex
      real*8 esave
      real*8 ensmd
      save
      end

      subroutine info_energy(rank)
      use energi
      use inform
      use iounit
      implicit none
      integer,intent(in):: rank
      real(8) ebonded

      if (rank.eq.0) then
 20   format ( 40('-'))
         print 20
 30   format (1x,A,F18.6)

         ebonded = 
     &   eb+ea+eba+eub+eaa+eid+eit+et+ept+ebt+ett+eat+eopb+eopd+eg+ex
         if (eb   /=real(0,8)) write(iout,30) 'eb     = ',eb
         if (ea   /=real(0,8)) write(iout,30) 'ea     = ',ea 
         if (eba  /=real(0,8)) write(iout,30) 'eba    = ',eba
         if (eub  /=real(0,8)) write(iout,30) 'eub    = ',eub
         if (eaa  /=real(0,8)) write(iout,30) 'eaa    = ',eaa
         if (eid  /=real(0,8)) write(iout,30) 'eid    = ',eid
         if (eit  /=real(0,8)) write(iout,30) 'eit    = ',eit
         if (et   /=real(0,8)) write(iout,30) 'et     = ',et
         if (ept  /=real(0,8)) write(iout,30) 'ept    = ',ept
         if (ebt  /=real(0,8)) write(iout,30) 'ebt    = ',ebt
         if (ett  /=real(0,8)) write(iout,30) 'ett    = ',ett
         if (eat  /=real(0,8)) write(iout,30) 'eat    = ',eat
         if (eopb /=real(0,8)) write(iout,30) 'eopb   = ',eopb
         if (eopd /=real(0,8)) write(iout,30) 'eopd   = ',eopd
         if (eg   /=real(0,8)) write(iout,30) 'eg     = ',eg
         if (ex   /=real(0,8)) write(iout,30) 'ex     = ',ex
         if (ebonded/=real(0,8)) write(iout,30) 'ebonded =',ebonded
         if (ec   /=real(0,8)) write(iout,30) 'ec     = ',ec
         if (ev   /=real(0,8)) write(iout,30) 'ev     = ',ev
         if (er   /=real(0,8)) write(iout,30) 'er     = ',er
         if (edsp /=real(0,8)) write(iout,30) 'edsp   = ',edsp
         if (em   /=real(0,8)) write(iout,30) 'em     = ',em
         if (ep   /=real(0,8)) write(iout,30) 'ep     = ',ep
         if (ect  /=real(0,8)) write(iout,30) 'ect    = ',ect
         if (esum/=real(0,8)) write(iout,30) 'esum    = ',esum
      end if
      end subroutine
