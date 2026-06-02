!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module energi  --  individual potential energy components  ##
!     ##                                                             ##
!     #################################################################
!
!
!     esum   !<total potential energy of the system
!     eb     !<bond stretch potential energy of the system
!     ea     !<angle bend potential energy of the system
!     eba    !<stretch-bend potential energy of the system
!     eub    !<Urey-Bradley potential energy of the system
!     eaa    !<angle-angle potential energy of the system
!     eopb   !<out-of-plane bend potential energy of the system
!     eopd   !<out-of-plane distance potential energy of the system
!     eid    !<improper dihedral potential energy of the system
!     eit    !<improper torsion potential energy of the system
!     et     !<torsional potential energy of the system
!     ept    !<pi-orbital torsion potential energy of the system
!     eat    !<angle-torsion potential energy of the system
!     ebt    !<stretch-torsion potential energy of the system
!     ett    !<torsion-torsion potential energy of the system
!     ev     !<van der Waals potential energy of the system
!     er     !<Pauli repulsion potential energy of the system
!     edsp   !<dispersion potential energy of the system
!     ec     !<charge-charge potential energy of the system
!     em     !<atomic multipole potential energy of the system
!     ep     !<polarization potential energy of the system
!     ect    !<charge transfer potential energy of the system
!     eg     !<geometric restraint potential energy of the system
!     ex     !<extra term potential energy of the system
!     esave  !<stored potential energy of the system
!     ensmd  !<extra term smd potential energy of the system
!
!
module energi
   implicit none
   real*8 :: esum
   real*8 :: eb
   real*8 :: ea
   real*8 :: eba
   real*8 :: eub
   real*8 :: eaa
   real*8 :: eopb
   real*8 :: eopd
   real*8 :: eid
   real*8 :: eit
   real*8 :: et
   real*8 :: ept
   real*8 :: eat
   real*8 :: ebt
   real*8 :: ett
   real*8 :: ev
   real*8 :: ec
   real*8 :: em
   real*8 ::ep
   real*8 :: er
   real*8 :: edsp
   real*8 :: ect
   real*8 :: eg
   real*8 :: ex
   real*8 :: esave
   real*8 :: ensmd
   save
end

!> @brief 
!> outputs the components of the current energy
!> @param no params
subroutine info_energy(rank)
   use energi
   use inform
   use iounit
   implicit none
   integer,intent(in):: rank
   real(8) ebonded

   if (rank.eq.0) then
20    format ( 40('-'))
      print 20
30    format (1x,A,F18.6)

      ebonded =&
      &eb+ea+eba+eub+eaa+eid+eit+et+ept+ebt+ett+eat+eopb+eopd+eg+ex
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
