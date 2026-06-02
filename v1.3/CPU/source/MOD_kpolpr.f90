!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  module kpolpr  --  special Thole forcefield parameters  ##
!     ##                                                          ##
!     ##############################################################
!
!     thlpr    !<Thole damping values for special polarization pairs
!     thdpr    !<Thole direct damping for special polarization pairs
!     kppr     !<string of atom types for special polarization pairs
!
module kpolpr
   use sizes
   implicit none
   integer maxnpp  !<maximum number of polarization pair parameter entries
   real*8 :: thlpr(maxtyp) !<Thole damping values for special polarization pairs
   real*8 :: thdpr(maxtyp) !<Thole direct damping for special polarization pairs
   parameter (maxnpp=100)
   character*8 :: kppr(maxnpp) !<string of atom types for special polarization pairs
   save
end
