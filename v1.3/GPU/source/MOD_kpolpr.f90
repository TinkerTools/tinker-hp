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
!     thlpr    Thole damping values for special polarization pairs
!     thdpr    Thole direct damping for special polarization pairs
!     kppr     string of atom types for special polarization pairs
!
#include "tinker_macro.h"
module kpolpr
   use sizes
   implicit none
   integer maxnpp
   real*8 :: thlpr(maxtyp)
   real*8 :: thdpr(maxtyp)
   parameter (maxnpp=100)
   character*8 :: kppr(maxnpp)
end
