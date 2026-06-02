!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kopbnd  --  forcefield parameters for out-of-plane bend  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     maxnopb   maximum number of out-of-plane bending entries
!
!     opbn      force constant parameters for out-of-plane bending
!     kopb      signature of atom classes for out-of-plane bending
!     kopb_sys  system signature of atom classes for out-of-plane bending
!
!
#include "tinker_macro.h"
module kopbnd
   implicit none
   integer maxnopb
   parameter (maxnopb=500)
   real(t_p) opbn(maxnopb)
   integer(8) kopb(maxnopb)
   integer(8) kopb_sys(0:maxnopb)
   !DIR$ ATTRIBUTES ALIGN:64:: jopb
   logical, allocatable :: jopb(:)
!$acc declare create(kopb,kopb_sys)
end
