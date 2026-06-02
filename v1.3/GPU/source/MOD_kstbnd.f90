!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module kstbnd  --  forcefield parameters for stretch-bend  ##
!     ##                                                             ##
!     #################################################################
!
!
!     maxnsb   maximum number of stretch-bend parameter entries
!
!     stbn     force constant parameters for stretch-bend terms
!     ksb      integer signature of atom classes for stretch-bend terms
!     ksb_sys  integer signature of atom classes for stretch-bend terms of the system simulated
!
!
#include "tinker_macro.h"
module kstbnd
   implicit none
   integer maxnsb
   parameter (maxnsb=2000)
   integer(8) ksb(maxnsb)
   integer(8) ksb_sys(0:maxnsb)
   real(t_p) stbn(2,maxnsb)
   save
!$acc declare create(ksb,ksb_sys)
end
