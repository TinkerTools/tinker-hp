!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module chgtrn  --  charge transfer for current structure  ##
!     ##                                                            ##
!     ################################################################
!
!
!     nct       total number of dispersion sites in the system
!     chgct     charge for charge transfer at each multipole site
!     dmpct     charge transfer damping factor at each multipole site
!     winchgct  window associated to chgct array
!     windmpct  window associated to dmpct array
!
!
#include "tinker_macro.h"
module chgtrn
   implicit none
   integer nct
   real(t_p), pointer :: chgct(:)
   real(t_p), pointer :: dmpct(:)
   integer :: winchgct,windmpct
   save
end
