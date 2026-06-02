!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module opbend  --  out-of-plane bends in the current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     nopbend   total number of out-of-plane bends in the system
!     nopbendloc   local number of out-of-plane bends in the system
!     iopb      bond angle numbers used in out-of-plane bending
!     winiopb    window object corresponding to iopb
!     npopbend  number of angle used in out-of-plane bending before each atom
!     winnbopbend    window object corresponding to nbopbend
!     opbk      force constant values for out-of-plane bending
!     winopbk    window object corresponding to opbk
!
!
#include "tinker_macro.h"
module opbend
   implicit none
   integer nopbend,nopbendloc
   integer, pointer :: iopb(:), nbopbend(:)
   real(t_p), pointer ::  opbk(:)
   integer :: winiopb,winnbopbend,winopbk
end
