!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #########################################################
!     ##                                                     ##
!     ##  function precise  --  determine machine precision  ##
!     ##                                                     ##
!     #########################################################
!
!
!     "precise" finds a machine precision value as selected by
!     the input argument: (1) the smallest positive floating
!     point value, (2) the smallest relative floating point
!     spacing, (3) the largest relative floating point spacing
!
!
#include "tinker_precision.h"
function precise (i)
   implicit none
   integer,parameter::ti_p=t_p
   integer i
   real(t_p) precise,value
   real(t_p) zero,one,delta
!
!
!     set values for zero, one and multiplicative factor
!
   parameter (&
      &zero  = 0.0_ti_p,&
      &one   = 1.0_ti_p,&
      &delta = 1.1_ti_p )
   precise = one
!
!     find the smallest positive floating point value;
!     minimum of 0.24x10-307 is a patch for some SGI's,
!     for Sparc cpu's under Linux, etc.
!
   if (i .eq. 1) then
!        do while (precise .ne. zero)
      do while (precise .ge. 0.24d-307)
         value = precise
         precise = precise / delta
      end do
      precise = value
!
!     find the smallest relative floating point spacing
!
   else if (i .eq. 2) then
      do while (one+precise .ne. one)
         value = precise
         precise = precise / delta
      end do
      precise = value
!
!     find the largest relative floating point spacing
!
   else if (i .eq. 3) then
      do while (one+precise .ne. precise)
         value = precise
         precise = precise * delta
      end do
      precise = value
   end if
   return
end
