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
!> @brief 
!> "precise" finds a machine precision value as selected by
!> the input argument: (1) the smallest positive floating
!> point value, (2) the smallest relative floating point
!> spacing, (3) the largest relative floating point spacing
!> @param no params
function precise (i)
   implicit none
   integer i
   real*8 precise,value
   real*8 zero,one,delta
!
!
!     set values for zero, one and multiplicative factor
!
   zero = 0.0d0
   one = 1.0d0
   delta = 1.1d0
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
