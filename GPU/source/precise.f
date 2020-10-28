c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #########################################################
c     ##                                                     ##
c     ##  function precise  --  determine machine precision  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "precise" finds a machine precision value as selected by
c     the input argument: (1) the smallest positive floating
c     point value, (2) the smallest relative floating point
c     spacing, (3) the largest relative floating point spacing
c
c
#include "tinker_precision.h"
      function precise (i)
      implicit none
      integer,parameter::ti_p=t_p
      integer i
      real(t_p) precise,value
      real(t_p) zero,one,delta
c
c
c     set values for zero, one and multiplicative factor
c
      parameter (
     &  zero  = 0.0_ti_p,
     &  one   = 1.0_ti_p,
     &  delta = 1.1_ti_p )
      precise = one
c
c     find the smallest positive floating point value;
c     minimum of 0.24x10-307 is a patch for some SGI's,
c     for Sparc cpu's under Linux, etc.
c
      if (i .eq. 1) then
c        do while (precise .ne. zero)
         do while (precise .ge. 0.24d-307)
            value = precise
            precise = precise / delta
         end do
         precise = value
c
c     find the smallest relative floating point spacing
c
      else if (i .eq. 2) then
         do while (one+precise .ne. one)
            value = precise
            precise = precise / delta
         end do
         precise = value
c
c     find the largest relative floating point spacing
c
      else if (i .eq. 3) then
         do while (one+precise .ne. precise)
           value = precise
           precise = precise * delta
         end do
         precise = value
      end if
      return
      end
