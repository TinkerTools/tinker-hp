!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine torphase  --  torsional amplitude and phase  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "torphase" sets the n-fold amplitude and phase values
!     for each torsion via sorting of the input parameters
!
!
#include "tinker_precision.h"
subroutine torphase (ft,vt,st)
   implicit none
   integer,parameter::ti_p=t_p
   integer i,k
   integer ft(6)
   real(t_p) vt(6),st(6)
   real(t_p) ampli(6),phase(6)
!
!
!     copy the input fold, amplitude and phase angles
!
   do i = 1, 6
      ampli(i) = vt(i)
      phase(i) = st(i)
      vt(i) = 0.0_ti_p
      st(i) = 0.0_ti_p
   end do
!
!     shift the phase angles into the standard range
!
   do i = 1, 6
      do while (phase(i) .lt. -180.0_ti_p)
         phase(i) = phase(i) + 360.0_ti_p
      end do
      do while (phase(i) .gt. 180.0_ti_p)
         phase(i) = phase(i) - 360.0_ti_p
      end do
   end do
!
!     convert input torsional parameters to storage format
!
   do i = 1, 6
      k = ft(i)
      if (k .eq. 0) then
         goto 10
      else if (k .le. 6) then
         vt(k) = ampli(i)
         st(k) = phase(i)
      end if
   end do
10 continue
   return
end
