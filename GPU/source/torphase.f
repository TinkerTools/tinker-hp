c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine torphase  --  torsional amplitude and phase  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torphase" sets the n-fold amplitude and phase values
c     for each torsion via sorting of the input parameters
c
c
#include "tinker_precision.h"
      subroutine torphase (ft,vt,st)
      implicit none
      integer,parameter::ti_p=t_p
      integer i,k
      integer ft(6)
      real(t_p) vt(6),st(6)
      real(t_p) ampli(6),phase(6)
c
c
c     copy the input fold, amplitude and phase angles
c
      do i = 1, 6
         ampli(i) = vt(i)
         phase(i) = st(i)
         vt(i) = 0.0_ti_p
         st(i) = 0.0_ti_p
      end do
c
c     shift the phase angles into the standard range
c
      do i = 1, 6
         do while (phase(i) .lt. -180.0_ti_p)
            phase(i) = phase(i) + 360.0_ti_p
         end do
         do while (phase(i) .gt. 180.0_ti_p)
            phase(i) = phase(i) - 360.0_ti_p
         end do
      end do
c
c     convert input torsional parameters to storage format
c
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
