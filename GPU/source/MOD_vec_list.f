c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  vec_list  --     Stuff for vectorization purpose     ##
c     ##                                                               ##
c     ###################################################################

      module vec_list
      use sizes
      use couple
      implicit none
#include "tinker_macro.h"
  
      contains 
  
      function  fmidpointimagevec(ncell_loc,ncell,
     &                            xk,yk,zk,
     &                            xr,yr,zr)
     &                     result(docompute)

      use domdec
      use utilvec
      implicit none
      integer, intent (in) :: ncell,ncell_loc
      real(t_p), contiguous, intent (in)  :: xk(:)
      real(t_p), contiguous, intent (in)  :: yk(:)
      real(t_p), contiguous, intent (in)  :: zk(:)
      real(t_p), contiguous, intent (in)  :: xr(:)
      real(t_p), contiguous, intent (in)  :: yr(:)
      real(t_p), contiguous, intent (in)  :: zr(:)
      logical  docompute(ncell)
      integer k
!DIR$ ATTRIBUTES ALIGN:64:: xrmid
      real(t_p) xrmid(ncell)
!DIR$ ATTRIBUTES ALIGN:64:: yrmid
      real(t_p) yrmid(ncell)
!DIR$ ATTRIBUTES ALIGN:64:: zrmid
      real(t_p) zrmid(ncell)
c
c     
c     definition of the middle point between i and k atoms
c
!DIR$ VECTOR ALIGNED
!DIR$ ASSUME (mod(ncell,16).eq.0)
!DIR$ SIMD
      do k = 1, ncell
         xrmid(k) = xk(k) + xr(k)/2.0_ti_p
         yrmid(k) = yk(k) + yr(k)/2.0_ti_p
         zrmid(k) = zk(k) + zr(k)/2.0_ti_p
      enddo
!DIR$ ASSUME (mod(ncell,16).eq.0)
        call  image3dvec(xrmid,yrmid,zrmid,ncell)
!DIR$ VECTOR ALIGNED
!DIR$ ASSUME (mod(ncell,16).eq.0)
!DIR$ SIMD
      do k = 1, ncell
         docompute(k) = (     zrmid(k).ge.zbegproc(rank+1)
     &                   .and.zrmid(k).lt.zendproc(rank+1)
     &                   .and.yrmid(k).ge.ybegproc(rank+1)
     &                   .and.yrmid(k).lt.yendproc(rank+1)
     &                   .and.xrmid(k).ge.xbegproc(rank+1)
     &                   .and.xrmid(k).lt.xendproc(rank+1)
     &                   .and.k.le.ncell_loc)
      enddo
      end function  fmidpointimagevec
c
c     compute the selection mask with the midpoint method
c     
      subroutine midpointimagevec(ncell_loc,ncell,
     &                            xk,yk,zk,
     &                            xr,yr,zr,
     &                            docompute)
      use domdec
      use utilvec
      implicit none
      integer, intent (in) :: ncell,ncell_loc
      real(t_p), contiguous, intent (in)  :: xk(:)
      real(t_p), contiguous, intent (in)  :: yk(:)
      real(t_p), contiguous, intent (in)  :: zk(:)
      real(t_p), contiguous, intent (in)  :: xr(:)
      real(t_p), contiguous, intent (in)  :: yr(:)
      real(t_p), contiguous, intent (in)  :: zr(:)
      logical, intent(out)  :: docompute(ncell)

      docompute = fmidpointimagevec(ncell_loc,ncell,xk,yk,zk,xr,yr,zr)

      end subroutine midpointimagevec
      end module
