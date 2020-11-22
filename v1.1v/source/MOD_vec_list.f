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
  
      contains 
  
      function  fmidpointimagevec(ncell_loc,ncell,
     &                            xk,yk,zk,
     &                            xr,yr,zr)
     &                     result(docompute)

      use domdec
      use utilvec
      implicit none
      integer, intent (in) :: ncell,ncell_loc
      real*8, contiguous, intent (in)  :: xk(:)
      real*8, contiguous, intent (in)  :: yk(:)
      real*8, contiguous, intent (in)  :: zk(:)
      real*8, contiguous, intent (in)  :: xr(:)
      real*8, contiguous, intent (in)  :: yr(:)
      real*8, contiguous, intent (in)  :: zr(:)
      logical  docompute(ncell)
      integer k
!DIR$ ATTRIBUTES ALIGN:64:: xrmid
      real*8 xrmid(ncell)
!DIR$ ATTRIBUTES ALIGN:64:: yrmid
      real*8 yrmid(ncell)
!DIR$ ATTRIBUTES ALIGN:64:: zrmid
      real*8 zrmid(ncell)
c
c     
c     definition of the middle point between i and k atoms
c
!DIR$ VECTOR ALIGNED
!DIR$ ASSUME (mod(ncell,16).eq.0)
!DIR$ SIMD
      do k = 1, ncell
         xrmid(k) = xk(k) + xr(k)/2.0d0
         yrmid(k) = yk(k) + yr(k)/2.0d0
         zrmid(k) = zk(k) + zr(k)/2.0d0
      enddo
!DIR$ ASSUME (mod(ncell,16).eq.0)
!!      call  imagevec(xrmid,ncell,1)
!!      call  imagevec(yrmid,ncell,2)
!!      call  imagevec(zrmid,ncell,3)
        call  image3dvec(xrmid,yrmid,zrmid,ncell)
!DIR$ VECTOR ALIGNED
!DIR$ ASSUME (mod(ncell,16).eq.0)
!DIR$ SIMD
      do k = 1, ncell
!!        docompute(k) = (     zrmid(k).ge.zbegproc(rank+1)
!!    &                   .and.zrmid(k).lt.zendproc(rank+1)
!!    &                   .and.yrmid(k).ge.ybegproc(rank+1)
!!    &                   .and.yrmid(k).lt.yendproc(rank+1)
!!    &                   .and.xrmid(k).ge.xbegproc(rank+1)
!!    &                   .and.xrmid(k).lt.xendproc(rank+1)
!!    &                   .and.k.le.ncell_loc)
         docompute(k) = .not.(    zrmid(k).lt.zbegproc(rank+1)
     &                        .or.zrmid(k).ge.zendproc(rank+1)
     &                        .or.yrmid(k).lt.ybegproc(rank+1)
     &                        .or.yrmid(k).ge.yendproc(rank+1)
     &                        .or.xrmid(k).lt.xbegproc(rank+1)
     &                        .or.xrmid(k).ge.xendproc(rank+1)
     &                        .or.k.gt.ncell_loc)
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
      real*8, contiguous, intent (in)  :: xk(:)
      real*8, contiguous, intent (in)  :: yk(:)
      real*8, contiguous, intent (in)  :: zk(:)
      real*8, contiguous, intent (in)  :: xr(:)
      real*8, contiguous, intent (in)  :: yr(:)
      real*8, contiguous, intent (in)  :: zr(:)
      logical, intent(out)  :: docompute(ncell)
      integer k
!DIR$ ATTRIBUTES ALIGN:64:: xrmid
      real*8 xrmid(ncell)
!DIR$ ATTRIBUTES ALIGN:64:: yrmid
      real*8 yrmid(ncell)
!DIR$ ATTRIBUTES ALIGN:64:: zrmid
      real*8 zrmid(ncell)
c
c     
c     definition of the middle point between i and k atoms
c
!DIR$ VECTOR ALIGNED
!DIR$ ASSUME (mod(ncell,16).eq.0)
!DIR$ SIMD
      do k = 1, ncell
         xrmid(k) = xk(k) + xr(k)/2.0d0
         yrmid(k) = yk(k) + yr(k)/2.0d0
         zrmid(k) = zk(k) + zr(k)/2.0d0
      enddo
!DIR$ ASSUME (mod(ncell,16).eq.0)
        call  image3dvec(xrmid,yrmid,zrmid,ncell)
!DIR$ VECTOR ALIGNED
!DIR$ ASSUME (mod(ncell,16).eq.0)
!DIR$ SIMD
      do k = 1, ncell
!!        docompute(k) = (     zrmid(k).ge.zbegproc(rank+1)
!!    &                   .and.zrmid(k).lt.zendproc(rank+1)
!!    &                   .and.yrmid(k).ge.ybegproc(rank+1)
!!    &                   .and.yrmid(k).lt.yendproc(rank+1)
!!    &                   .and.xrmid(k).ge.xbegproc(rank+1)
!!    &                   .and.xrmid(k).lt.xendproc(rank+1)
!!    &                   .and.k.le.ncell_loc)
         docompute(k) = .not.(    zrmid(k).lt.zbegproc(rank+1)
     &                        .or.zrmid(k).ge.zendproc(rank+1)
     &                        .or.yrmid(k).lt.ybegproc(rank+1)
     &                        .or.yrmid(k).ge.yendproc(rank+1)
     &                        .or.xrmid(k).lt.xbegproc(rank+1)
     &                        .or.xrmid(k).ge.xendproc(rank+1)
     &                        .or.k.gt.ncell_loc)
      enddo
c     docompute = fmidpointimagevec(ncell_loc,ncell,xk,yk,zk,xr,yr,zr)

      end subroutine midpointimagevec
      end module
