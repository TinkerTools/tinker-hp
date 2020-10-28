c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torquevec  --  convert single site torque to force  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque2vec" takes the torque values on multiple sites defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     npolelocnl and npolelocnlloop came from mpole module
c
c     INPUT  : ivec is indices of sites 
c              trqvec are trqs of sites 
c              devecx,devecy and devecz are deps coming from epolar
c              or empole
c
c     OUTPUT : dvvecx, dwvecx and duvecx, etc... are forces 
c              devecx,devecy and devecz are deps sent back to epolar
c              or empole
#include "tinker_precision.h"
      subroutine torquevec2 (ivec,trqvecx,trqvecy,trqvecz,
     &                      dvvecx,dvvecy,dvvecz,
     &                      dwvecx,dwvecy,dwvecz,
     &                      duvecx,duvecy,duvecz,
     &                      devecx,devecy,devecz)
      use atoms
      use deriv
      use domdec
      use mpole
      use sizes
      use mpi
      use timestat
      use tinheader ,only: ti_p,re_p
      implicit none
      integer i,j,k,kk
      integer npolelocnl1
      integer ivec(npolelocnlloop)
      real(t_p) sqrt3over2
      real(t_p) half
      real(t_p) trqvecx(npolelocnlloop)
      real(t_p) trqvecy(npolelocnlloop)
      real(t_p) trqvecz(npolelocnlloop)
      real(t_p) dvvecx(npolelocnlloop)
      real(t_p) dvvecy(npolelocnlloop)
      real(t_p) dvvecz(npolelocnlloop)
      real(t_p) dwvecx(npolelocnlloop)
      real(t_p) dwvecy(npolelocnlloop)
      real(t_p) dwvecz(npolelocnlloop)
      real(t_p) duvecx(npolelocnlloop)
      real(t_p) duvecy(npolelocnlloop)
      real(t_p) duvecz(npolelocnlloop)
      real(t_p) devecx(*)
      real(t_p) devecy(*)
      real(t_p) devecz(*)
!DIR$ ATTRIBUTES ALIGN:64::axetypvec
      integer axetypvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec1
      integer ivec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::iavec,ibvec
      integer iavec(npolelocnlloop),ibvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::icvec,idvec
      integer icvec(npolelocnlloop),idvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ialocvec,iblocvec
      integer ialocvec(npolelocnlloop),iblocvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::iclocvec,idlocvec
      integer iclocvec(npolelocnlloop),idlocvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: delvecx,delvecy
      real(t_p) delvecx(npolelocnlloop)
      real(t_p) delvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: delvecz
      real(t_p) delvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::epsvecx,epsvecy,epsvecz
      real(t_p) epsvecx(npolelocnlloop)
      real(t_p) epsvecy(npolelocnlloop)
      real(t_p) epsvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invdelsizvec,epssizvec
      real(t_p) invdelsizvec(npolelocnlloop)
      real(t_p) epssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::pvecx,pvecy,pvecz,invpsizvec
      real(t_p) pvecx(npolelocnlloop)
      real(t_p) pvecy(npolelocnlloop)
      real(t_p) pvecz(npolelocnlloop)
      real(t_p) invpsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rvecx,rvecy,rvecz
      real(t_p) rvecx(npolelocnlloop)
      real(t_p) rvecy(npolelocnlloop)
      real(t_p) rvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invrsizvec
      real(t_p) invrsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::svecx,svecy,svecz
      real(t_p) svecx(npolelocnlloop)
      real(t_p) svecy(npolelocnlloop)
      real(t_p) svecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invssizvec
      real(t_p) invssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::invusizvec
      real(t_p) invusizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uvecx,uvecy,uvecz
      real(t_p) uvecx(npolelocnlloop)
      real(t_p) uvecy(npolelocnlloop)
      real(t_p) uvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::invvsizvec
      real(t_p) invvsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vvecx,vvecy,vvecz
      real(t_p) vvecx(npolelocnlloop)
      real(t_p) vvecy(npolelocnlloop)
      real(t_p) vvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::invwsizvec
      real(t_p) invwsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wvecx,wvecy,wvecz
      real(t_p) wvecx(npolelocnlloop)
      real(t_p) wvecy(npolelocnlloop)
      real(t_p) wvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uvvecx,uvvecy,uvvecz
      real(t_p) uvvecx(npolelocnlloop)
      real(t_p) uvvecy(npolelocnlloop)
      real(t_p) uvvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invuvsizvec
      real(t_p) invuvsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uwvecx,uwvecy,uwvecz
      real(t_p) uwvecx(npolelocnlloop)
      real(t_p) uwvecy(npolelocnlloop)
      real(t_p) uwvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invuwsizvec
      real(t_p) invuwsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vwvecx,vwvecy,vwvecz
      real(t_p) vwvecx(npolelocnlloop)
      real(t_p) vwvecy(npolelocnlloop)
      real(t_p) vwvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invvwsizvec
      real(t_p) invvwsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::usvecx,usvecy,usvecz
      real(t_p) usvecx(npolelocnlloop)
      real(t_p) usvecy(npolelocnlloop)
      real(t_p) usvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invussizvec
      real(t_p) invussizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vsvecx,vsvecy,vsvecz
      real(t_p) vsvecx(npolelocnlloop)
      real(t_p) vsvecy(npolelocnlloop)
      real(t_p) vsvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invvssizvec
      real(t_p) invvssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wsvecx,wsvecy,wsvecz
      real(t_p) wsvecx(npolelocnlloop)
      real(t_p) wsvecy(npolelocnlloop)
      real(t_p) wsvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invwssizvec
      real(t_p) invwssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::urvecx,urvecy,urvecz
      real(t_p) urvecx(npolelocnlloop)
      real(t_p) urvecy(npolelocnlloop)
      real(t_p) urvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invursizvec
      real(t_p) invursizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::t1vecx,t1vecy,t1vecz
      real(t_p) t1vecx(npolelocnlloop)
      real(t_p) t1vecy(npolelocnlloop)
      real(t_p) t1vecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invt1sizvec
      real(t_p) invt1sizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::t2vecx,t2vecy,t2vecz
      real(t_p) t2vecx(npolelocnlloop)
      real(t_p) t2vecy(npolelocnlloop)
      real(t_p) t2vecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invt2sizvec
      real(t_p) invt2sizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rucosvec,invrusinvec
      real(t_p) rucosvec(npolelocnlloop)
      real(t_p) invrusinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rvcosvec,invrvsinvec
      real(t_p) rvcosvec(npolelocnlloop)
      real(t_p) invrvsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rwcosvec,invrwsinvec
      real(t_p) rwcosvec(npolelocnlloop)
      real(t_p) invrwsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uvcosvec,invuvsinvec
      real(t_p) uvcosvec(npolelocnlloop)
      real(t_p) invuvsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uwcosvec,uwsinvec
      real(t_p) uwcosvec(npolelocnlloop)
      real(t_p) uwsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vwcosvec,vwsinvec
      real(t_p) vwcosvec(npolelocnlloop)
      real(t_p) vwsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::upcosvec
      real(t_p) upcosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::urcosvec,invursinvec
      real(t_p) urcosvec(npolelocnlloop)
      real(t_p) invursinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uscosvec
      real(t_p) uscosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vpcosvec
      real(t_p) vpcosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vscosvec,vssinvec
      real(t_p) vscosvec(npolelocnlloop)
      real(t_p) vssinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wpcosvec
      real(t_p) wpcosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wscosvec,wssinvec
      real(t_p) wscosvec(npolelocnlloop)
      real(t_p) wssinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ut1cosvec,ut1sinvec
      real(t_p) ut1cosvec(npolelocnlloop)
      real(t_p) ut1sinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ut2cosvec,ut2sinvec,invutsinvec
      real(t_p) ut2cosvec(npolelocnlloop)
      real(t_p) ut2sinvec(npolelocnlloop)
      real(t_p) invutsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphiduvec,dphidvvec
      real(t_p) dphiduvec(npolelocnlloop)
      real(t_p) dphidvvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphidwvec,dphidrvec
      real(t_p) dphidwvec(npolelocnlloop)
      real(t_p) dphidrvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphiddelvec
      real(t_p) dphiddelvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphidsvec
      real(t_p) dphidsvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xiavec,yiavec,ziavec
      real(t_p) xiavec(npolelocnlloop)
      real(t_p) yiavec(npolelocnlloop)
      real(t_p) ziavec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xibvec,yibvec,zibvec
      real(t_p) xibvec(npolelocnlloop)
      real(t_p) yibvec(npolelocnlloop)
      real(t_p) zibvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xicvec,yicvec,zicvec
      real(t_p) xicvec(npolelocnlloop)
      real(t_p) yicvec(npolelocnlloop)
      real(t_p) zicvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xidvec,yidvec,zidvec
      real(t_p) xidvec(npolelocnlloop)
      real(t_p) yidvec(npolelocnlloop)
      real(t_p) zidvec(npolelocnlloop)
c     3-Fold
!DIR$ ATTRIBUTES ALIGN:64:: c3foldvec
      real(t_p) c3foldvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecx
      real(t_p) du3foldvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecy
      real(t_p) du3foldvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecz
      real(t_p) du3foldvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecx
      real(t_p) dv3foldvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecy
      real(t_p) dv3foldvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecz
      real(t_p) dv3foldvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecx
      real(t_p) dw3foldvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecy
      real(t_p) dw3foldvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecz
      real(t_p) dw3foldvecz(npolelocnlloop)
c     Bisector
!DIR$ ATTRIBUTES ALIGN:64:: cbisectorvec
      real(t_p) cbisectorvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecx
      real(t_p) dubisectorvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecy
      real(t_p) dubisectorvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecz
      real(t_p) dubisectorvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecx
      real(t_p) dvbisectorvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecy
      real(t_p) dvbisectorvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecz
      real(t_p) dvbisectorvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecx
      real(t_p) dwbisectorvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecy
      real(t_p) dwbisectorvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecz
      real(t_p) dwbisectorvecz(npolelocnlloop)
c     Z-Bisect
!DIR$ ATTRIBUTES ALIGN:64:: czbisectvec
      real(t_p) czbisectvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecx
      real(t_p) duzbisectvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecy
      real(t_p) duzbisectvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecz
      real(t_p) duzbisectvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecx
      real(t_p) dvzbisectvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecy
      real(t_p) dvzbisectvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecz
      real(t_p) dvzbisectvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecx
      real(t_p) dwzbisectvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecy
      real(t_p) dwzbisectvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecz
      real(t_p) dwzbisectvecz(npolelocnlloop)
C     Z-Only
!DIR$ ATTRIBUTES ALIGN:64:: czonlyvec
      real(t_p) czonlyvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecx
      real(t_p) duzonlyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecy
      real(t_p) duzonlyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecz
      real(t_p) duzonlyvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecx
      real(t_p) dvzonlyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecy
      real(t_p) dvzonlyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecz
      real(t_p) dvzonlyvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecx
      real(t_p) dwzonlyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecy
      real(t_p) dwzonlyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecz
      real(t_p) dwzonlyvecz(npolelocnlloop)
C     Z-then-X
!DIR$ ATTRIBUTES ALIGN:64:: czthenxvec
      real(t_p) czthenxvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecx
      real(t_p) duzthenxvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecy
      real(t_p) duzthenxvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecz
      real(t_p) duzthenxvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecx
      real(t_p) dvzthenxvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecy
      real(t_p) dvzthenxvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecz
      real(t_p) dvzthenxvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecx
      real(t_p) dwzthenxvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecy
      real(t_p) dwzthenxvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecz
      real(t_p) dwzthenxvecz(npolelocnlloop)
      character*8 paxe(npolelocnlloop)

c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'torquevec2'

      call timer_enter( timer_torque )
      sqrt3over2 = sqrt(3.0_ti_p)/2.0_ti_p
      half       = 0.5_ti_p
c
c     zero out force components on local frame-defining atoms
c
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         iavec(k)     = 1
         ibvec(k)     = 1
         icvec(k)     = 1
         idvec(k)     = 1
         ialocvec(k)  = 1
         iblocvec(k)  = 1
         iclocvec(k)  = 1
         idlocvec(k)  = 1
         ivec1(k)     = 1
         axetypvec(k) = 0
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         c3foldvec(k)   = 0.0_ti_p
         du3foldvecx(k) = 0.0_ti_p
         du3foldvecy(k) = 0.0_ti_p
         du3foldvecz(k) = 0.0_ti_p
         dv3foldvecx(k) = 0.0_ti_p
         dv3foldvecy(k) = 0.0_ti_p
         dv3foldvecz(k) = 0.0_ti_p
         dw3foldvecx(k) = 0.0_ti_p
         dw3foldvecy(k) = 0.0_ti_p
         dw3foldvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         cbisectorvec(k)   = 0.0_ti_p
         dubisectorvecx(k) = 0.0_ti_p
         dubisectorvecy(k) = 0.0_ti_p
         dubisectorvecz(k) = 0.0_ti_p
         dvbisectorvecx(k) = 0.0_ti_p
         dvbisectorvecy(k) = 0.0_ti_p
         dvbisectorvecz(k) = 0.0_ti_p
         dwbisectorvecx(k) = 0.0_ti_p
         dwbisectorvecy(k) = 0.0_ti_p
         dwbisectorvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         czbisectvec(k)   = 0.0_ti_p
         duzbisectvecx(k) = 0.0_ti_p
         duzbisectvecy(k) = 0.0_ti_p
         duzbisectvecz(k) = 0.0_ti_p
         dvzbisectvecx(k) = 0.0_ti_p
         dvzbisectvecy(k) = 0.0_ti_p
         dvzbisectvecz(k) = 0.0_ti_p
         dwzbisectvecx(k) = 0.0_ti_p
         dwzbisectvecy(k) = 0.0_ti_p
         dwzbisectvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         czonlyvec(k)   = 0.0_ti_p
         duzonlyvecx(k) = 0.0_ti_p
         duzonlyvecy(k) = 0.0_ti_p
         duzonlyvecz(k) = 0.0_ti_p
         dvzonlyvecx(k) = 0.0_ti_p
         dvzonlyvecy(k) = 0.0_ti_p
         dvzonlyvecz(k) = 0.0_ti_p
         dwzonlyvecx(k) = 0.0_ti_p
         dwzonlyvecy(k) = 0.0_ti_p
         dwzonlyvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         czthenxvec(k)   = 0.0_ti_p
         duzthenxvecx(k) = 0.0_ti_p
         duzthenxvecy(k) = 0.0_ti_p
         duzthenxvecz(k) = 0.0_ti_p
         dvzthenxvecx(k) = 0.0_ti_p
         dvzthenxvecy(k) = 0.0_ti_p
         dvzthenxvecz(k) = 0.0_ti_p
         dwzthenxvecx(k) = 0.0_ti_p
         dwzthenxvecy(k) = 0.0_ti_p
         dwzthenxvecz(k) = 0.0_ti_p
      enddo

c
c     get the local frame type and the frame-defining atoms
c
      do k = 1, npolelocnl
         select case (polaxe(ivec(k)))
                case ('3-Fold'  ) ;  axetypvec(k) = 1
                                     c3foldvec(k)= 1.0_ti_p
                case ('Bisector') ;  axetypvec(k) = 2
                                     cbisectorvec(k)= 1.0_ti_p
                case ('Z-Bisect') ;  axetypvec(k) = 3
                                     czbisectvec(k)= 1.0_ti_p
                case ('Z-Only'  ) ;  axetypvec(k) = 4
                                     czonlyvec(k)= 1.0_ti_p
                case ('Z-then-X') ;  axetypvec(k) = 5
                                     czthenxvec(k)= 1.0_ti_p
                case ('None'    ) ;  axetypvec(k) = 0
                case default      ;  axetypvec(k) = 0
         endselect
      enddo

      kk = 0
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
      do k = 1,npolelocnlloop
         if (axetypvec(k).ne.0) then
            kk = kk +1
            ivec1 (kk) = ivec (k)
         endif
      enddo
      npolelocnl1 = kk
      if (npolelocnl1.eq.0) return ! All axetyp are 0, so return
    
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         iavec(k) = zaxis(ivec1(k))
         if(iavec(k).gt.0) then
            ialocvec(k) = loc(iavec(k))
         else
            iavec(k) = 1
         endif
         if(ialocvec(k).eq.0)  ialocvec(k) = 1
!!        ibvec(k) = ipole(ivec1(k))
!!        icvec(k) = xaxis(ivec1(k))
!!        idvec(k) = yaxis(ivec1(k))
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         ibvec(k) = ipole(ivec1(k))
         if(ibvec(k).gt.0) then
            iblocvec(k) = loc(ibvec(k))
         else
            ibvec(k) = 1
         endif
         if(iblocvec(k).eq.0)  iblocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         icvec(k) = xaxis(ivec1(k))
         if(icvec(k).gt.0) then
            iclocvec(k) = loc(icvec(k))
         else
            icvec(k) = 1
         endif
         if(iclocvec(k).eq.0)  iclocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         idvec(k) = yaxis(ivec1(k))
         if(idvec(k).gt.0) then
            idlocvec(k) = loc(idvec(k))
         else
            idvec(k) = 1
         endif
         if(idlocvec(k).eq.0)  idlocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         xiavec(k)=x(iavec(k))
         yiavec(k)=y(iavec(k))
         ziavec(k)=z(iavec(k))
         xibvec(k)=x(ibvec(k))
         yibvec(k)=y(ibvec(k))
         zibvec(k)=z(ibvec(k))
         xicvec(k)=x(icvec(k))
         yicvec(k)=y(icvec(k))
         zicvec(k)=z(icvec(k))
         xidvec(k)=x(idvec(k))
         yidvec(k)=y(idvec(k))
         zidvec(k)=z(idvec(k))
      enddo
c
c     construct the three rotation axes for the local frame
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvecx(k)  = xiavec(k)- xibvec(k)
         uvecy(k)  = yiavec(k)- yibvec(k)
         uvecz(k)  = ziavec(k)- zibvec(k)
         invusizvec(k) = (  uvecx(k)**2 + uvecy(k)**2 + uvecz(k)**2
     &                   ) ** ( - half )
         uvecx(k) = uvecx(k) * invusizvec(k)
         uvecy(k) = uvecy(k) * invusizvec(k)
         uvecz(k) = uvecz(k) * invusizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vvecx(k)  = xicvec(k) - xibvec(k)
         vvecy(k)  = yicvec(k) - yibvec(k)
         vvecz(k)  = zicvec(k) - zibvec(k)
         invvsizvec(k) = (  vvecx(k)**2 + vvecy(k)**2 + vvecz(k)**2
     &                   ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         if (axetypvec(k).eq.4) then ! 'Z-Only' 
            if (abs(uvecx(k)).gt.sqrt3over2) then
               vvecx(k)  = 0.0_ti_p
            else
               vvecx(k)  = 1.0_ti_p
            endif
            vvecy(k) = 1.0_ti_p - vvecx(k)
            vvecz(k)  = 0.0_ti_p
            invvsizvec(k) = 1.0_ti_p
         endif
      enddo
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         if (    axetypvec(k).eq.1  ! '3-Fold'
     &       .or.axetypvec(k).eq.3) ! 'Z-Bisect'
     &   then
            wvecx(k) = xidvec(k) - xibvec(k)
            wvecy(k) = yidvec(k) - yibvec(k)
            wvecz(k) = zidvec(k) - zibvec(k)
         else
           wvecx(k) =  uvecy(k) * vvecz(k) - uvecz(k) * vvecy(k)
           wvecy(k) =  uvecz(k) * vvecx(k) - uvecx(k) * vvecz(k)
           wvecz(k) =  uvecx(k) * vvecy(k) - uvecy(k) * vvecx(k)
         endif
         invwsizvec(k) = (  wvecx(k)**2 + wvecy(k)**2
     &                    + wvecz(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vvecx(k) = vvecx(k) * invvsizvec(k)
         vvecy(k) = vvecy(k) * invvsizvec(k)
         vvecz(k) = vvecz(k) * invvsizvec(k)
         wvecx(k) = wvecx(k) * invwsizvec(k)
         wvecy(k) = wvecy(k) * invwsizvec(k)
         wvecz(k) = wvecz(k) * invwsizvec(k)
      enddo
c
c     build some additional axes needed for the Z-Bisect method
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)  = vvecx(k) + wvecx(k)
         rvecy(k)  = vvecy(k) + wvecy(k)
         rvecz(k)  = vvecz(k) + wvecz(k)
         svecx(k) =  uvecy(k) * rvecz(k) - uvecz(k) * rvecy(k)
         svecy(k) =  uvecz(k) * rvecx(k) - uvecx(k) * rvecz(k)
         svecz(k) =  uvecx(k) * rvecy(k) - uvecy(k) * rvecx(k)
         invrsizvec(k) = (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                   ) ** ( - half )
         invssizvec(k) = (  svecx(k)**2 + svecy(k)**2 + svecz(k)**2
     &                   ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k) = rvecx(k) * invrsizvec(k)
         svecx(k) = svecx(k) * invssizvec(k)
         rvecy(k) = rvecy(k) * invrsizvec(k)
         svecy(k) = svecy(k) * invssizvec(k)
         rvecz(k) = rvecz(k) * invrsizvec(k)
         svecz(k) = svecz(k) * invssizvec(k)
      enddo
c     find the perpendicular and angle for each pair of axes
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvvecx(k)  =  vvecy(k) * uvecz(k) - vvecz(k) * uvecy(k)
         uvvecy(k)  =  vvecz(k) * uvecx(k) - vvecx(k) * uvecz(k)
         uvvecz(k)  =  vvecx(k) * uvecy(k) - vvecy(k) * uvecx(k)
         invuvsizvec(k) = (  uvvecx(k)**2 + uvvecy(k)**2 + uvvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uwvecx(k)  =  wvecy(k) * uvecz(k) - wvecz(k) * uvecy(k)
         uwvecy(k)  =  wvecz(k) * uvecx(k) - wvecx(k) * uvecz(k)
         uwvecz(k)  =  wvecx(k) * uvecy(k) - wvecy(k) * uvecx(k)
         invuwsizvec(k) = (  uwvecx(k)**2 + uwvecy(k)**2 + uwvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vwvecx(k)  =  wvecy(k) * vvecz(k) - wvecz(k) * vvecy(k)
         vwvecy(k)  =  wvecz(k) * vvecx(k) - wvecx(k) * vvecz(k)
         vwvecz(k)  =  wvecx(k) * vvecy(k) - wvecy(k) * vvecx(k)
         invvwsizvec(k) = (  vwvecx(k)**2 + vwvecy(k)**2 + vwvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvvecx(k) = uvvecx(k) * invuvsizvec(k)
         uvvecy(k) = uvvecy(k) * invuvsizvec(k)
         uvvecz(k) = uvvecz(k) * invuvsizvec(k)
         uwvecx(k) = uwvecx(k) * invuwsizvec(k)
         uwvecy(k) = uwvecy(k) * invuwsizvec(k)
         uwvecz(k) = uwvecz(k) * invuwsizvec(k)
         vwvecx(k) = vwvecx(k) * invvwsizvec(k)
         vwvecy(k) = vwvecy(k) * invvwsizvec(k)
         vwvecz(k) = vwvecz(k) * invvwsizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         urvecx(k)  =  rvecy(k) * uvecz(k) - rvecz(k) * uvecy(k)
         urvecy(k)  =  rvecz(k) * uvecx(k) - rvecx(k) * uvecz(k)
         urvecz(k)  =  rvecx(k) * uvecy(k) - rvecy(k) * uvecx(k)
         invursizvec(k) = (  urvecx(k)**2 + urvecy(k)**2 + urvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         usvecx(k)  =  svecy(k) * uvecz(k) - svecz(k) * uvecy(k)
         usvecy(k)  =  svecz(k) * uvecx(k) - svecx(k) * uvecz(k)
         usvecz(k)  =  svecx(k) * uvecy(k) - svecy(k) * uvecx(k)
         invussizvec(k) = (  usvecx(k)**2 + usvecy(k)**2 + usvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vsvecx(k)  =  svecy(k) * vvecz(k) - svecz(k) * vvecy(k)
         vsvecy(k)  =  svecz(k) * vvecx(k) - svecx(k) * vvecz(k)
         vsvecz(k)  =  svecx(k) * vvecy(k) - svecy(k) * vvecx(k)
         invvssizvec(k) = (  vsvecx(k)**2 + vsvecy(k)**2 + vsvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         wsvecx(k)  =  svecy(k) * wvecz(k) - svecz(k) * wvecy(k)
         wsvecy(k)  =  svecz(k) * wvecx(k) - svecx(k) * wvecz(k)
         wsvecz(k)  =  svecx(k) * wvecy(k) - svecy(k) * wvecx(k)
         invwssizvec(k) = (  wsvecx(k)**2 + wsvecy(k)**2 + wsvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         urvecx(k)  = urvecx(k) * invursizvec(k)
         usvecx(k)  = usvecx(k) * invussizvec(k)
         vsvecx(k)  = vsvecx(k) * invvssizvec(k)
         wsvecx(k)  = wsvecx(k) * invwssizvec(k)
         urvecy(k)  = urvecy(k) * invursizvec(k)
         usvecy(k)  = usvecy(k) * invussizvec(k)
         vsvecy(k)  = vsvecy(k) * invvssizvec(k)
         wsvecy(k)  = wsvecy(k) * invwssizvec(k)
         urvecz(k)  = urvecz(k) * invursizvec(k)
         usvecz(k)  = usvecz(k) * invussizvec(k)
         vsvecz(k)  = vsvecz(k) * invvssizvec(k)
         wsvecz(k)  = wsvecz(k) * invwssizvec(k)
      enddo
c
c     get sine and cosine of angles between the rotation axes
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvcosvec(k) =  uvecx(k) * vvecx(k) + uvecy(k) * vvecy(k)
     &                + uvecz(k) * vvecz(k)
         uwcosvec(k) =  uvecx(k) * wvecx(k) + uvecy(k) * wvecy(k)
     &                + uvecz(k) * wvecz(k)
         vwcosvec(k) =  vvecx(k) * wvecx(k) + vvecy(k) * wvecy(k)
     &                + vvecz(k) * wvecz(k)
         invuvsinvec(k) = (1.0_ti_p - uvcosvec(k)**2) ** ( - half )
      enddo
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         urcosvec(k) =  uvecx(k) * rvecx(k) + uvecy(k) * rvecy(k)
     &                + uvecz(k) * rvecz(k)
         uscosvec(k) =  uvecx(k) * svecx(k) + uvecy(k) * svecy(k)
     &                + uvecz(k) * svecz(k)
         vscosvec(k) =  vvecx(k) * svecx(k) + vvecy(k) * svecy(k)
     &                + vvecz(k) * svecz(k)
         wscosvec(k) =  wvecx(k) * svecx(k) + wvecy(k) * svecy(k)
     &                + wvecz(k) * svecz(k)
      enddo
c
c     compute the projection of v and w onto the ru-plane
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         t1vecx(k) =  vvecx(k) - svecx(k) * vscosvec(k)
         t1vecy(k) =  vvecy(k) - svecy(k) * vscosvec(k)
         t1vecz(k) =  vvecz(k) - svecz(k) * vscosvec(k)
         t2vecx(k) =  wvecx(k) - svecx(k) * wscosvec(k)
         t2vecy(k) =  wvecy(k) - svecy(k) * wscosvec(k)
         t2vecz(k) =  wvecz(k) - svecz(k) * wscosvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         invt1sizvec(k) = (  t1vecx(k)**2 + t1vecy(k)**2 + t1vecz(k)**2
     &                    ) ** ( - half )
         invt2sizvec(k) = (  t2vecx(k)**2 + t2vecy(k)**2 + t2vecz(k)**2
     &                    ) ** ( - half )
         t1vecx(k) = t1vecx(k) * invt1sizvec(k)
         t1vecy(k) = t1vecy(k) * invt1sizvec(k)
         t1vecz(k) = t1vecz(k) * invt1sizvec(k)
         t2vecx(k) = t2vecx(k) * invt2sizvec(k)
         t2vecy(k) = t2vecy(k) * invt2sizvec(k)
         t2vecz(k) = t2vecz(k) * invt2sizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         ut1cosvec(k) =  uvecx(k) * t1vecx(k) + uvecy(k) * t1vecy(k)
     &                 + uvecz(k) * t1vecz(k)
         ut2cosvec(k) =  uvecx(k) * t2vecx(k) + uvecy(k) * t2vecy(k)
     &                 + uvecz(k) * t2vecz(k)
         invutsinvec(k) = (  (1.0_ti_p - ut1cosvec(k)**2) ** half
     &                     + (1.0_ti_p - ut2cosvec(k)**2) ** half
     &                    ) **(-1.0_ti_p)
      enddo
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphiduvec(k) = - trqvecx(k) * uvecx(k) - trqvecy(k) * uvecy(k)
     &                  - trqvecz(k) * uvecz(k)
         dphidvvec(k) = - trqvecx(k) * vvecx(k) - trqvecy(k) * vvecy(k)
     &                  - trqvecz(k) * vvecz(k)
         dphidwvec(k) = - trqvecx(k) * wvecx(k) - trqvecy(k) * wvecy(k)
     &                  - trqvecz(k) * wvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphidrvec(k) = - trqvecx(k) * rvecx(k) - trqvecy(k) * rvecy(k)
     &                  - trqvecz(k) * rvecz(k)
         dphidsvec(k) = - trqvecx(k) * svecx(k) - trqvecy(k) * svecy(k)
     &                  - trqvecz(k) * svecz(k)
      enddo
c
c     force distribution for the Z-Only local coordinate method
c
          
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duzonlyvecx(k) =   uvvecx(k)     * dphidvvec(k) * invusizvec(k)
     &                    * invuvsinvec(k)
     &                   +  uwvecx(k)     * dphidwvec(k) * invusizvec(k)
         duzonlyvecy(k) =   uvvecy(k)     * dphidvvec(k) * invusizvec(k)
     &                    * invuvsinvec(k)
     &                   +  uwvecy(k)     * dphidwvec(k) * invusizvec(k)
         duzonlyvecz(k) =   uvvecz(k)     * dphidvvec(k) * invusizvec(k)
     &                    * invuvsinvec(k)
     &                   +  uwvecz(k)     * dphidwvec(k) * invusizvec(k)

         duzonlyvecx(k) = duzonlyvecx(k)  * czonlyvec(k)
         duzonlyvecy(k) = duzonlyvecy(k)  * czonlyvec(k)
         duzonlyvecz(k) = duzonlyvecz(k)  * czonlyvec(k)
      enddo
c
c     force distribution for the Z-then-X local coordinate method
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duzthenxvecx(k) =   uvvecx(k) * dphidvvec(k)    * invusizvec(k)
     &                                 * invuvsinvec(k)
     &                    +  uwvecx(k) * dphidwvec(k)    * invusizvec(k)
         duzthenxvecy(k) =   uvvecy(k) * dphidvvec(k)    * invusizvec(k)
     &                                 * invuvsinvec(k)
     &                    +  uwvecy(k) * dphidwvec(k)    * invusizvec(k)
         duzthenxvecz(k) =   uvvecz(k) * dphidvvec(k)    * invusizvec(k)
     &                                 * invuvsinvec(k)
     &                    +  uwvecz(k) * dphidwvec(k)    * invusizvec(k)

         dvzthenxvecx(k) = -  uvvecx(k) * dphiduvec(k)   * invvsizvec(k)
     &                                  * invuvsinvec(k)
         dvzthenxvecy(k) = -  uvvecy(k) * dphiduvec(k)   * invvsizvec(k)
     &                                  * invuvsinvec(k)
         dvzthenxvecz(k) = -  uvvecz(k) * dphiduvec(k)   * invvsizvec(k)
     &                                  * invuvsinvec(k)

         duzthenxvecx(k) = duzthenxvecx(k) * czthenxvec(k)
         duzthenxvecy(k) = duzthenxvecy(k) * czthenxvec(k)
         duzthenxvecz(k) = duzthenxvecz(k) * czthenxvec(k)
         dvzthenxvecx(k) = dvzthenxvecx(k) * czthenxvec(k)
         dvzthenxvecy(k) = dvzthenxvecy(k) * czthenxvec(k)
         dvzthenxvecz(k) = dvzthenxvecz(k) * czthenxvec(k)
      enddo

c
c     force distribution for the Bisector local coordinate method
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dubisectorvecx(k) =           uvvecx(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * uwvecx(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dubisectorvecy(k) =           uvvecy(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * uwvecy(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dubisectorvecz(k) =           uvvecz(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * uwvecz(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dvbisectorvecx(k) = -         uvvecx(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * vwvecx(k)     * dphidwvec(k)
     &                               * invvsizvec(k)
         dvbisectorvecy(k) = -         uvvecy(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * vwvecy(k)     * dphidwvec(k)
     &                               * invvsizvec(k)
         dvbisectorvecz(k) = -         uvvecz(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * vwvecz(k)     * dphidwvec(k)
     &                               * invvsizvec(k)

         dubisectorvecx(k) = dubisectorvecx(k) * cbisectorvec(k)
         dubisectorvecy(k) = dubisectorvecy(k) * cbisectorvec(k)
         dubisectorvecz(k) = dubisectorvecz(k) * cbisectorvec(k)
         dvbisectorvecx(k) = dvbisectorvecx(k) * cbisectorvec(k)
         dvbisectorvecy(k) = dvbisectorvecy(k) * cbisectorvec(k)
         dvbisectorvecz(k) = dvbisectorvecz(k) * cbisectorvec(k)
       enddo
c   
c     force distribution for the Z-Bisect local coordinate method
c   
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         invursinvec(k) = (1.0_ti_p - urcosvec(k)**2) ** ( - half )
         vssinvec(k)    = (1.0_ti_p - vscosvec(k)**2) **   half
         wssinvec(k)    = (1.0_ti_p - wscosvec(k)**2) **   half
         duzbisectvecx(k) =  urvecx(k) * dphidrvec(k) * invusizvec(k)
     &                                 * invursinvec(k)
     &                    +  usvecx(k) * dphidsvec(k) * invusizvec(k)
         duzbisectvecy(k) =  urvecy(k) * dphidrvec(k) * invusizvec(k)
     &                                 * invursinvec(k)
     &                    +  usvecy(k) * dphidsvec(k) * invusizvec(k)
         duzbisectvecz(k) =  urvecz(k) * dphidrvec(k) * invusizvec(k)
     &                                 * invursinvec(k)
     &                    +  usvecz(k) * dphidsvec(k) * invusizvec(k)
       enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dvzbisectvecx(k) =  (  vssinvec(k) * svecx(k)
     &                        - vscosvec(k) * t1vecx(k))
     &                     * dphiduvec(k) * invvsizvec(k)
     &                     * invutsinvec(k)
         dvzbisectvecy(k) =  (  vssinvec(k) * svecy(k)
     &                        - vscosvec(k) * t1vecy(k))
     &                     * dphiduvec(k) * invvsizvec(k)
     &                     * invutsinvec(k)
         dvzbisectvecz(k) =  (  vssinvec(k) * svecz(k)
     &                        - vscosvec(k) * t1vecz(k))
     &                     * dphiduvec(k) * invvsizvec(k)
     &                     * invutsinvec(k)
       enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dwzbisectvecx(k) =  (  wssinvec(k) * svecx(k)
     &                        - wscosvec(k) * t2vecx(k))
     &                     * dphiduvec(k) * invwsizvec(k)
     &                     * invutsinvec(k)
         dwzbisectvecy(k) =  (  wssinvec(k) * svecy(k)
     &                        - wscosvec(k) * t2vecy(k))
     &                     * dphiduvec(k) * invwsizvec(k)
     &                     * invutsinvec(k)
         dwzbisectvecz(k) =  (  wssinvec(k) * svecz(k)
     &                        - wscosvec(k) * t2vecz(k))
     &                     * dphiduvec(k) * invwsizvec(k)
     &                     * invutsinvec(k)
       enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duzbisectvecx(k) = duzbisectvecx(k) *  czbisectvec(k)
         duzbisectvecy(k) = duzbisectvecy(k) *  czbisectvec(k)
         duzbisectvecz(k) = duzbisectvecz(k) *  czbisectvec(k)
         dvzbisectvecx(k) = dvzbisectvecx(k) *  czbisectvec(k)
         dvzbisectvecy(k) = dvzbisectvecy(k) *  czbisectvec(k)
         dvzbisectvecz(k) = dvzbisectvecz(k) *  czbisectvec(k)
         dwzbisectvecx(k) = dwzbisectvecx(k) *  czbisectvec(k)
         dwzbisectvecy(k) = dwzbisectvecy(k) *  czbisectvec(k)
         dwzbisectvecz(k) = dwzbisectvecz(k) *  czbisectvec(k)
      enddo

c
c     force distribution for the 3-Fold local coordinate method
c

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         pvecx(k)        =  uvecx(k) + vvecx(k) + wvecx(k)
         pvecy(k)        =  uvecy(k) + vvecy(k) + wvecy(k)
         pvecz(k)        =  uvecz(k) + vvecz(k) + wvecz(k)
         invpsizvec(k)   =  (  pvecx(k)**2 + pvecy(k)**2 + pvecz(k)**2
     &                      ) ** ( - half )

         pvecx(k)        =  pvecx(k) * invpsizvec(k)
         pvecy(k)        =  pvecy(k) * invpsizvec(k)
         pvecz(k)        =  pvecz(k) * invpsizvec(k)
         upcosvec(k)     =  uvecx(k) * pvecx(k) + uvecy(k) * pvecy(k)
     &                    + uvecz(k) * pvecz(k)
         vpcosvec(k)     =  vvecx(k) * pvecx(k) + vvecy(k) * pvecy(k)
     &                    + vvecz(k) * pvecz(k)
         wpcosvec(k)     =  wvecx(k) * pvecx(k) + wvecy(k) * pvecy(k)
     &                    + wvecz(k) * pvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)        =  uvecx(k) + vvecx(k)
         rvecy(k)        =  uvecy(k) + vvecy(k)
         rvecz(k)        =  uvecz(k) + vvecz(k)
         invrsizvec(k)   =  (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                      ) ** ( - half )

         rvecx(k)        =  rvecx(k) * invrsizvec(k)
         rvecy(k)        =  rvecy(k) * invrsizvec(k)
         rvecz(k)        =  rvecz(k) * invrsizvec(k)
         rwcosvec(k)     =  rvecx(k) * wvecx(k) + rvecy(k) * wvecy(k)
     &                    + rvecz(k) * wvecz(k)
         invrwsinvec(k)  = (1.0_ti_p - rwcosvec(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphidrvec(k)    = - trqvecx(k) * rvecx(k)
     &                     - trqvecy(k) * rvecy(k)
     &                     - trqvecz(k) * rvecz(k)
         delvecx(k)      =   rvecy(k) * wvecz(k) - rvecz(k) * wvecy(k)
         delvecy(k)      =   rvecz(k) * wvecx(k) - rvecx(k) * wvecz(k)
         delvecz(k)      =   rvecx(k) * wvecy(k) - rvecy(k) * wvecx(k)
         invdelsizvec(k) =   (  delvecx(k)**2 + delvecy(k)**2
     &                        + delvecz(k)**2 ) ** ( - half )
         delvecx(k)      =   delvecx(k) * invdelsizvec(k)
         delvecy(k)      =   delvecy(k) * invdelsizvec(k)
         delvecz(k)      =   delvecz(k) * invdelsizvec(k)
         dphiddelvec(k)  = - trqvecx(k) * delvecx(k)
     &                     - trqvecy(k) * delvecy(k)
     &                     - trqvecz(k) * delvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         epsvecx(k)     = delvecy(k) * wvecz(k) - delvecz(k) * wvecy(k)
         epsvecy(k)     = delvecz(k) * wvecx(k) - delvecx(k) * wvecz(k)
         epsvecz(k)     = delvecx(k) * wvecy(k) - delvecy(k) * wvecx(k)
         dw3foldvecx(k) =   delvecx(k)    * dphidrvec(k)
     &                    * invwsizvec(k) * invrwsinvec(k)
     &                   +  epsvecx(k)    * dphiddelvec(k) * wpcosvec(k)
     &                    * invwsizvec(k) * invpsizvec(k) 
         dw3foldvecy(k) =   delvecy(k)    * dphidrvec  (k)
     &                    * invwsizvec(k) * invrwsinvec(k)
     &                   +  epsvecy(k)    * dphiddelvec(k) * wpcosvec(k)
     &                    * invwsizvec(k) * invpsizvec(k) 
         dw3foldvecz(k) =   delvecz(k)    * dphidrvec(k)  
     &                    * invwsizvec(k) * invrwsinvec(k)
     &                   +  epsvecz(k)    * dphiddelvec(k) * wpcosvec(k)
     &                    * invwsizvec(k) * invpsizvec(k) 
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)       =    vvecx(k) + wvecx(k)
         rvecy(k)       =    vvecy(k) + wvecy(k)
         rvecz(k)       =    vvecz(k) + wvecz(k)
         invrsizvec(k)  =   (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                      ) ** ( - half )
 
         rvecx(k)       =   rvecx(k) * invrsizvec(k)
         rvecy(k)       =   rvecy(k) * invrsizvec(k)
         rvecz(k)       =   rvecz(k) * invrsizvec(k)
         rucosvec(k)    =   rvecx(k) * uvecx(k) + rvecy(k) * uvecy(k)
     &                    + rvecz(k) * uvecz(k)
         invrusinvec(k) = (1.0_ti_p - rucosvec(k)**2) ** ( - half )

         dphidrvec(k)   = - trqvecx(k) * rvecx(k)
     &                    - trqvecy(k) * rvecy(k)
     &                    - trqvecz(k) * rvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop

         delvecx(k)      =   rvecy(k) * uvecz(k) - rvecz(k) * uvecy(k)
         delvecy(k)      =   rvecz(k) * uvecx(k) - rvecx(k) * uvecz(k)
         delvecz(k)      =   rvecx(k) * uvecy(k) - rvecy(k) * uvecx(k)
         invdelsizvec(k) =  (  delvecx(k)**2 + delvecy(k)**2
     &                       + delvecz(k)**2 ) ** ( - half )
         delvecx(k)      =   delvecx(k) * invdelsizvec(k)
         delvecy(k)      =   delvecy(k) * invdelsizvec(k)
         delvecz(k)      =   delvecz(k) * invdelsizvec(k)
         dphiddelvec(k)  = - trqvecx(k) * delvecx(k)
     &                     - trqvecy(k) * delvecy(k)
     &                     - trqvecz(k) * delvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         epsvecx(k)     = delvecy(k) * uvecz(k) - delvecz(k) * uvecy(k)
         epsvecy(k)     = delvecz(k) * uvecx(k) - delvecx(k) * uvecz(k)
         epsvecz(k)     = delvecx(k) * uvecy(k) - delvecy(k) * uvecx(k)

         du3foldvecx(k) =   delvecx(k)    * dphidrvec(k)
     &                    * invusizvec(k) * invrusinvec(k)
     &                   +  epsvecx(k)    * dphiddelvec(k) * upcosvec(k)
     &                    * invusizvec(k) * invpsizvec(k) 
         du3foldvecy(k) =   delvecy(k)    * dphidrvec(k)
     &                    * invusizvec(k) * invrusinvec(k)
     &                   +  epsvecy(k)    * dphiddelvec(k) * upcosvec(k)
     &                    * invusizvec(k) * invpsizvec(k)
         du3foldvecz(k) =   delvecz(k)    * dphidrvec(k)
     &                    * invusizvec(k) * invrusinvec(k)
     &                   +  epsvecz(k)    * dphiddelvec(k) * upcosvec(k)
     &                    * invusizvec(k) * invpsizvec(k) 

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)       =  uvecx(k) + wvecx(k)
         rvecy(k)       =  uvecy(k) + wvecy(k)
         rvecz(k)       =  uvecz(k) + wvecz(k)
         invrsizvec(k)  = (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                    ) ** ( - half )

         rvecx(k)       =  rvecx(k) * invrsizvec(k)
         rvecy(k)       =  rvecy(k) * invrsizvec(k)
         rvecz(k)       =  rvecz(k) * invrsizvec(k)
         rvcosvec(k)    =  rvecx(k) * vvecx(k) + rvecy(k) * vvecy(k)
     &                   + rvecz(k) * vvecz(k)
         invrvsinvec(k) = (1.0_ti_p - rvcosvec(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphidrvec(k)    = - trqvecx(k) * rvecx(k)
     &                     - trqvecy(k) * rvecy(k)
     &                     - trqvecz(k) * rvecz(k)

         delvecx(k)      =  rvecy(k) * vvecz(k) - rvecz(k) * vvecy(k)
         delvecy(k)      =  rvecz(k) * vvecx(k) - rvecx(k) * vvecz(k)
         delvecz(k)      =  rvecx(k) * vvecy(k) - rvecy(k) * vvecx(k)

         invdelsizvec(k) = (  delvecx(k)**2 + delvecy(k)**2
     &                      + delvecz(k)**2 ) ** ( - half )
         delvecx(k)      =  delvecx(k) * invdelsizvec(k)
         delvecy(k)      =  delvecy(k) * invdelsizvec(k)
         delvecz(k)      =  delvecz(k) * invdelsizvec(k)

         dphiddelvec(k)  = - trqvecx(k) * delvecx(k)
     &                     - trqvecy(k) * delvecy(k)
     &                     - trqvecz(k) * delvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         epsvecx(k)     = delvecy(k) * vvecz(k) - delvecz(k) * vvecy(k)
         epsvecy(k)     = delvecz(k) * vvecx(k) - delvecx(k) * vvecz(k)
         epsvecz(k)     = delvecx(k) * vvecy(k) - delvecy(k) * vvecx(k)

         dv3foldvecx(k) =   delvecx(k)    * dphidrvec(k)
     &                    * invvsizvec(k) * invrvsinvec(k)
     &                   +  epsvecx(k)    * dphiddelvec(k) * vpcosvec(k)
     &                    * invvsizvec(k) * invpsizvec(k)
         dv3foldvecy(k) =   delvecy(k)    * dphidrvec(k)
     &                    * invvsizvec(k) * invrvsinvec(k)
     &                   +  epsvecy(k)    * dphiddelvec(k) * vpcosvec(k)
     &                    * invvsizvec(k) * invpsizvec(k)
         dv3foldvecz(k) =   delvecz(k)    * dphidrvec(k)
     &                    * invvsizvec(k) * invrvsinvec(k)
     &                   +  epsvecz(k)    * dphiddelvec(k) * vpcosvec(k)
     &                    * invvsizvec(k) * invpsizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         du3foldvecx(k) = du3foldvecx(k) * c3foldvec(k)
         du3foldvecy(k) = du3foldvecy(k) * c3foldvec(k)
         du3foldvecz(k) = du3foldvecz(k) * c3foldvec(k)
         dv3foldvecx(k) = dv3foldvecx(k) * c3foldvec(k)
         dv3foldvecy(k) = dv3foldvecy(k) * c3foldvec(k)
         dv3foldvecz(k) = dv3foldvecz(k) * c3foldvec(k)
         dw3foldvecx(k) = dw3foldvecx(k) * c3foldvec(k)
         dw3foldvecy(k) = dw3foldvecy(k) * c3foldvec(k)
         dw3foldvecz(k) = dw3foldvecz(k) * c3foldvec(k)
      enddo
!DIR$ ASSUME_ALIGNED duvecx : 64
!DIR$ ASSUME_ALIGNED duvecy : 64
!DIR$ ASSUME_ALIGNED duvecz : 64
!DIR$ ASSUME_ALIGNED dvvecx : 64
!DIR$ ASSUME_ALIGNED dvvecy : 64
!DIR$ ASSUME_ALIGNED dvvecz : 64
!DIR$ ASSUME_ALIGNED dwvecx : 64
!DIR$ ASSUME_ALIGNED dwvecy : 64
!DIR$ ASSUME_ALIGNED dwvecz : 64
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duvecx(k) =  duzonlyvecx(k)    + duzthenxvecx(k) 
     &              + dubisectorvecx(k) + duzbisectvecx(k) 
     &              + du3foldvecx(k)   

         duvecy(k) =  duzonlyvecy(k)    + duzthenxvecy(k)
     &              + dubisectorvecy(k) + duzbisectvecy(k) 
     &              + du3foldvecy(k)   

         duvecz(k) =  duzonlyvecz(k)    + duzthenxvecz(k)
     &              + dubisectorvecz(k) + duzbisectvecz(k)
     &              + du3foldvecz(k) 
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dvvecx(k) =  dvzonlyvecx(k)    + dvzthenxvecx(k)   
     &              + dvbisectorvecx(k) + dvzbisectvecx(k)  
     &              + dv3foldvecx(k)    

         dvvecy(k) =  dvzonlyvecy(k)    + dvzthenxvecy(k)   
     &              + dvbisectorvecy(k) + dvzbisectvecy(k)  
     &           + dv3foldvecy(k)    

         dvvecz(k) =  dvzonlyvecz(k)    + dvzthenxvecz(k)   
     &              + dvbisectorvecz(k) + dvzbisectvecz(k)  
     &              + dv3foldvecz(k)    
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dwvecx(k) =  dwzonlyvecx(k)    + dwzthenxvecx(k)   
     &              + dwbisectorvecx(k) + dwzbisectvecx(k)  
     &              + dw3foldvecx(k)    

         dwvecy(k) =  dwzonlyvecy(k)    + dwzthenxvecy(k)   
     &              + dwbisectorvecy(k) + dwzbisectvecy(k)  
     &              + dw3foldvecy(k)    

         dwvecz(k) =  dwzonlyvecz(k)    + dwzthenxvecz(k)   
     &              + dwbisectorvecz(k) + dwzbisectvecz(k)  
     &              + dw3foldvecz(k)    
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(ialocvec(k)) =  devecx(ialocvec(k)) + duvecx(k)
         devecy(ialocvec(k)) =  devecy(ialocvec(k)) + duvecy(k)
         devecz(ialocvec(k)) =  devecz(ialocvec(k)) + duvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(iblocvec(k)) =  devecx(iblocvec(k))
     &                          - duvecx(k) -dvvecx(k) -dwvecx(k)
         devecy(iblocvec(k)) =  devecy(iblocvec(k))
     &                          - duvecy(k) -dvvecy(k) -dwvecy(k)
         devecz(iblocvec(k)) =  devecz(iblocvec(k))
     &                          - duvecz(k) -dvvecz(k) -dwvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(iclocvec(k)) =  devecx(iclocvec(k)) + dvvecx(k)
         devecy(iclocvec(k)) =  devecy(iclocvec(k)) + dvvecy(k)
         devecz(iclocvec(k)) =  devecz(iclocvec(k)) + dvvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(idlocvec(k)) =  devecx(idlocvec(k)) + dwvecx(k)
         devecy(idlocvec(k)) =  devecy(idlocvec(k)) + dwvecy(k)
         devecz(idlocvec(k)) =  devecz(idlocvec(k)) + dwvecz(k)
      enddo
      call timer_exit( timer_torque )
      end
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine torquevecrec -- convert single site torque to force  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "torque2vec" takes the torque values on multiple sites defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     npolelocnl and npolelocnlloop came from mpole module
c
c     INPUT  : ivec is indices of sites 
c              trqvec are trqs of sites 
c              devecx,devecy and devecz are deps coming from epolar
c              or empole
c
c     OUTPUT : dvvecx, dwvecx and duvecx, etc... are forces 
c              devecx,devecy and devecz are deps sent back to epolar
c              or empole
      subroutine torquevec2_rec (ivec,trqvecx,trqvecy,trqvecz,
     &                           dvvecx,dvvecy,dvvecz,
     &                           dwvecx,dwvecy,dwvecz,
     &                           duvecx,duvecy,duvecz,
     &                           devecx,devecy,devecz)
      use atoms
      use deriv
      use domdec
      use mpole
      use sizes
      use mpi
      use timestat
      use tinheader ,only: ti_p,re_p
      implicit none
      integer i,j,k,kk
      integer npolelocnl1
      integer ivec(npolelocnlloop)
      real(t_p) sqrt3over2
      real(t_p) half
      real(t_p) trqvecx(npolelocnlloop)
      real(t_p) trqvecy(npolelocnlloop)
      real(t_p) trqvecz(npolelocnlloop)
      real(t_p) dvvecx(npolelocnlloop)
      real(t_p) dvvecy(npolelocnlloop)
      real(t_p) dvvecz(npolelocnlloop)
      real(t_p) dwvecx(npolelocnlloop)
      real(t_p) dwvecy(npolelocnlloop)
      real(t_p) dwvecz(npolelocnlloop)
      real(t_p) duvecx(npolelocnlloop)
      real(t_p) duvecy(npolelocnlloop)
      real(t_p) duvecz(npolelocnlloop)
      real(t_p) devecx(*)
      real(t_p) devecy(*)
      real(t_p) devecz(*)
!DIR$ ATTRIBUTES ALIGN:64::axetypvec
      integer axetypvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec1
      integer ivec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::iavec,ibvec
      integer iavec(npolelocnlloop),ibvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::icvec,idvec
      integer icvec(npolelocnlloop),idvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ialocvec,iblocvec
      integer ialocvec(npolelocnlloop),iblocvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::iclocvec,idlocvec
      integer iclocvec(npolelocnlloop),idlocvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: delvecx,delvecy
      real(t_p) delvecx(npolelocnlloop)
      real(t_p) delvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: delvecz
      real(t_p) delvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::epsvecx,epsvecy,epsvecz
      real(t_p) epsvecx(npolelocnlloop)
      real(t_p) epsvecy(npolelocnlloop)
      real(t_p) epsvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invdelsizvec,epssizvec
      real(t_p) invdelsizvec(npolelocnlloop)
      real(t_p) epssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::pvecx,pvecy,pvecz,invpsizvec
      real(t_p) pvecx(npolelocnlloop)
      real(t_p) pvecy(npolelocnlloop)
      real(t_p) pvecz(npolelocnlloop)
      real(t_p) invpsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rvecx,rvecy,rvecz
      real(t_p) rvecx(npolelocnlloop)
      real(t_p) rvecy(npolelocnlloop)
      real(t_p) rvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invrsizvec
      real(t_p) invrsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::svecx,svecy,svecz
      real(t_p) svecx(npolelocnlloop)
      real(t_p) svecy(npolelocnlloop)
      real(t_p) svecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invssizvec
      real(t_p) invssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::invusizvec
      real(t_p) invusizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uvecx,uvecy,uvecz
      real(t_p) uvecx(npolelocnlloop)
      real(t_p) uvecy(npolelocnlloop)
      real(t_p) uvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::invvsizvec
      real(t_p) invvsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vvecx,vvecy,vvecz
      real(t_p) vvecx(npolelocnlloop)
      real(t_p) vvecy(npolelocnlloop)
      real(t_p) vvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::invwsizvec
      real(t_p) invwsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wvecx,wvecy,wvecz
      real(t_p) wvecx(npolelocnlloop)
      real(t_p) wvecy(npolelocnlloop)
      real(t_p) wvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uvvecx,uvvecy,uvvecz
      real(t_p) uvvecx(npolelocnlloop)
      real(t_p) uvvecy(npolelocnlloop)
      real(t_p) uvvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invuvsizvec
      real(t_p) invuvsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uwvecx,uwvecy,uwvecz
      real(t_p) uwvecx(npolelocnlloop)
      real(t_p) uwvecy(npolelocnlloop)
      real(t_p) uwvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invuwsizvec
      real(t_p) invuwsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vwvecx,vwvecy,vwvecz
      real(t_p) vwvecx(npolelocnlloop)
      real(t_p) vwvecy(npolelocnlloop)
      real(t_p) vwvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invvwsizvec
      real(t_p) invvwsizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::usvecx,usvecy,usvecz
      real(t_p) usvecx(npolelocnlloop)
      real(t_p) usvecy(npolelocnlloop)
      real(t_p) usvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invussizvec
      real(t_p) invussizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vsvecx,vsvecy,vsvecz
      real(t_p) vsvecx(npolelocnlloop)
      real(t_p) vsvecy(npolelocnlloop)
      real(t_p) vsvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invvssizvec
      real(t_p) invvssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wsvecx,wsvecy,wsvecz
      real(t_p) wsvecx(npolelocnlloop)
      real(t_p) wsvecy(npolelocnlloop)
      real(t_p) wsvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invwssizvec
      real(t_p) invwssizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::urvecx,urvecy,urvecz
      real(t_p) urvecx(npolelocnlloop)
      real(t_p) urvecy(npolelocnlloop)
      real(t_p) urvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invursizvec
      real(t_p) invursizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::t1vecx,t1vecy,t1vecz
      real(t_p) t1vecx(npolelocnlloop)
      real(t_p) t1vecy(npolelocnlloop)
      real(t_p) t1vecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invt1sizvec
      real(t_p) invt1sizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::t2vecx,t2vecy,t2vecz
      real(t_p) t2vecx(npolelocnlloop)
      real(t_p) t2vecy(npolelocnlloop)
      real(t_p) t2vecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: invt2sizvec
      real(t_p) invt2sizvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rucosvec,invrusinvec
      real(t_p) rucosvec(npolelocnlloop)
      real(t_p) invrusinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rvcosvec,invrvsinvec
      real(t_p) rvcosvec(npolelocnlloop)
      real(t_p) invrvsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::rwcosvec,invrwsinvec
      real(t_p) rwcosvec(npolelocnlloop)
      real(t_p) invrwsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uvcosvec,invuvsinvec
      real(t_p) uvcosvec(npolelocnlloop)
      real(t_p) invuvsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uwcosvec,uwsinvec
      real(t_p) uwcosvec(npolelocnlloop)
      real(t_p) uwsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vwcosvec,vwsinvec
      real(t_p) vwcosvec(npolelocnlloop)
      real(t_p) vwsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::upcosvec
      real(t_p) upcosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::urcosvec,invursinvec
      real(t_p) urcosvec(npolelocnlloop)
      real(t_p) invursinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::uscosvec
      real(t_p) uscosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vpcosvec
      real(t_p) vpcosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::vscosvec,vssinvec
      real(t_p) vscosvec(npolelocnlloop)
      real(t_p) vssinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wpcosvec
      real(t_p) wpcosvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::wscosvec,wssinvec
      real(t_p) wscosvec(npolelocnlloop)
      real(t_p) wssinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ut1cosvec,ut1sinvec
      real(t_p) ut1cosvec(npolelocnlloop)
      real(t_p) ut1sinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::ut2cosvec,ut2sinvec,invutsinvec
      real(t_p) ut2cosvec(npolelocnlloop)
      real(t_p) ut2sinvec(npolelocnlloop)
      real(t_p) invutsinvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphiduvec,dphidvvec
      real(t_p) dphiduvec(npolelocnlloop)
      real(t_p) dphidvvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphidwvec,dphidrvec
      real(t_p) dphidwvec(npolelocnlloop)
      real(t_p) dphidrvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphiddelvec
      real(t_p) dphiddelvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::dphidsvec
      real(t_p) dphidsvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xiavec,yiavec,ziavec
      real(t_p) xiavec(npolelocnlloop)
      real(t_p) yiavec(npolelocnlloop)
      real(t_p) ziavec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xibvec,yibvec,zibvec
      real(t_p) xibvec(npolelocnlloop)
      real(t_p) yibvec(npolelocnlloop)
      real(t_p) zibvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xicvec,yicvec,zicvec
      real(t_p) xicvec(npolelocnlloop)
      real(t_p) yicvec(npolelocnlloop)
      real(t_p) zicvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::xidvec,yidvec,zidvec
      real(t_p) xidvec(npolelocnlloop)
      real(t_p) yidvec(npolelocnlloop)
      real(t_p) zidvec(npolelocnlloop)
c     3-Fold
!DIR$ ATTRIBUTES ALIGN:64:: c3foldvec
      real(t_p) c3foldvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecx
      real(t_p) du3foldvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecy
      real(t_p) du3foldvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecz
      real(t_p) du3foldvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecx
      real(t_p) dv3foldvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecy
      real(t_p) dv3foldvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecz
      real(t_p) dv3foldvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecx
      real(t_p) dw3foldvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecy
      real(t_p) dw3foldvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecz
      real(t_p) dw3foldvecz(npolelocnlloop)
c     Bisector
!DIR$ ATTRIBUTES ALIGN:64:: cbisectorvec
      real(t_p) cbisectorvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecx
      real(t_p) dubisectorvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecy
      real(t_p) dubisectorvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecz
      real(t_p) dubisectorvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecx
      real(t_p) dvbisectorvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecy
      real(t_p) dvbisectorvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecz
      real(t_p) dvbisectorvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecx
      real(t_p) dwbisectorvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecy
      real(t_p) dwbisectorvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecz
      real(t_p) dwbisectorvecz(npolelocnlloop)
c     Z-Bisect
!DIR$ ATTRIBUTES ALIGN:64:: czbisectvec
      real(t_p) czbisectvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecx
      real(t_p) duzbisectvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecy
      real(t_p) duzbisectvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecz
      real(t_p) duzbisectvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecx
      real(t_p) dvzbisectvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecy
      real(t_p) dvzbisectvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecz
      real(t_p) dvzbisectvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecx
      real(t_p) dwzbisectvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecy
      real(t_p) dwzbisectvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecz
      real(t_p) dwzbisectvecz(npolelocnlloop)
C     Z-Only
!DIR$ ATTRIBUTES ALIGN:64:: czonlyvec
      real(t_p) czonlyvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecx
      real(t_p) duzonlyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecy
      real(t_p) duzonlyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecz
      real(t_p) duzonlyvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecx
      real(t_p) dvzonlyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecy
      real(t_p) dvzonlyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecz
      real(t_p) dvzonlyvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecx
      real(t_p) dwzonlyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecy
      real(t_p) dwzonlyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecz
      real(t_p) dwzonlyvecz(npolelocnlloop)
C     Z-then-X
!DIR$ ATTRIBUTES ALIGN:64:: czthenxvec
      real(t_p) czthenxvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecx
      real(t_p) duzthenxvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecy
      real(t_p) duzthenxvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecz
      real(t_p) duzthenxvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecx
      real(t_p) dvzthenxvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecy
      real(t_p) dvzthenxvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecz
      real(t_p) dvzthenxvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecx
      real(t_p) dwzthenxvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecy
      real(t_p) dwzthenxvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecz
      real(t_p) dwzthenxvecz(npolelocnlloop)
      character*8 paxe(npolelocnlloop)

c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'torquevec2rec'
      call timer_exit( timer_torque )
      sqrt3over2 = sqrt(3.0_ti_p)/2.0_ti_p
      half       = 0.5_ti_p
c
c     zero out force components on local frame-defining atoms
c
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         iavec(k)     = 1
         ibvec(k)     = 1
         icvec(k)     = 1
         idvec(k)     = 1
         ialocvec(k)  = 1
         iblocvec(k)  = 1
         iclocvec(k)  = 1
         idlocvec(k)  = 1
         ivec1(k)     = 1
         axetypvec(k) = 0
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         c3foldvec(k)   = 0.0_ti_p
         du3foldvecx(k) = 0.0_ti_p
         du3foldvecy(k) = 0.0_ti_p
         du3foldvecz(k) = 0.0_ti_p
         dv3foldvecx(k) = 0.0_ti_p
         dv3foldvecy(k) = 0.0_ti_p
         dv3foldvecz(k) = 0.0_ti_p
         dw3foldvecx(k) = 0.0_ti_p
         dw3foldvecy(k) = 0.0_ti_p
         dw3foldvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         cbisectorvec(k)   = 0.0_ti_p
         dubisectorvecx(k) = 0.0_ti_p
         dubisectorvecy(k) = 0.0_ti_p
         dubisectorvecz(k) = 0.0_ti_p
         dvbisectorvecx(k) = 0.0_ti_p
         dvbisectorvecy(k) = 0.0_ti_p
         dvbisectorvecz(k) = 0.0_ti_p
         dwbisectorvecx(k) = 0.0_ti_p
         dwbisectorvecy(k) = 0.0_ti_p
         dwbisectorvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         czbisectvec(k)   = 0.0_ti_p
         duzbisectvecx(k) = 0.0_ti_p
         duzbisectvecy(k) = 0.0_ti_p
         duzbisectvecz(k) = 0.0_ti_p
         dvzbisectvecx(k) = 0.0_ti_p
         dvzbisectvecy(k) = 0.0_ti_p
         dvzbisectvecz(k) = 0.0_ti_p
         dwzbisectvecx(k) = 0.0_ti_p
         dwzbisectvecy(k) = 0.0_ti_p
         dwzbisectvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         czonlyvec(k)   = 0.0_ti_p
         duzonlyvecx(k) = 0.0_ti_p
         duzonlyvecy(k) = 0.0_ti_p
         duzonlyvecz(k) = 0.0_ti_p
         dvzonlyvecx(k) = 0.0_ti_p
         dvzonlyvecy(k) = 0.0_ti_p
         dvzonlyvecz(k) = 0.0_ti_p
         dwzonlyvecx(k) = 0.0_ti_p
         dwzonlyvecy(k) = 0.0_ti_p
         dwzonlyvecz(k) = 0.0_ti_p
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         czthenxvec(k)   = 0.0_ti_p
         duzthenxvecx(k) = 0.0_ti_p
         duzthenxvecy(k) = 0.0_ti_p
         duzthenxvecz(k) = 0.0_ti_p
         dvzthenxvecx(k) = 0.0_ti_p
         dvzthenxvecy(k) = 0.0_ti_p
         dvzthenxvecz(k) = 0.0_ti_p
         dwzthenxvecx(k) = 0.0_ti_p
         dwzthenxvecy(k) = 0.0_ti_p
         dwzthenxvecz(k) = 0.0_ti_p
      enddo

c
c     get the local frame type and the frame-defining atoms
c
      do k = 1, npolelocnl
         select case (polaxe(ivec(k)))
                case ('3-Fold'  ) ;  axetypvec(k) = 1
                                     c3foldvec(k)= 1.0_ti_p
                case ('Bisector') ;  axetypvec(k) = 2
                                     cbisectorvec(k)= 1.0_ti_p
                case ('Z-Bisect') ;  axetypvec(k) = 3
                                     czbisectvec(k)= 1.0_ti_p
                case ('Z-Only'  ) ;  axetypvec(k) = 4
                                     czonlyvec(k)= 1.0_ti_p
                case ('Z-then-X') ;  axetypvec(k) = 5
                                     czthenxvec(k)= 1.0_ti_p
                case ('None'    ) ;  axetypvec(k) = 0
                case default      ;  axetypvec(k) = 0
         endselect
      enddo

      kk = 0
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
      do k = 1,npolelocnlloop
         if (axetypvec(k).ne.0) then
            kk = kk +1
            ivec1 (kk) = ivec (k)
         endif
      enddo
      npolelocnl1 = kk
      if (npolelocnl1.eq.0) return ! All axetyp are 0, so return
    
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         iavec(k) = zaxis(ivec1(k))
         if(iavec(k).gt.0) then
            ialocvec(k) = locrec1(iavec(k))
         else
            iavec(k) = 1
         endif
         if(ialocvec(k).eq.0)  ialocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         ibvec(k) = ipole(ivec1(k))
         if(ibvec(k).gt.0) then
            iblocvec(k) = locrec1(ibvec(k))
         else
            ibvec(k) = 1
         endif
         if(iblocvec(k).eq.0)  iblocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         icvec(k) = xaxis(ivec1(k))
         if(icvec(k).gt.0) then
            iclocvec(k) = locrec1(icvec(k))
         else
            icvec(k) = 1
         endif
         if(iclocvec(k).eq.0)  iclocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         idvec(k) = yaxis(ivec1(k))
         if(idvec(k).gt.0) then
            idlocvec(k) = locrec1(idvec(k))
         else
            idvec(k) = 1
         endif
         if(idlocvec(k).eq.0)  idlocvec(k) = 1
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         xiavec(k)=x(iavec(k))
         yiavec(k)=y(iavec(k))
         ziavec(k)=z(iavec(k))
         xibvec(k)=x(ibvec(k))
         yibvec(k)=y(ibvec(k))
         zibvec(k)=z(ibvec(k))
         xicvec(k)=x(icvec(k))
         yicvec(k)=y(icvec(k))
         zicvec(k)=z(icvec(k))
         xidvec(k)=x(idvec(k))
         yidvec(k)=y(idvec(k))
         zidvec(k)=z(idvec(k))
      enddo
c
c     construct the three rotation axes for the local frame
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvecx(k)  = xiavec(k)- xibvec(k)
         uvecy(k)  = yiavec(k)- yibvec(k)
         uvecz(k)  = ziavec(k)- zibvec(k)
         invusizvec(k) = (  uvecx(k)**2 + uvecy(k)**2 + uvecz(k)**2
     &                   ) ** ( - half )
         uvecx(k) = uvecx(k) * invusizvec(k)
         uvecy(k) = uvecy(k) * invusizvec(k)
         uvecz(k) = uvecz(k) * invusizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vvecx(k)  = xicvec(k) - xibvec(k)
         vvecy(k)  = yicvec(k) - yibvec(k)
         vvecz(k)  = zicvec(k) - zibvec(k)
         invvsizvec(k) = (  vvecx(k)**2 + vvecy(k)**2 + vvecz(k)**2
     &                   ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         if (axetypvec(k).eq.4) then ! 'Z-Only' 
            if (abs(uvecx(k)).gt.sqrt3over2) then
               vvecx(k)  = 0.0_ti_p
            else
               vvecx(k)  = 1.0_ti_p
            endif
            vvecy(k) = 1.0_ti_p - vvecx(k)
            vvecz(k)  = 0.0_ti_p
            invvsizvec(k) = 1.0_ti_p
         endif
      enddo
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         if (    axetypvec(k).eq.1  ! '3-Fold'
     &       .or.axetypvec(k).eq.3) ! 'Z-Bisect'
     &   then
            wvecx(k) = xidvec(k) - xibvec(k)
            wvecy(k) = yidvec(k) - yibvec(k)
            wvecz(k) = zidvec(k) - zibvec(k)
         else
           wvecx(k) =  uvecy(k) * vvecz(k) - uvecz(k) * vvecy(k)
           wvecy(k) =  uvecz(k) * vvecx(k) - uvecx(k) * vvecz(k)
           wvecz(k) =  uvecx(k) * vvecy(k) - uvecy(k) * vvecx(k)
         endif
         invwsizvec(k) = (  wvecx(k)**2 + wvecy(k)**2
     &                    + wvecz(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vvecx(k) = vvecx(k) * invvsizvec(k)
         vvecy(k) = vvecy(k) * invvsizvec(k)
         vvecz(k) = vvecz(k) * invvsizvec(k)
         wvecx(k) = wvecx(k) * invwsizvec(k)
         wvecy(k) = wvecy(k) * invwsizvec(k)
         wvecz(k) = wvecz(k) * invwsizvec(k)
      enddo
c
c     build some additional axes needed for the Z-Bisect method
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)  = vvecx(k) + wvecx(k)
         rvecy(k)  = vvecy(k) + wvecy(k)
         rvecz(k)  = vvecz(k) + wvecz(k)
         svecx(k) =  uvecy(k) * rvecz(k) - uvecz(k) * rvecy(k)
         svecy(k) =  uvecz(k) * rvecx(k) - uvecx(k) * rvecz(k)
         svecz(k) =  uvecx(k) * rvecy(k) - uvecy(k) * rvecx(k)
         invrsizvec(k) = (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                   ) ** ( - half )
         invssizvec(k) = (  svecx(k)**2 + svecy(k)**2 + svecz(k)**2
     &                   ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k) = rvecx(k) * invrsizvec(k)
         svecx(k) = svecx(k) * invssizvec(k)
         rvecy(k) = rvecy(k) * invrsizvec(k)
         svecy(k) = svecy(k) * invssizvec(k)
         rvecz(k) = rvecz(k) * invrsizvec(k)
         svecz(k) = svecz(k) * invssizvec(k)
      enddo
c     find the perpendicular and angle for each pair of axes
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvvecx(k)  =  vvecy(k) * uvecz(k) - vvecz(k) * uvecy(k)
         uvvecy(k)  =  vvecz(k) * uvecx(k) - vvecx(k) * uvecz(k)
         uvvecz(k)  =  vvecx(k) * uvecy(k) - vvecy(k) * uvecx(k)
         invuvsizvec(k) = (  uvvecx(k)**2 + uvvecy(k)**2 + uvvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uwvecx(k)  =  wvecy(k) * uvecz(k) - wvecz(k) * uvecy(k)
         uwvecy(k)  =  wvecz(k) * uvecx(k) - wvecx(k) * uvecz(k)
         uwvecz(k)  =  wvecx(k) * uvecy(k) - wvecy(k) * uvecx(k)
         invuwsizvec(k) = (  uwvecx(k)**2 + uwvecy(k)**2 + uwvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vwvecx(k)  =  wvecy(k) * vvecz(k) - wvecz(k) * vvecy(k)
         vwvecy(k)  =  wvecz(k) * vvecx(k) - wvecx(k) * vvecz(k)
         vwvecz(k)  =  wvecx(k) * vvecy(k) - wvecy(k) * vvecx(k)
         invvwsizvec(k) = (  vwvecx(k)**2 + vwvecy(k)**2 + vwvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvvecx(k) = uvvecx(k) * invuvsizvec(k)
         uvvecy(k) = uvvecy(k) * invuvsizvec(k)
         uvvecz(k) = uvvecz(k) * invuvsizvec(k)
         uwvecx(k) = uwvecx(k) * invuwsizvec(k)
         uwvecy(k) = uwvecy(k) * invuwsizvec(k)
         uwvecz(k) = uwvecz(k) * invuwsizvec(k)
         vwvecx(k) = vwvecx(k) * invvwsizvec(k)
         vwvecy(k) = vwvecy(k) * invvwsizvec(k)
         vwvecz(k) = vwvecz(k) * invvwsizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         urvecx(k)  =  rvecy(k) * uvecz(k) - rvecz(k) * uvecy(k)
         urvecy(k)  =  rvecz(k) * uvecx(k) - rvecx(k) * uvecz(k)
         urvecz(k)  =  rvecx(k) * uvecy(k) - rvecy(k) * uvecx(k)
         invursizvec(k) = (  urvecx(k)**2 + urvecy(k)**2 + urvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         usvecx(k)  =  svecy(k) * uvecz(k) - svecz(k) * uvecy(k)
         usvecy(k)  =  svecz(k) * uvecx(k) - svecx(k) * uvecz(k)
         usvecz(k)  =  svecx(k) * uvecy(k) - svecy(k) * uvecx(k)
         invussizvec(k) = (  usvecx(k)**2 + usvecy(k)**2 + usvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vsvecx(k)  =  svecy(k) * vvecz(k) - svecz(k) * vvecy(k)
         vsvecy(k)  =  svecz(k) * vvecx(k) - svecx(k) * vvecz(k)
         vsvecz(k)  =  svecx(k) * vvecy(k) - svecy(k) * vvecx(k)
         invvssizvec(k) = (  vsvecx(k)**2 + vsvecy(k)**2 + vsvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         wsvecx(k)  =  svecy(k) * wvecz(k) - svecz(k) * wvecy(k)
         wsvecy(k)  =  svecz(k) * wvecx(k) - svecx(k) * wvecz(k)
         wsvecz(k)  =  svecx(k) * wvecy(k) - svecy(k) * wvecx(k)
         invwssizvec(k) = (  wsvecx(k)**2 + wsvecy(k)**2 + wsvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         urvecx(k)  = urvecx(k) * invursizvec(k)
         usvecx(k)  = usvecx(k) * invussizvec(k)
         vsvecx(k)  = vsvecx(k) * invvssizvec(k)
         wsvecx(k)  = wsvecx(k) * invwssizvec(k)
         urvecy(k)  = urvecy(k) * invursizvec(k)
         usvecy(k)  = usvecy(k) * invussizvec(k)
         vsvecy(k)  = vsvecy(k) * invvssizvec(k)
         wsvecy(k)  = wsvecy(k) * invwssizvec(k)
         urvecz(k)  = urvecz(k) * invursizvec(k)
         usvecz(k)  = usvecz(k) * invussizvec(k)
         vsvecz(k)  = vsvecz(k) * invvssizvec(k)
         wsvecz(k)  = wsvecz(k) * invwssizvec(k)
      enddo
c
c     get sine and cosine of angles between the rotation axes
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         uvcosvec(k) =  uvecx(k) * vvecx(k) + uvecy(k) * vvecy(k)
     &                + uvecz(k) * vvecz(k)
         uwcosvec(k) =  uvecx(k) * wvecx(k) + uvecy(k) * wvecy(k)
     &                + uvecz(k) * wvecz(k)
         vwcosvec(k) =  vvecx(k) * wvecx(k) + vvecy(k) * wvecy(k)
     &                + vvecz(k) * wvecz(k)
         invuvsinvec(k) = (1.0_ti_p - uvcosvec(k)**2) ** ( - half )
      enddo
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         urcosvec(k) =  uvecx(k) * rvecx(k) + uvecy(k) * rvecy(k)
     &                + uvecz(k) * rvecz(k)
         uscosvec(k) =  uvecx(k) * svecx(k) + uvecy(k) * svecy(k)
     &                + uvecz(k) * svecz(k)
         vscosvec(k) =  vvecx(k) * svecx(k) + vvecy(k) * svecy(k)
     &                + vvecz(k) * svecz(k)
         wscosvec(k) =  wvecx(k) * svecx(k) + wvecy(k) * svecy(k)
     &                + wvecz(k) * svecz(k)
      enddo
c
c     compute the projection of v and w onto the ru-plane
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         t1vecx(k) =  vvecx(k) - svecx(k) * vscosvec(k)
         t1vecy(k) =  vvecy(k) - svecy(k) * vscosvec(k)
         t1vecz(k) =  vvecz(k) - svecz(k) * vscosvec(k)
         t2vecx(k) =  wvecx(k) - svecx(k) * wscosvec(k)
         t2vecy(k) =  wvecy(k) - svecy(k) * wscosvec(k)
         t2vecz(k) =  wvecz(k) - svecz(k) * wscosvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         invt1sizvec(k) = (  t1vecx(k)**2 + t1vecy(k)**2 + t1vecz(k)**2
     &                    ) ** ( - half )
         invt2sizvec(k) = (  t2vecx(k)**2 + t2vecy(k)**2 + t2vecz(k)**2
     &                    ) ** ( - half )
         t1vecx(k) = t1vecx(k) * invt1sizvec(k)
         t1vecy(k) = t1vecy(k) * invt1sizvec(k)
         t1vecz(k) = t1vecz(k) * invt1sizvec(k)
         t2vecx(k) = t2vecx(k) * invt2sizvec(k)
         t2vecy(k) = t2vecy(k) * invt2sizvec(k)
         t2vecz(k) = t2vecz(k) * invt2sizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         ut1cosvec(k) =  uvecx(k) * t1vecx(k) + uvecy(k) * t1vecy(k)
     &                 + uvecz(k) * t1vecz(k)
         ut2cosvec(k) =  uvecx(k) * t2vecx(k) + uvecy(k) * t2vecy(k)
     &                 + uvecz(k) * t2vecz(k)
         invutsinvec(k) = (  (1.0_ti_p - ut1cosvec(k)**2) ** half
     &                     + (1.0_ti_p - ut2cosvec(k)**2) ** half
     &                    ) **(-1.0_ti_p)
      enddo
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphiduvec(k) = - trqvecx(k) * uvecx(k) - trqvecy(k) * uvecy(k)
     &                  - trqvecz(k) * uvecz(k)
         dphidvvec(k) = - trqvecx(k) * vvecx(k) - trqvecy(k) * vvecy(k)
     &                  - trqvecz(k) * vvecz(k)
         dphidwvec(k) = - trqvecx(k) * wvecx(k) - trqvecy(k) * wvecy(k)
     &                  - trqvecz(k) * wvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphidrvec(k) = - trqvecx(k) * rvecx(k) - trqvecy(k) * rvecy(k)
     &                  - trqvecz(k) * rvecz(k)
         dphidsvec(k) = - trqvecx(k) * svecx(k) - trqvecy(k) * svecy(k)
     &                  - trqvecz(k) * svecz(k)
      enddo
c
c     force distribution for the Z-Only local coordinate method
c
          
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duzonlyvecx(k) =   uvvecx(k)     * dphidvvec(k) * invusizvec(k)
     &                    * invuvsinvec(k)
     &                   +  uwvecx(k)     * dphidwvec(k) * invusizvec(k)
         duzonlyvecy(k) =   uvvecy(k)     * dphidvvec(k) * invusizvec(k)
     &                    * invuvsinvec(k)
     &                   +  uwvecy(k)     * dphidwvec(k) * invusizvec(k)
         duzonlyvecz(k) =   uvvecz(k)     * dphidvvec(k) * invusizvec(k)
     &                    * invuvsinvec(k)
     &                   +  uwvecz(k)     * dphidwvec(k) * invusizvec(k)

         duzonlyvecx(k) = duzonlyvecx(k)  * czonlyvec(k)
         duzonlyvecy(k) = duzonlyvecy(k)  * czonlyvec(k)
         duzonlyvecz(k) = duzonlyvecz(k)  * czonlyvec(k)
      enddo
c
c     force distribution for the Z-then-X local coordinate method
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duzthenxvecx(k) =   uvvecx(k) * dphidvvec(k)    * invusizvec(k)
     &                                 * invuvsinvec(k)
     &                    +  uwvecx(k) * dphidwvec(k)    * invusizvec(k)
         duzthenxvecy(k) =   uvvecy(k) * dphidvvec(k)    * invusizvec(k)
     &                                 * invuvsinvec(k)
     &                    +  uwvecy(k) * dphidwvec(k)    * invusizvec(k)
         duzthenxvecz(k) =   uvvecz(k) * dphidvvec(k)    * invusizvec(k)
     &                                 * invuvsinvec(k)
     &                    +  uwvecz(k) * dphidwvec(k)    * invusizvec(k)

         dvzthenxvecx(k) = -  uvvecx(k) * dphiduvec(k)   * invvsizvec(k)
     &                                  * invuvsinvec(k)
         dvzthenxvecy(k) = -  uvvecy(k) * dphiduvec(k)   * invvsizvec(k)
     &                                  * invuvsinvec(k)
         dvzthenxvecz(k) = -  uvvecz(k) * dphiduvec(k)   * invvsizvec(k)
     &                                  * invuvsinvec(k)

         duzthenxvecx(k) = duzthenxvecx(k) * czthenxvec(k)
         duzthenxvecy(k) = duzthenxvecy(k) * czthenxvec(k)
         duzthenxvecz(k) = duzthenxvecz(k) * czthenxvec(k)
         dvzthenxvecx(k) = dvzthenxvecx(k) * czthenxvec(k)
         dvzthenxvecy(k) = dvzthenxvecy(k) * czthenxvec(k)
         dvzthenxvecz(k) = dvzthenxvecz(k) * czthenxvec(k)
      enddo

c
c     force distribution for the Bisector local coordinate method
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dubisectorvecx(k) =           uvvecx(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * uwvecx(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dubisectorvecy(k) =           uvvecy(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * uwvecy(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dubisectorvecz(k) =           uvvecz(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * uwvecz(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dvbisectorvecx(k) = -         uvvecx(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * vwvecx(k)     * dphidwvec(k)
     &                               * invvsizvec(k)
         dvbisectorvecy(k) = -         uvvecy(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * vwvecy(k)     * dphidwvec(k)
     &                               * invvsizvec(k)
         dvbisectorvecz(k) = -         uvvecz(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5_ti_p * vwvecz(k)     * dphidwvec(k)
     &                               * invvsizvec(k)

         dubisectorvecx(k) = dubisectorvecx(k) * cbisectorvec(k)
         dubisectorvecy(k) = dubisectorvecy(k) * cbisectorvec(k)
         dubisectorvecz(k) = dubisectorvecz(k) * cbisectorvec(k)
         dvbisectorvecx(k) = dvbisectorvecx(k) * cbisectorvec(k)
         dvbisectorvecy(k) = dvbisectorvecy(k) * cbisectorvec(k)
         dvbisectorvecz(k) = dvbisectorvecz(k) * cbisectorvec(k)
       enddo
c   
c     force distribution for the Z-Bisect local coordinate method
c   
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         invursinvec(k) = (1.0_ti_p - urcosvec(k)**2) ** ( - half )
         vssinvec(k)    = (1.0_ti_p - vscosvec(k)**2) **   half
         wssinvec(k)    = (1.0_ti_p - wscosvec(k)**2) **   half
         duzbisectvecx(k) =  urvecx(k) * dphidrvec(k) * invusizvec(k)
     &                                 * invursinvec(k)
     &                    +  usvecx(k) * dphidsvec(k) * invusizvec(k)
         duzbisectvecy(k) =  urvecy(k) * dphidrvec(k) * invusizvec(k)
     &                                 * invursinvec(k)
     &                    +  usvecy(k) * dphidsvec(k) * invusizvec(k)
         duzbisectvecz(k) =  urvecz(k) * dphidrvec(k) * invusizvec(k)
     &                                 * invursinvec(k)
     &                    +  usvecz(k) * dphidsvec(k) * invusizvec(k)
       enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dvzbisectvecx(k) =  (  vssinvec(k) * svecx(k)
     &                        - vscosvec(k) * t1vecx(k))
     &                     * dphiduvec(k) * invvsizvec(k)
     &                     * invutsinvec(k)
         dvzbisectvecy(k) =  (  vssinvec(k) * svecy(k)
     &                        - vscosvec(k) * t1vecy(k))
     &                     * dphiduvec(k) * invvsizvec(k)
     &                     * invutsinvec(k)
         dvzbisectvecz(k) =  (  vssinvec(k) * svecz(k)
     &                        - vscosvec(k) * t1vecz(k))
     &                     * dphiduvec(k) * invvsizvec(k)
     &                     * invutsinvec(k)
       enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dwzbisectvecx(k) =  (  wssinvec(k) * svecx(k)
     &                        - wscosvec(k) * t2vecx(k))
     &                     * dphiduvec(k) * invwsizvec(k)
     &                     * invutsinvec(k)
         dwzbisectvecy(k) =  (  wssinvec(k) * svecy(k)
     &                        - wscosvec(k) * t2vecy(k))
     &                     * dphiduvec(k) * invwsizvec(k)
     &                     * invutsinvec(k)
         dwzbisectvecz(k) =  (  wssinvec(k) * svecz(k)
     &                        - wscosvec(k) * t2vecz(k))
     &                     * dphiduvec(k) * invwsizvec(k)
     &                     * invutsinvec(k)
       enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duzbisectvecx(k) = duzbisectvecx(k) *  czbisectvec(k)
         duzbisectvecy(k) = duzbisectvecy(k) *  czbisectvec(k)
         duzbisectvecz(k) = duzbisectvecz(k) *  czbisectvec(k)
         dvzbisectvecx(k) = dvzbisectvecx(k) *  czbisectvec(k)
         dvzbisectvecy(k) = dvzbisectvecy(k) *  czbisectvec(k)
         dvzbisectvecz(k) = dvzbisectvecz(k) *  czbisectvec(k)
         dwzbisectvecx(k) = dwzbisectvecx(k) *  czbisectvec(k)
         dwzbisectvecy(k) = dwzbisectvecy(k) *  czbisectvec(k)
         dwzbisectvecz(k) = dwzbisectvecz(k) *  czbisectvec(k)
      enddo

c
c     force distribution for the 3-Fold local coordinate method
c

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         pvecx(k)        =  uvecx(k) + vvecx(k) + wvecx(k)
         pvecy(k)        =  uvecy(k) + vvecy(k) + wvecy(k)
         pvecz(k)        =  uvecz(k) + vvecz(k) + wvecz(k)
         invpsizvec(k)   =  (  pvecx(k)**2 + pvecy(k)**2 + pvecz(k)**2
     &                      ) ** ( - half )

         pvecx(k)        =  pvecx(k) * invpsizvec(k)
         pvecy(k)        =  pvecy(k) * invpsizvec(k)
         pvecz(k)        =  pvecz(k) * invpsizvec(k)
         upcosvec(k)     =  uvecx(k) * pvecx(k) + uvecy(k) * pvecy(k)
     &                    + uvecz(k) * pvecz(k)
         vpcosvec(k)     =  vvecx(k) * pvecx(k) + vvecy(k) * pvecy(k)
     &                    + vvecz(k) * pvecz(k)
         wpcosvec(k)     =  wvecx(k) * pvecx(k) + wvecy(k) * pvecy(k)
     &                    + wvecz(k) * pvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)        =  uvecx(k) + vvecx(k)
         rvecy(k)        =  uvecy(k) + vvecy(k)
         rvecz(k)        =  uvecz(k) + vvecz(k)
         invrsizvec(k)   =  (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                      ) ** ( - half )

         rvecx(k)        =  rvecx(k) * invrsizvec(k)
         rvecy(k)        =  rvecy(k) * invrsizvec(k)
         rvecz(k)        =  rvecz(k) * invrsizvec(k)
         rwcosvec(k)     =  rvecx(k) * wvecx(k) + rvecy(k) * wvecy(k)
     &                    + rvecz(k) * wvecz(k)
         invrwsinvec(k)  = (1.0_ti_p - rwcosvec(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphidrvec(k)    = - trqvecx(k) * rvecx(k)
     &                     - trqvecy(k) * rvecy(k)
     &                     - trqvecz(k) * rvecz(k)
         delvecx(k)      =   rvecy(k) * wvecz(k) - rvecz(k) * wvecy(k)
         delvecy(k)      =   rvecz(k) * wvecx(k) - rvecx(k) * wvecz(k)
         delvecz(k)      =   rvecx(k) * wvecy(k) - rvecy(k) * wvecx(k)
         invdelsizvec(k) =   (  delvecx(k)**2 + delvecy(k)**2
     &                        + delvecz(k)**2 ) ** ( - half )
         delvecx(k)      =   delvecx(k) * invdelsizvec(k)
         delvecy(k)      =   delvecy(k) * invdelsizvec(k)
         delvecz(k)      =   delvecz(k) * invdelsizvec(k)
         dphiddelvec(k)  = - trqvecx(k) * delvecx(k)
     &                     - trqvecy(k) * delvecy(k)
     &                     - trqvecz(k) * delvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         epsvecx(k)     = delvecy(k) * wvecz(k) - delvecz(k) * wvecy(k)
         epsvecy(k)     = delvecz(k) * wvecx(k) - delvecx(k) * wvecz(k)
         epsvecz(k)     = delvecx(k) * wvecy(k) - delvecy(k) * wvecx(k)
         dw3foldvecx(k) =   delvecx(k)    * dphidrvec(k)
     &                    * invwsizvec(k) * invrwsinvec(k)
     &                   +  epsvecx(k)    * dphiddelvec(k) * wpcosvec(k)
     &                    * invwsizvec(k) * invpsizvec(k) 
         dw3foldvecy(k) =   delvecy(k)    * dphidrvec  (k)
     &                    * invwsizvec(k) * invrwsinvec(k)
     &                   +  epsvecy(k)    * dphiddelvec(k) * wpcosvec(k)
     &                    * invwsizvec(k) * invpsizvec(k) 
         dw3foldvecz(k) =   delvecz(k)    * dphidrvec(k)  
     &                    * invwsizvec(k) * invrwsinvec(k)
     &                   +  epsvecz(k)    * dphiddelvec(k) * wpcosvec(k)
     &                    * invwsizvec(k) * invpsizvec(k) 
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)       =    vvecx(k) + wvecx(k)
         rvecy(k)       =    vvecy(k) + wvecy(k)
         rvecz(k)       =    vvecz(k) + wvecz(k)
         invrsizvec(k)  =   (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                      ) ** ( - half )
 
         rvecx(k)       =   rvecx(k) * invrsizvec(k)
         rvecy(k)       =   rvecy(k) * invrsizvec(k)
         rvecz(k)       =   rvecz(k) * invrsizvec(k)
         rucosvec(k)    =   rvecx(k) * uvecx(k) + rvecy(k) * uvecy(k)
     &                    + rvecz(k) * uvecz(k)
         invrusinvec(k) = (1.0_ti_p - rucosvec(k)**2) ** ( - half )

         dphidrvec(k)   = - trqvecx(k) * rvecx(k)
     &                    - trqvecy(k) * rvecy(k)
     &                    - trqvecz(k) * rvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop

         delvecx(k)      =   rvecy(k) * uvecz(k) - rvecz(k) * uvecy(k)
         delvecy(k)      =   rvecz(k) * uvecx(k) - rvecx(k) * uvecz(k)
         delvecz(k)      =   rvecx(k) * uvecy(k) - rvecy(k) * uvecx(k)
         invdelsizvec(k) =  (  delvecx(k)**2 + delvecy(k)**2
     &                       + delvecz(k)**2 ) ** ( - half )
         delvecx(k)      =   delvecx(k) * invdelsizvec(k)
         delvecy(k)      =   delvecy(k) * invdelsizvec(k)
         delvecz(k)      =   delvecz(k) * invdelsizvec(k)
         dphiddelvec(k)  = - trqvecx(k) * delvecx(k)
     &                     - trqvecy(k) * delvecy(k)
     &                     - trqvecz(k) * delvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         epsvecx(k)     = delvecy(k) * uvecz(k) - delvecz(k) * uvecy(k)
         epsvecy(k)     = delvecz(k) * uvecx(k) - delvecx(k) * uvecz(k)
         epsvecz(k)     = delvecx(k) * uvecy(k) - delvecy(k) * uvecx(k)

         du3foldvecx(k) =   delvecx(k)    * dphidrvec(k)
     &                    * invusizvec(k) * invrusinvec(k)
     &                   +  epsvecx(k)    * dphiddelvec(k) * upcosvec(k)
     &                    * invusizvec(k) * invpsizvec(k) 
         du3foldvecy(k) =   delvecy(k)    * dphidrvec(k)
     &                    * invusizvec(k) * invrusinvec(k)
     &                   +  epsvecy(k)    * dphiddelvec(k) * upcosvec(k)
     &                    * invusizvec(k) * invpsizvec(k)
         du3foldvecz(k) =   delvecz(k)    * dphidrvec(k)
     &                    * invusizvec(k) * invrusinvec(k)
     &                   +  epsvecz(k)    * dphiddelvec(k) * upcosvec(k)
     &                    * invusizvec(k) * invpsizvec(k) 

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         rvecx(k)       =  uvecx(k) + wvecx(k)
         rvecy(k)       =  uvecy(k) + wvecy(k)
         rvecz(k)       =  uvecz(k) + wvecz(k)
         invrsizvec(k)  = (  rvecx(k)**2 + rvecy(k)**2 + rvecz(k)**2
     &                    ) ** ( - half )

         rvecx(k)       =  rvecx(k) * invrsizvec(k)
         rvecy(k)       =  rvecy(k) * invrsizvec(k)
         rvecz(k)       =  rvecz(k) * invrsizvec(k)
         rvcosvec(k)    =  rvecx(k) * vvecx(k) + rvecy(k) * vvecy(k)
     &                   + rvecz(k) * vvecz(k)
         invrvsinvec(k) = (1.0_ti_p - rvcosvec(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dphidrvec(k)    = - trqvecx(k) * rvecx(k)
     &                     - trqvecy(k) * rvecy(k)
     &                     - trqvecz(k) * rvecz(k)

         delvecx(k)      =  rvecy(k) * vvecz(k) - rvecz(k) * vvecy(k)
         delvecy(k)      =  rvecz(k) * vvecx(k) - rvecx(k) * vvecz(k)
         delvecz(k)      =  rvecx(k) * vvecy(k) - rvecy(k) * vvecx(k)

         invdelsizvec(k) = (  delvecx(k)**2 + delvecy(k)**2
     &                      + delvecz(k)**2 ) ** ( - half )
         delvecx(k)      =  delvecx(k) * invdelsizvec(k)
         delvecy(k)      =  delvecy(k) * invdelsizvec(k)
         delvecz(k)      =  delvecz(k) * invdelsizvec(k)

         dphiddelvec(k)  = - trqvecx(k) * delvecx(k)
     &                     - trqvecy(k) * delvecy(k)
     &                     - trqvecz(k) * delvecz(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         epsvecx(k)     = delvecy(k) * vvecz(k) - delvecz(k) * vvecy(k)
         epsvecy(k)     = delvecz(k) * vvecx(k) - delvecx(k) * vvecz(k)
         epsvecz(k)     = delvecx(k) * vvecy(k) - delvecy(k) * vvecx(k)

         dv3foldvecx(k) =   delvecx(k)    * dphidrvec(k)
     &                    * invvsizvec(k) * invrvsinvec(k)
     &                   +  epsvecx(k)    * dphiddelvec(k) * vpcosvec(k)
     &                    * invvsizvec(k) * invpsizvec(k)
         dv3foldvecy(k) =   delvecy(k)    * dphidrvec(k)
     &                    * invvsizvec(k) * invrvsinvec(k)
     &                   +  epsvecy(k)    * dphiddelvec(k) * vpcosvec(k)
     &                    * invvsizvec(k) * invpsizvec(k)
         dv3foldvecz(k) =   delvecz(k)    * dphidrvec(k)
     &                    * invvsizvec(k) * invrvsinvec(k)
     &                   +  epsvecz(k)    * dphiddelvec(k) * vpcosvec(k)
     &                    * invvsizvec(k) * invpsizvec(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         du3foldvecx(k) = du3foldvecx(k) * c3foldvec(k)
         du3foldvecy(k) = du3foldvecy(k) * c3foldvec(k)
         du3foldvecz(k) = du3foldvecz(k) * c3foldvec(k)
         dv3foldvecx(k) = dv3foldvecx(k) * c3foldvec(k)
         dv3foldvecy(k) = dv3foldvecy(k) * c3foldvec(k)
         dv3foldvecz(k) = dv3foldvecz(k) * c3foldvec(k)
         dw3foldvecx(k) = dw3foldvecx(k) * c3foldvec(k)
         dw3foldvecy(k) = dw3foldvecy(k) * c3foldvec(k)
         dw3foldvecz(k) = dw3foldvecz(k) * c3foldvec(k)
      enddo
!DIR$ ASSUME_ALIGNED duvecx : 64
!DIR$ ASSUME_ALIGNED duvecy : 64
!DIR$ ASSUME_ALIGNED duvecz : 64
!DIR$ ASSUME_ALIGNED dvvecx : 64
!DIR$ ASSUME_ALIGNED dvvecy : 64
!DIR$ ASSUME_ALIGNED dvvecz : 64
!DIR$ ASSUME_ALIGNED dwvecx : 64
!DIR$ ASSUME_ALIGNED dwvecy : 64
!DIR$ ASSUME_ALIGNED dwvecz : 64
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         duvecx(k) =  duzonlyvecx(k)    + duzthenxvecx(k) 
     &              + dubisectorvecx(k) + duzbisectvecx(k) 
     &              + du3foldvecx(k)   

         duvecy(k) =  duzonlyvecy(k)    + duzthenxvecy(k)
     &              + dubisectorvecy(k) + duzbisectvecy(k) 
     &              + du3foldvecy(k)   

         duvecz(k) =  duzonlyvecz(k)    + duzthenxvecz(k)
     &              + dubisectorvecz(k) + duzbisectvecz(k)
     &              + du3foldvecz(k) 
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dvvecx(k) =  dvzonlyvecx(k)    + dvzthenxvecx(k)   
     &              + dvbisectorvecx(k) + dvzbisectvecx(k)  
     &              + dv3foldvecx(k)    

         dvvecy(k) =  dvzonlyvecy(k)    + dvzthenxvecy(k)   
     &              + dvbisectorvecy(k) + dvzbisectvecy(k)  
     &           + dv3foldvecy(k)    

         dvvecz(k) =  dvzonlyvecz(k)    + dvzthenxvecz(k)   
     &              + dvbisectorvecz(k) + dvzbisectvecz(k)  
     &              + dv3foldvecz(k)    
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         dwvecx(k) =  dwzonlyvecx(k)    + dwzthenxvecx(k)   
     &              + dwbisectorvecx(k) + dwzbisectvecx(k)  
     &              + dw3foldvecx(k)    

         dwvecy(k) =  dwzonlyvecy(k)    + dwzthenxvecy(k)   
     &              + dwbisectorvecy(k) + dwzbisectvecy(k)  
     &              + dw3foldvecy(k)    

         dwvecz(k) =  dwzonlyvecz(k)    + dwzthenxvecz(k)   
     &              + dwbisectorvecz(k) + dwzbisectvecz(k)  
     &              + dw3foldvecz(k)    
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(ialocvec(k)) =  devecx(ialocvec(k)) + duvecx(k)
         devecy(ialocvec(k)) =  devecy(ialocvec(k)) + duvecy(k)
         devecz(ialocvec(k)) =  devecz(ialocvec(k)) + duvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(iblocvec(k)) =  devecx(iblocvec(k))
     &                          - duvecx(k) -dvvecx(k) -dwvecx(k)
         devecy(iblocvec(k)) =  devecy(iblocvec(k))
     &                          - duvecy(k) -dvvecy(k) -dwvecy(k)
         devecz(iblocvec(k)) =  devecz(iblocvec(k))
     &                          - duvecz(k) -dvvecz(k) -dwvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(iclocvec(k)) =  devecx(iclocvec(k)) + dvvecx(k)
         devecy(iclocvec(k)) =  devecy(iclocvec(k)) + dvvecy(k)
         devecz(iclocvec(k)) =  devecz(iclocvec(k)) + dvvecz(k)

      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ SIMD
      do k = 1, npolelocnlloop
         devecx(idlocvec(k)) =  devecx(idlocvec(k)) + dwvecx(k)
         devecy(idlocvec(k)) =  devecy(idlocvec(k)) + dwvecy(k)
         devecz(idlocvec(k)) =  devecz(idlocvec(k)) + dwvecz(k)
      enddo
      call timer_exit( timer_torque )
      end
