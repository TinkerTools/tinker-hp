c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine torquevec  --  convert single site torque to force  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "torque2vec" takes the torque values on multiple sites defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     npolelocnl and npolelocnlloop came from mpole module
c     npoleloc and npolelocloop came from mpole module
c
c     They are transmitted to torquevec2 in nt and ntloop
c
c     INPUT  : ivec is indices of sites 
c              local is loc or locrec1
c              trqvec are trqs of sites 
c              devecx,devecy and devecz are deps coming from epolar
c              or empole
c
c     OUTPUT : dvvecx, dwvecx and duvecx, etc... are forces 
c              devecx,devecy and devecz are deps sent back to epolar
c              or empole
      subroutine torquevec2( nt,ntloop,ivec,local,
     &                       trqvecx,trqvecy,trqvecz,
     &                       dvvecx,dvvecy,dvvecz,
     &                       dwvecx,dwvecy,dwvecz,
     &                       duvecx,duvecy,duvecz,
     &                       devecx,devecy,devecz)
      use atoms
      use deriv
      use domdec
      use mpole
      use sizes
      use mpi
      implicit none
      integer i,j,k,kk
      integer nt1
      integer nt
      integer ntloop
      integer, intent(IN) ::  ivec(*)
      integer, intent(IN) ::  local(*)
      real*8 sqrt3over2,time0,time1
      real*8 half
      real*8 trqvecx(*)
      real*8 trqvecy(*)
      real*8 trqvecz(*)
      real*8 dvvecx(*)
      real*8 dvvecy(*)
      real*8 dvvecz(*)
      real*8 dwvecx(*)
      real*8 dwvecy(*)
      real*8 dwvecz(*)
      real*8 duvecx(*)
      real*8 duvecy(*)
      real*8 duvecz(*)
      real*8 devecx(*)
      real*8 devecy(*)
      real*8 devecz(*)
!DIR$ ATTRIBUTES ALIGN:64::axetypvec
      integer axetypvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec1
      integer ivec1(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::iavec,ibvec
      integer iavec(ntloop),ibvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::icvec,idvec
      integer icvec(ntloop),idvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::ialocvec,iblocvec
      integer ialocvec(ntloop),iblocvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::iclocvec,idlocvec
      integer iclocvec(ntloop),idlocvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: delvecx,delvecy
      real*8 delvecx(ntloop)
      real*8 delvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: delvecz
      real*8 delvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::epsvecx,epsvecy,epsvecz
      real*8 epsvecx(ntloop)
      real*8 epsvecy(ntloop)
      real*8 epsvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invdelsizvec,epssizvec
      real*8 invdelsizvec(ntloop)
      real*8 epssizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::pvecx,pvecy,pvecz,invpsizvec
      real*8 pvecx(ntloop)
      real*8 pvecy(ntloop)
      real*8 pvecz(ntloop)
      real*8 invpsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::rvecx,rvecy,rvecz
      real*8 rvecx(ntloop)
      real*8 rvecy(ntloop)
      real*8 rvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invrsizvec
      real*8 invrsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::svecx,svecy,svecz
      real*8 svecx(ntloop)
      real*8 svecy(ntloop)
      real*8 svecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invssizvec
      real*8 invssizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::invusizvec
      real*8 invusizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::uvecx,uvecy,uvecz
      real*8 uvecx(ntloop)
      real*8 uvecy(ntloop)
      real*8 uvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::invvsizvec
      real*8 invvsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::vvecx,vvecy,vvecz
      real*8 vvecx(ntloop)
      real*8 vvecy(ntloop)
      real*8 vvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::invwsizvec
      real*8 invwsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::wvecx,wvecy,wvecz
      real*8 wvecx(ntloop)
      real*8 wvecy(ntloop)
      real*8 wvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::uvvecx,uvvecy,uvvecz
      real*8 uvvecx(ntloop)
      real*8 uvvecy(ntloop)
      real*8 uvvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invuvsizvec
      real*8 invuvsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::uwvecx,uwvecy,uwvecz
      real*8 uwvecx(ntloop)
      real*8 uwvecy(ntloop)
      real*8 uwvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invuwsizvec
      real*8 invuwsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::vwvecx,vwvecy,vwvecz
      real*8 vwvecx(ntloop)
      real*8 vwvecy(ntloop)
      real*8 vwvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invvwsizvec
      real*8 invvwsizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::usvecx,usvecy,usvecz
      real*8 usvecx(ntloop)
      real*8 usvecy(ntloop)
      real*8 usvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invussizvec
      real*8 invussizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::vsvecx,vsvecy,vsvecz
      real*8 vsvecx(ntloop)
      real*8 vsvecy(ntloop)
      real*8 vsvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invvssizvec
      real*8 invvssizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::wsvecx,wsvecy,wsvecz
      real*8 wsvecx(ntloop)
      real*8 wsvecy(ntloop)
      real*8 wsvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invwssizvec
      real*8 invwssizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::urvecx,urvecy,urvecz
      real*8 urvecx(ntloop)
      real*8 urvecy(ntloop)
      real*8 urvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invursizvec
      real*8 invursizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::t1vecx,t1vecy,t1vecz
      real*8 t1vecx(ntloop)
      real*8 t1vecy(ntloop)
      real*8 t1vecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invt1sizvec
      real*8 invt1sizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::t2vecx,t2vecy,t2vecz
      real*8 t2vecx(ntloop)
      real*8 t2vecy(ntloop)
      real*8 t2vecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: invt2sizvec
      real*8 invt2sizvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::rucosvec,invrusinvec
      real*8 rucosvec(ntloop)
      real*8 invrusinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::rvcosvec,invrvsinvec
      real*8 rvcosvec(ntloop)
      real*8 invrvsinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::rwcosvec,invrwsinvec
      real*8 rwcosvec(ntloop)
      real*8 invrwsinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::uvcosvec,invuvsinvec
      real*8 uvcosvec(ntloop)
      real*8 invuvsinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::uwcosvec,uwsinvec
      real*8 uwcosvec(ntloop)
      real*8 uwsinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::vwcosvec,vwsinvec
      real*8 vwcosvec(ntloop)
      real*8 vwsinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::upcosvec
      real*8 upcosvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::urcosvec,invursinvec
      real*8 urcosvec(ntloop)
      real*8 invursinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::uscosvec
      real*8 uscosvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::vpcosvec
      real*8 vpcosvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::vscosvec,vssinvec
      real*8 vscosvec(ntloop)
      real*8 vssinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::wpcosvec
      real*8 wpcosvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::wscosvec,wssinvec
      real*8 wscosvec(ntloop)
      real*8 wssinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::ut1cosvec,ut1sinvec
      real*8 ut1cosvec(ntloop)
      real*8 ut1sinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::ut2cosvec,ut2sinvec,invutsinvec
      real*8 ut2cosvec(ntloop)
      real*8 ut2sinvec(ntloop)
      real*8 invutsinvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::dphiduvec,dphidvvec
      real*8 dphiduvec(ntloop)
      real*8 dphidvvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::dphidwvec,dphidrvec
      real*8 dphidwvec(ntloop)
      real*8 dphidrvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::dphiddelvec
      real*8 dphiddelvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::dphidsvec
      real*8 dphidsvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::xiavec,yiavec,ziavec
      real*8 xiavec(ntloop)
      real*8 yiavec(ntloop)
      real*8 ziavec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::xibvec,yibvec,zibvec
      real*8 xibvec(ntloop)
      real*8 yibvec(ntloop)
      real*8 zibvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::xicvec,yicvec,zicvec
      real*8 xicvec(ntloop)
      real*8 yicvec(ntloop)
      real*8 zicvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64::xidvec,yidvec,zidvec
      real*8 xidvec(ntloop)
      real*8 yidvec(ntloop)
      real*8 zidvec(ntloop)
c     3-Fold
!DIR$ ATTRIBUTES ALIGN:64:: c3foldvec
      real*8 c3foldvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecx
      real*8 du3foldvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecy
      real*8 du3foldvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: du3foldvecz
      real*8 du3foldvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecx
      real*8 dv3foldvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecy
      real*8 dv3foldvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dv3foldvecz
      real*8 dv3foldvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecx
      real*8 dw3foldvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecy
      real*8 dw3foldvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dw3foldvecz
      real*8 dw3foldvecz(ntloop)
c     Bisector
!DIR$ ATTRIBUTES ALIGN:64:: cbisectorvec
      real*8 cbisectorvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecx
      real*8 dubisectorvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecy
      real*8 dubisectorvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dubisectorvecz
      real*8 dubisectorvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecx
      real*8 dvbisectorvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecy
      real*8 dvbisectorvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvbisectorvecz
      real*8 dvbisectorvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecx
      real*8 dwbisectorvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecy
      real*8 dwbisectorvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwbisectorvecz
      real*8 dwbisectorvecz(ntloop)
c     Z-Bisect
!DIR$ ATTRIBUTES ALIGN:64:: czbisectvec
      real*8 czbisectvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecx
      real*8 duzbisectvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecy
      real*8 duzbisectvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzbisectvecz
      real*8 duzbisectvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecx
      real*8 dvzbisectvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecy
      real*8 dvzbisectvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzbisectvecz
      real*8 dvzbisectvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecx
      real*8 dwzbisectvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecy
      real*8 dwzbisectvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzbisectvecz
      real*8 dwzbisectvecz(ntloop)
C     Z-Only
!DIR$ ATTRIBUTES ALIGN:64:: czonlyvec
      real*8 czonlyvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecx
      real*8 duzonlyvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecy
      real*8 duzonlyvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzonlyvecz
      real*8 duzonlyvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecx
      real*8 dvzonlyvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecy
      real*8 dvzonlyvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzonlyvecz
      real*8 dvzonlyvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecx
      real*8 dwzonlyvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecy
      real*8 dwzonlyvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzonlyvecz
      real*8 dwzonlyvecz(ntloop)
C     Z-then-X
!DIR$ ATTRIBUTES ALIGN:64:: czthenxvec
      real*8 czthenxvec(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecx
      real*8 duzthenxvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecy
      real*8 duzthenxvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: duzthenxvecz
      real*8 duzthenxvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecx
      real*8 dvzthenxvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecy
      real*8 dvzthenxvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dvzthenxvecz
      real*8 dvzthenxvecz(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecx
      real*8 dwzthenxvecx(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecy
      real*8 dwzthenxvecy(ntloop)
!DIR$ ATTRIBUTES ALIGN:64:: dwzthenxvecz
      real*8 dwzthenxvecz(ntloop)
      character*8 paxe(ntloop)

c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'torquevec2'

      time0 = mpi_wtime()
      sqrt3over2 = sqrt(3.0d0)/2.0d0
      half       = 0.5d0
c
c     zero out force components on local frame-defining atoms
c
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         iavec(k)     = -1
         ibvec(k)     = -1
         icvec(k)     = -1
         idvec(k)     = -1
         ialocvec(k)  = 1
         iblocvec(k)  = 1
         iclocvec(k)  = 1
         idlocvec(k)  = 1
         ivec1(k)     = 1
         axetypvec(k) = 0
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         c3foldvec(k)   = 0.0d0
         du3foldvecx(k) = 0.0d0
         du3foldvecy(k) = 0.0d0
         du3foldvecz(k) = 0.0d0
         dv3foldvecx(k) = 0.0d0
         dv3foldvecy(k) = 0.0d0
         dv3foldvecz(k) = 0.0d0
         dw3foldvecx(k) = 0.0d0
         dw3foldvecy(k) = 0.0d0
         dw3foldvecz(k) = 0.0d0
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         cbisectorvec(k)   = 0.0d0
         dubisectorvecx(k) = 0.0d0
         dubisectorvecy(k) = 0.0d0
         dubisectorvecz(k) = 0.0d0
         dvbisectorvecx(k) = 0.0d0
         dvbisectorvecy(k) = 0.0d0
         dvbisectorvecz(k) = 0.0d0
         dwbisectorvecx(k) = 0.0d0
         dwbisectorvecy(k) = 0.0d0
         dwbisectorvecz(k) = 0.0d0
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         czbisectvec(k)   = 0.0d0
         duzbisectvecx(k) = 0.0d0
         duzbisectvecy(k) = 0.0d0
         duzbisectvecz(k) = 0.0d0
         dvzbisectvecx(k) = 0.0d0
         dvzbisectvecy(k) = 0.0d0
         dvzbisectvecz(k) = 0.0d0
         dwzbisectvecx(k) = 0.0d0
         dwzbisectvecy(k) = 0.0d0
         dwzbisectvecz(k) = 0.0d0
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         czonlyvec(k)   = 0.0d0
         duzonlyvecx(k) = 0.0d0
         duzonlyvecy(k) = 0.0d0
         duzonlyvecz(k) = 0.0d0
         dvzonlyvecx(k) = 0.0d0
         dvzonlyvecy(k) = 0.0d0
         dvzonlyvecz(k) = 0.0d0
         dwzonlyvecx(k) = 0.0d0
         dwzonlyvecy(k) = 0.0d0
         dwzonlyvecz(k) = 0.0d0
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         czthenxvec(k)   = 0.0d0
         duzthenxvecx(k) = 0.0d0
         duzthenxvecy(k) = 0.0d0
         duzthenxvecz(k) = 0.0d0
         dvzthenxvecx(k) = 0.0d0
         dvzthenxvecy(k) = 0.0d0
         dvzthenxvecz(k) = 0.0d0
         dwzthenxvecx(k) = 0.0d0
         dwzthenxvecy(k) = 0.0d0
         dwzthenxvecz(k) = 0.0d0
      enddo
c
c     get the local frame type and the frame-defining atoms
c
      do k = 1, ntloop
         if(k.le.nt) then
            select case (polaxe(ivec(k)))
                   case ('3-Fold'  ) ;  axetypvec(k) = 1
                                        c3foldvec(k)= 1.0d0
                   case ('Bisector') ;  axetypvec(k) = 2
                                        cbisectorvec(k)= 1.0d0
                   case ('Z-Bisect') ;  axetypvec(k) = 3
                                        czbisectvec(k)= 1.0d0
                   case ('Z-Only'  ) ;  axetypvec(k) = 4
                                        czonlyvec(k)= 1.0d0
                   case ('Z-then-X') ;  axetypvec(k) = 5
                                        czthenxvec(k)= 1.0d0
                   case ('None'    ) ;  axetypvec(k) = 0
                   case default      ;  axetypvec(k) = 0
            endselect
         else
            axetypvec(k) = 0
         endif
      enddo

      kk = 0
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1,ntloop
         if (axetypvec(k).ne.0) then
            kk = kk +1
            ivec1 (kk) = ivec (k)
         endif
      enddo
      nt1 = kk
      if (nt1.eq.0) return ! All axetyp are 0, so return
    
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, ntloop
         iavec(k) = zaxis(ivec1(k))
         ibvec(k) = ipole(ivec1(k))
         icvec(k) = xaxis(ivec1(k))
         idvec(k) = yaxis(ivec1(k))
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, ntloop
         if(iavec(k).gt.0) ialocvec(k) = local(iavec(k))
         if(ibvec(k).gt.0) iblocvec(k) = local(ibvec(k))
         if(icvec(k).gt.0) iclocvec(k) = local(icvec(k))
         if(idvec(k).gt.0) idlocvec(k) = local(idvec(k))
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, ntloop
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, ntloop
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         uvecx(k)  = xiavec(k)- xibvec(k)
         uvecy(k)  = yiavec(k)- yibvec(k)
         uvecz(k)  = ziavec(k)- zibvec(k)
         invusizvec(k) = (  uvecx(k)**2 + uvecy(k)**2 + uvecz(k)**2
     &                   ) ** ( - half )
         uvecx(k) = uvecx(k) * invusizvec(k)
         uvecy(k) = uvecy(k) * invusizvec(k)
         uvecz(k) = uvecz(k) * invusizvec(k)
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         vvecx(k)  = xicvec(k) - xibvec(k)
         vvecy(k)  = yicvec(k) - yibvec(k)
         vvecz(k)  = zicvec(k) - zibvec(k)
         invvsizvec(k) = (  vvecx(k)**2 + vvecy(k)**2 + vvecz(k)**2
     &                   ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         if (axetypvec(k).eq.4) then ! 'Z-Only' 
            if (abs(uvecx(k)).gt.sqrt3over2) then
               vvecx(k)  = 0.0d0
            else
               vvecx(k)  = 1.0d0
            endif
            vvecy(k) = 1.0d0 - vvecx(k)
            vvecz(k)  = 0.0d0
            invvsizvec(k) = 1.0d0
         endif
      enddo
c
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         rvecx(k) = rvecx(k) * invrsizvec(k)
         svecx(k) = svecx(k) * invssizvec(k)
         rvecy(k) = rvecy(k) * invrsizvec(k)
         svecy(k) = svecy(k) * invssizvec(k)
         rvecz(k) = rvecz(k) * invrsizvec(k)
         svecz(k) = svecz(k) * invssizvec(k)
      enddo
c     find the perpendicular and angle for each pair of axes
c
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         uvvecx(k)  =  vvecy(k) * uvecz(k) - vvecz(k) * uvecy(k)
         uvvecy(k)  =  vvecz(k) * uvecx(k) - vvecx(k) * uvecz(k)
         uvvecz(k)  =  vvecx(k) * uvecy(k) - vvecy(k) * uvecx(k)
         invuvsizvec(k) = (  uvvecx(k)**2 + uvvecy(k)**2 + uvvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         uwvecx(k)  =  wvecy(k) * uvecz(k) - wvecz(k) * uvecy(k)
         uwvecy(k)  =  wvecz(k) * uvecx(k) - wvecx(k) * uvecz(k)
         uwvecz(k)  =  wvecx(k) * uvecy(k) - wvecy(k) * uvecx(k)
         invuwsizvec(k) = (  uwvecx(k)**2 + uwvecy(k)**2 + uwvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         vwvecx(k)  =  wvecy(k) * vvecz(k) - wvecz(k) * vvecy(k)
         vwvecy(k)  =  wvecz(k) * vvecx(k) - wvecx(k) * vvecz(k)
         vwvecz(k)  =  wvecx(k) * vvecy(k) - wvecy(k) * vvecx(k)
         invvwsizvec(k) = (  vwvecx(k)**2 + vwvecy(k)**2 + vwvecz(k)**2
     &                    )** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         urvecx(k)  =  rvecy(k) * uvecz(k) - rvecz(k) * uvecy(k)
         urvecy(k)  =  rvecz(k) * uvecx(k) - rvecx(k) * uvecz(k)
         urvecz(k)  =  rvecx(k) * uvecy(k) - rvecy(k) * uvecx(k)
         invursizvec(k) = (  urvecx(k)**2 + urvecy(k)**2 + urvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         usvecx(k)  =  svecy(k) * uvecz(k) - svecz(k) * uvecy(k)
         usvecy(k)  =  svecz(k) * uvecx(k) - svecx(k) * uvecz(k)
         usvecz(k)  =  svecx(k) * uvecy(k) - svecy(k) * uvecx(k)
         invussizvec(k) = (  usvecx(k)**2 + usvecy(k)**2 + usvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         vsvecx(k)  =  svecy(k) * vvecz(k) - svecz(k) * vvecy(k)
         vsvecy(k)  =  svecz(k) * vvecx(k) - svecx(k) * vvecz(k)
         vsvecz(k)  =  svecx(k) * vvecy(k) - svecy(k) * vvecx(k)
         invvssizvec(k) = (  vsvecx(k)**2 + vsvecy(k)**2 + vsvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, ntloop
         wsvecx(k)  =  svecy(k) * wvecz(k) - svecz(k) * wvecy(k)
         wsvecy(k)  =  svecz(k) * wvecx(k) - svecx(k) * wvecz(k)
         wsvecz(k)  =  svecx(k) * wvecy(k) - svecy(k) * wvecx(k)
         invwssizvec(k) = (  wsvecx(k)**2 + wsvecy(k)**2 + wsvecz(k)**2
     &                    ) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         uvcosvec(k) =  uvecx(k) * vvecx(k) + uvecy(k) * vvecy(k)
     &                + uvecz(k) * vvecz(k)
         uwcosvec(k) =  uvecx(k) * wvecx(k) + uvecy(k) * wvecy(k)
     &                + uvecz(k) * wvecz(k)
         vwcosvec(k) =  vvecx(k) * wvecx(k) + vvecy(k) * wvecy(k)
     &                + vvecz(k) * wvecz(k)
         invuvsinvec(k) = (1.0d0 - uvcosvec(k)**2) ** ( - half )
      enddo
c
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         t1vecx(k) =  vvecx(k) - svecx(k) * vscosvec(k)
         t1vecy(k) =  vvecy(k) - svecy(k) * vscosvec(k)
         t1vecz(k) =  vvecz(k) - svecz(k) * vscosvec(k)
         t2vecx(k) =  wvecx(k) - svecx(k) * wscosvec(k)
         t2vecy(k) =  wvecy(k) - svecy(k) * wscosvec(k)
         t2vecz(k) =  wvecz(k) - svecz(k) * wscosvec(k)
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         ut1cosvec(k) =  uvecx(k) * t1vecx(k) + uvecy(k) * t1vecy(k)
     &                 + uvecz(k) * t1vecz(k)
         ut2cosvec(k) =  uvecx(k) * t2vecx(k) + uvecy(k) * t2vecy(k)
     &                 + uvecz(k) * t2vecz(k)
         invutsinvec(k) = (  (1.0d0 - ut1cosvec(k)**2) ** half
     &                     + (1.0d0 - ut2cosvec(k)**2) ** half
     &                    ) **-1.0d0
      enddo
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         dphiduvec(k) = - trqvecx(k) * uvecx(k) - trqvecy(k) * uvecy(k)
     &                  - trqvecz(k) * uvecz(k)
         dphidvvec(k) = - trqvecx(k) * vvecx(k) - trqvecy(k) * vvecy(k)
     &                  - trqvecz(k) * vvecz(k)
         dphidwvec(k) = - trqvecx(k) * wvecx(k) - trqvecy(k) * wvecy(k)
     &                  - trqvecz(k) * wvecz(k)
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         dphidrvec(k) = - trqvecx(k) * rvecx(k) - trqvecy(k) * rvecy(k)
     &                  - trqvecz(k) * rvecz(k)
         dphidsvec(k) = - trqvecx(k) * svecx(k) - trqvecy(k) * svecy(k)
     &                  - trqvecz(k) * svecz(k)
      enddo
c
c     force distribution for the Z-Only local coordinate method
c
          
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         dubisectorvecx(k) =           uvvecx(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5d0 * uwvecx(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dubisectorvecy(k) =           uvvecy(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5d0 * uwvecy(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dubisectorvecz(k) =           uvvecz(k)     * dphidvvec(k)
     &                               * invusizvec(k) * invuvsinvec(k)
     &                      +  0.5d0 * uwvecz(k)     * dphidwvec(k)
     &                               * invusizvec(k)
         dvbisectorvecx(k) = -         uvvecx(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5d0 * vwvecx(k)     * dphidwvec(k)
     &                               * invvsizvec(k)
         dvbisectorvecy(k) = -         uvvecy(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5d0 * vwvecy(k)     * dphidwvec(k)
     &                               * invvsizvec(k)
         dvbisectorvecz(k) = -         uvvecz(k)     * dphiduvec(k)
     &                               * invvsizvec(k) * invuvsinvec(k)
     &                      +  0.5d0 * vwvecz(k)     * dphidwvec(k)
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
         invursinvec(k) = (1.0d0 - urcosvec(k)**2) ** ( - half )
         vssinvec(k)    = (1.0d0 - vscosvec(k)**2) **   half
         wssinvec(k)    = (1.0d0 - wscosvec(k)**2) **   half
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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

!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
         invrwsinvec(k)  = (1.0d0 - rwcosvec(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
         invrusinvec(k) = (1.0d0 - rucosvec(k)**2) ** ( - half )

         dphidrvec(k)   = - trqvecx(k) * rvecx(k)
     &                    - trqvecy(k) * rvecy(k)
     &                    - trqvecz(k) * rvecz(k)
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop

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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
         invrvsinvec(k) = (1.0d0 - rvcosvec(k)**2) ** ( - half )
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, ntloop
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
!     if(rank.eq.0.and.tinkerdebug) write(*,*)'ntloop,sizeduvec'
!    &                                 ,ntloop,size(duvecx)
!    &                                 ,nblocloop
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
      do k = 1, ntloop
         if(k.le.nt) then
            devecx(ialocvec(k)) =  devecx(ialocvec(k)) + duvecx(k)
            devecy(ialocvec(k)) =  devecy(ialocvec(k)) + duvecy(k)
            devecz(ialocvec(k)) =  devecz(ialocvec(k)) + duvecz(k)
         endif
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!!DIR$ SIMD
      do k = 1, ntloop
         if(k.le.nt) then
            devecx(iblocvec(k)) =   devecx(iblocvec(k))
     &                            - duvecx(k) - dvvecx(k) - dwvecx(k)
            devecy(iblocvec(k)) =   devecy(iblocvec(k))
     &                            - duvecy(k) - dvvecy(k) - dwvecy(k)
            devecz(iblocvec(k)) =   devecz(iblocvec(k))
     &                            - duvecz(k) - dvvecz(k) - dwvecz(k)
         endif
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!!DIR$ SIMD
      do k = 1, ntloop
         if(k.le.nt) then
            devecx(iclocvec(k)) =  devecx(iclocvec(k)) + dvvecx(k)
            devecy(iclocvec(k)) =  devecy(iclocvec(k)) + dvvecy(k)
            devecz(iclocvec(k)) =  devecz(iclocvec(k)) + dvvecz(k)
         endif
      enddo
!DIR$ ASSUME (mod(ntloop,16).eq.0)
!!DIR$ SIMD
      do k = 1, ntloop
         if(k.le.nt) then
            devecx(idlocvec(k)) =  devecx(idlocvec(k)) + dwvecx(k)
            devecy(idlocvec(k)) =  devecy(idlocvec(k)) + dwvecy(k)
            devecz(idlocvec(k)) =  devecz(idlocvec(k)) + dwvecz(k)
         endif
      enddo
      time1 = mpi_wtime()
      if(rank.eq.0.and.tinkertime) write(*,*) 'time torquevec2',
     &                 time1 - time0
      return
      end
