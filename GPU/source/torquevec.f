c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torquevec  --  convert single site torque to force  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torquevec" takes the torque values on multiple sites defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     INPUT  : nploc is number of sites
c              ivec is indices of sites 
c              trqvec are trqs of sites 
c              devec are deps coming from epolar
c
c     OUTPUT : frcxvec, frcyvec and frczvec are forces 
c              devec are deps sent back to epolar
#include "tinker_precision.h"
      subroutine torquevec (nploc,ivec,trqvec,
     &                      frcxvec,frcyvec,frczvec,devec)
      use atoms
      use deriv
      use domdec
      use mpole
      use sizes
      implicit none
      integer i,j,nploc,nploc1,nploop
      integer ivec(*)
      real(t_p) trqvec(nploc,3)
      real(t_p) frcxvec(nploc,3),frcyvec(nploc,3),frczvec(nploc,3)
      real(t_p) devec(3,*)
      !DIR$ ATTRIBUTES ALIGN:64::mask
      logical,allocatable :: mask(:)
      !DIR$ ATTRIBUTES ALIGN:64::ivec1
      integer,allocatable :: ivec1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: frcxvec1,frcyvec1,frczvec1
      real(t_p),allocatable :: frcxvec1(:,:),frcyvec1(:,:),frczvec1(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: devec1
      real(t_p),allocatable :: devec1(:,:)
      !DIR$ ATTRIBUTES ALIGN:64::iavec,ibvec,icvec,idvec
      integer,allocatable :: iavec(:),ibvec(:),icvec(:),idvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ialocvec,iblocvec
      integer,allocatable :: ialocvec(:),iblocvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::iclocvec,idlocvec
      integer,allocatable :: iclocvec(:),idlocvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::delvec,delsizvec
      real(t_p),allocatable :: delvec(:,:),delsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::epsvec,epssizvec
      real(t_p),allocatable :: epsvec(:,:),epssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::pvec,psizvec
      real(t_p),allocatable :: pvec(:,:),psizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rvec,rsizvec
      real(t_p),allocatable :: rvec(:,:),rsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::svec,ssizvec
      real(t_p),allocatable :: svec(:,:),ssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uvec,usizvec
      real(t_p),allocatable :: uvec(:,:),usizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vvec,vsizvec
      real(t_p),allocatable :: vvec(:,:),vsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wvec,wsizvec
      real(t_p),allocatable :: wvec(:,:),wsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uvvec,uvsizvec
      real(t_p),allocatable :: uvvec(:,:),uvsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uwvec,uwsizvec
      real(t_p),allocatable :: uwvec(:,:),uwsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vwvec,vwsizvec
      real(t_p),allocatable :: vwvec(:,:),vwsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::usvec,ussizvec
      real(t_p),allocatable :: usvec(:,:),ussizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vsvec,vssizvec
      real(t_p),allocatable :: vsvec(:,:),vssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wsvec,wssizvec
      real(t_p),allocatable :: wsvec(:,:),wssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::urvec,ursizvec
      real(t_p),allocatable :: urvec(:,:),ursizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::t1vec,t1sizvec
      real(t_p),allocatable :: t1vec(:,:),t1sizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::t2vec,t2sizvec
      real(t_p),allocatable :: t2vec(:,:),t2sizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rucosvec,rusinvec
      real(t_p),allocatable :: rucosvec(:),rusinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rvcosvec,rvsinvec
      real(t_p),allocatable :: rvcosvec(:),rvsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rwcosvec,rwsinvec
      real(t_p),allocatable :: rwcosvec(:),rwsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uvcosvec,uvsinvec
      real(t_p),allocatable :: uvcosvec(:),uvsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uwcosvec,uwsinvec
      real(t_p),allocatable :: uwcosvec(:),uwsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vwcosvec,vwsinvec
      real(t_p),allocatable :: vwcosvec(:),vwsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::upcosvec,upsinvec
      real(t_p),allocatable :: upcosvec(:),upsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::urcosvec,ursinvec
      real(t_p),allocatable :: urcosvec(:),ursinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uscosvec,ussinvec
      real(t_p),allocatable :: uscosvec(:),ussinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vpcosvec,vpsinvec
      real(t_p),allocatable :: vpcosvec(:),vpsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vscosvec,vssinvec
      real(t_p),allocatable :: vscosvec(:),vssinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wpcosvec,wpsinvec
      real(t_p),allocatable :: wpcosvec(:),wpsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wscosvec,wssinvec
      real(t_p),allocatable :: wscosvec(:),wssinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ut1cosvec,ut1sinvec
      real(t_p),allocatable :: ut1cosvec(:),ut1sinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ut2cosvec,ut2sinvec
      real(t_p),allocatable :: ut2cosvec(:),ut2sinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphiduvec,dphidvvec
      real(t_p),allocatable :: dphiduvec(:),dphidvvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphidwvec,dphidrvec
      real(t_p),allocatable :: dphidwvec(:),dphidrvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphiddelvec
      real(t_p),allocatable :: dphiddelvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphidsvec
      real(t_p),allocatable :: dphidsvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::duvec,dvvec,dwvec
      real(t_p),allocatable :: duvec(:,:),dvvec(:,:),dwvec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64::axetypvec
      character*8,allocatable :: axetypvec(:)
      
      allocate (ivec1(nploc))
      allocate (mask(nploc))
      allocate (iavec(nploc))
      allocate (ibvec(nploc))
      allocate (icvec(nploc))
      allocate (idvec(nploc))
      allocate (ialocvec(nploc))
      allocate (iblocvec(nploc))
      allocate (iclocvec(nploc))
      allocate (idlocvec(nploc))
      allocate (frcxvec1(nploc,3))
      allocate (frcyvec1(nploc,3))
      allocate (frczvec1(nploc,3))
      allocate (rvec(nploc,3))
      allocate (rsizvec(nploc))
      allocate (svec(nploc,3))
      allocate (ssizvec(nploc))
      allocate (uvec(nploc,3))
      allocate (usizvec(nploc))
      allocate (uvvec(nploc,3))
      allocate (uvsizvec(nploc))
      allocate (vvec(nploc,3))
      allocate (vsizvec(nploc))
      allocate (wvec(nploc,3))
      allocate (wsizvec(nploc))
      allocate (devec1(3,nploc))
      allocate (axetypvec(nploc))
      allocate (uwvec(nploc,3),uwsizvec(nploc))
      allocate (vwvec(nploc,3),vwsizvec(nploc))
      allocate (usvec(nploc,3),ussizvec(nploc))
      allocate (vsvec(nploc,3),vssizvec(nploc))
      allocate (wsvec(nploc,3),wssizvec(nploc))
      allocate (urvec(nploc,3),ursizvec(nploc))
      allocate (t1vec(nploc,3),t1sizvec(nploc))
      allocate (t2vec(nploc,3),t2sizvec(nploc))
      allocate (rucosvec(nploc),rusinvec(nploc))
      allocate (rvcosvec(nploc),rvsinvec(nploc))
      allocate (rwcosvec(nploc),rwsinvec(nploc))
      allocate (uvcosvec(nploc),uvsinvec(nploc))
      allocate (uwcosvec(nploc),uwsinvec(nploc))
      allocate (vwcosvec(nploc),vwsinvec(nploc))
      allocate (upcosvec(nploc),upsinvec(nploc))
      allocate (urcosvec(nploc),ursinvec(nploc))
      allocate (uscosvec(nploc),ussinvec(nploc))
      allocate (vpcosvec(nploc),vpsinvec(nploc))
      allocate (vscosvec(nploc),vssinvec(nploc))
      allocate (wpcosvec(nploc),wpsinvec(nploc))
      allocate (wscosvec(nploc),wssinvec(nploc))
      allocate (ut1cosvec(nploc),ut1sinvec(nploc))
      allocate (ut2cosvec(nploc),ut2sinvec(nploc))
      allocate (dphiduvec(nploc),dphidvvec(nploc))
      allocate (dphidwvec(nploc),dphidrvec(nploc))
      allocate (dphiddelvec(nploc))
      allocate (dphidsvec(nploc))
      allocate (duvec(nploc,3),dvvec(nploc,3),dwvec(nploc,3))
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'torquevec'
c
c     zero out force components on local frame-defining atoms
c
      frcxvec (1:nploc,:)  = 0.0_ti_p
      frcyvec (1:nploc,:)  = 0.0_ti_p
      frczvec (1:nploc,:)  = 0.0_ti_p
      frcxvec1(1:nploc,:)  = 0.0_ti_p
      frcyvec1(1:nploc,:)  = 0.0_ti_p
      frczvec1(1:nploc,:)  = 0.0_ti_p
c
c     get the local frame type and the frame-defining atoms
c

      axetypvec(1:nploc)  = polaxe(ivec(1:nploc))

      mask(1:nploc) = (axetypvec(1:nploc) /= 'None') !mask for axetyp
                                                     !keep if axetyp is not None
      nploc1        = count(mask(1:nploc))

      if (nploc1.eq.0) return ! All axetyp are 'None', so return
cdeb     if(rank.eq.0.and.tinkerdebug) then
cdeb        write(*,*) 'TORQUEVEC IN nploc, nploc1, nploop,trqvec'
cdeb        write(*,*) nploc,nploc1,nploop
cdeb        write(*,*)
cdeb &      'trq1 ',trqvec(1:nploc,1),
cdeb &      'trq2 ',trqvec(1:nploc,2),
cdeb &      'trq3 ',trqvec(1:nploc,3)
cdeb        write(*,*) 'dep',
cdeb &      'devec1 ',devec(1,1:nploc),
cdeb &      'devec2 ',devec(2,1:nploc),
cdeb &      'devec3 ',devec(3,1:nploc)
cdeb        write(*,*) 'ivec' ,ivec(1:nploc)
cdeb        write(*,*) 'axetypvec',axetypvec(1:nploc)
cdeb     endif

      ivec1 (1:nploc1)      = pack (ivec (1:nploc)  , mask(1:nploc))
      iavec(1:nploc1) = zaxis(ivec1(1:nploc1))
      where (iavec(1:nploc1) > 0)
            ialocvec = loc(iavec)
      endwhere
      ibvec(1:nploc1) = ipole(ivec1(1:nploc1))
      iblocvec(1:nploc1) = loc(ibvec(1:nploc1))
      icvec(1:nploc1) = xaxis(ivec1(1:nploc1))
      where (icvec(1:nploc1) > 0)
            iclocvec = loc(icvec)
      endwhere
      idvec(1:nploc1) = yaxis(ivec1(1:nploc1))
      where (idvec(1:nploc1) > 0)
            idlocvec = loc(idvec)
      endwhere
c
c     construct the three rotation axes for the local frame
c
      uvec(1:nploc1,1) = x(iavec(1:nploc1)) - x(ibvec(1:nploc1))
      uvec(1:nploc1,2) = y(iavec(1:nploc1)) - y(ibvec(1:nploc1))
      uvec(1:nploc1,3) = z(iavec(1:nploc1)) - z(ibvec(1:nploc1))
      usizvec(1:nploc1) = sqrt(  uvec(1:nploc1,1)**2
     &                         + uvec(1:nploc1,2)**2
     &                         + uvec(1:nploc1,3)**2
     &                        )
      where (axetypvec(1:nploc1) /= 'Z-Only')
           vvec(:,1) = x(icvec(:)) - x(ibvec(:))
           vvec(:,2) = y(icvec(:)) - y(ibvec(:))
           vvec(:,3) = z(icvec(:)) - z(ibvec(:))
           vsizvec   = sqrt(vvec(:,1)**2 + vvec(:,2)**2 + vvec(:,3)**2)
      elsewhere
           vvec(:,1) = merge (0.0_ti_p, 1.0_ti_p,
     &           (abs(uvec(:,1) / usizvec)) > 0.866_ti_p)
           vvec(:,2) = merge (1.0_ti_p, 0.0_ti_p,
     &           (abs(uvec(:,1) / usizvec)) > 0.866_ti_p)
           vvec(:,3) = 0.0_ti_p
           vsizvec   = 1.0_ti_p
      endwhere
c
      where (axetypvec(1:nploc1) == 'Z-Bisect')
         wvec(:,1) = x(idvec) - x(ibvec)
         wvec(:,2) = y(idvec) - y(ibvec)
         wvec(:,3) = z(idvec) - z(ibvec)
      elsewhere
         wvec(:,1) =  uvec(:,2) * vvec(:,3) - uvec(:,3) * vvec(:,2)
         wvec(:,2) =  uvec(:,3) * vvec(:,1) - uvec(:,1) * vvec(:,3)
         wvec(:,3) =  uvec(:,1) * vvec(:,2) - uvec(:,2) * vvec(:,1)
      endwhere
c
      where (axetypvec(1:nploc1) == '3-Fold')
         wvec(:,1) = x(idvec) - x(ibvec)
         wvec(:,2) = y(idvec) - y(ibvec)
         wvec(:,3) = z(idvec) - z(ibvec)
      elsewhere
         wvec(:,1) =  uvec(:,2) * vvec(:,3) - uvec(:,3) * vvec(:,2)
         wvec(:,2) =  uvec(:,3) * vvec(:,1) - uvec(:,1) * vvec(:,3)
         wvec(:,3) =  uvec(:,1) * vvec(:,2) - uvec(:,2) * vvec(:,1)
      endwhere

      wsizvec(1:nploc1) = sqrt(  wvec(1:nploc1,1)**2
     &                         + wvec(1:nploc1,2)**2
     &                         + wvec(1:nploc1,3)**2
     &                        )
      do j = 1, 3
         uvec(1:nploc1,j) = uvec(1:nploc1,j) / usizvec(1:nploc1)
         vvec(1:nploc1,j) = vvec(1:nploc1,j) / vsizvec(1:nploc1)
         wvec(1:nploc1,j) = wvec(1:nploc1,j) / wsizvec(1:nploc1)
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      where (axetypvec(1:nploc1) == 'Z-Bisect')
         rvec(:,1)  = vvec(:,1) + wvec(:,1)
         rvec(:,2)  = vvec(:,2) + wvec(:,2)
         rvec(:,3)  = vvec(:,3) + wvec(:,3)
         rsizvec = sqrt( rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2 )

         svec(:,1) =  uvec(:,2) * rvec(:,3) - uvec(:,3) * rvec(:,2)
         svec(:,2) =  uvec(:,3) * rvec(:,1) - uvec(:,1) * rvec(:,3)
         svec(:,3) =  uvec(:,1) * rvec(:,2) - uvec(:,2) * rvec(:,1)
         ssizvec = sqrt( svec(:,1)**2 + svec(:,2)**2 + svec(:,3)**2 )
         rvec(:,1) = rvec(:,1) / rsizvec
         svec(:,1) = svec(:,1) / ssizvec
         rvec(:,2) = rvec(:,2) / rsizvec
         svec(:,2) = svec(:,2) / ssizvec
         rvec(:,3) = rvec(:,3) / rsizvec
         svec(:,3) = svec(:,3) / ssizvec
      endwhere
c
c     find the perpendicular and angle for each pair of axes
c
         uvvec(1:nploc1,1)  =  vvec(1:nploc1,2) * uvec(1:nploc1,3)
     &                       - vvec(1:nploc1,3) * uvec(1:nploc1,2)
         uvvec(1:nploc1,2)  =  vvec(1:nploc1,3) * uvec(1:nploc1,1)
     &                       - vvec(1:nploc1,1) * uvec(1:nploc1,3)
         uvvec(1:nploc1,3)  =  vvec(1:nploc1,1) * uvec(1:nploc1,2)
     &                       - vvec(1:nploc1,2) * uvec(1:nploc1,1)

         uwvec(1:nploc1,1)  =  wvec(1:nploc1,2) * uvec(1:nploc1,3)
     &                       - wvec(1:nploc1,3) * uvec(1:nploc1,2)
         uwvec(1:nploc1,2)  =  wvec(1:nploc1,3) * uvec(1:nploc1,1)
     &                       - wvec(1:nploc1,1) * uvec(1:nploc1,3)
         uwvec(1:nploc1,3)  =  wvec(1:nploc1,1) * uvec(1:nploc1,2)
     &                       - wvec(1:nploc1,2) * uvec(1:nploc1,1)

         vwvec(1:nploc1,1)  =  wvec(1:nploc1,2) * vvec(1:nploc1,3)
     &                       - wvec(1:nploc1,3) * vvec(1:nploc1,2)
         vwvec(1:nploc1,2)  =  wvec(1:nploc1,3) * vvec(1:nploc1,1)
     &                       - wvec(1:nploc1,1) * vvec(1:nploc1,3)
         vwvec(1:nploc1,3)  =  wvec(1:nploc1,1) * vvec(1:nploc1,2)
     &                       - wvec(1:nploc1,2) * vvec(1:nploc1,1)

         uvsizvec(1:nploc1) = sqrt(  uvvec(1:nploc1,1)**2
     &                             + uvvec(1:nploc1,2)**2
     &                             + uvvec(1:nploc1,3)**2
     &                            )
         uwsizvec(1:nploc1) = sqrt(  uwvec(1:nploc1,1)**2
     &                             + uwvec(1:nploc1,2)**2
     &                             + uwvec(1:nploc1,3)**2
     &                            )
         vwsizvec(1:nploc1) = sqrt(  vwvec(1:nploc1,1)**2
     &                             + vwvec(1:nploc1,2)**2
     &                             + vwvec(1:nploc1,3)**2
     &                            )
      do j = 1, 3
         uvvec(1:nploc1,j) = uvvec(1:nploc1,j) / uvsizvec(1:nploc1)
         uwvec(1:nploc1,j) = uwvec(1:nploc1,j) / uwsizvec(1:nploc1)
         vwvec(1:nploc1,j) = vwvec(1:nploc1,j) / vwsizvec(1:nploc1)
      end do
      where (axetypvec(1:nploc1) == 'Z-Bisect') 
         urvec(:,1)  =  rvec(:,2) * uvec(:,3) - rvec(:,3) * uvec(:,2)
         urvec(:,2)  =  rvec(:,3) * uvec(:,1) - rvec(:,1) * uvec(:,3)
         urvec(:,3)  =  rvec(:,1) * uvec(:,2) - rvec(:,2) * uvec(:,1)

         usvec(:,1)  =  svec(:,2) * uvec(:,3) - svec(:,3) * uvec(:,2)
         usvec(:,2)  =  svec(:,3) * uvec(:,1) - svec(:,1) * uvec(:,3)
         usvec(:,3)  =  svec(:,1) * uvec(:,2) - svec(:,2) * uvec(:,1)

         vsvec(:,1)  =  svec(:,2) * vvec(:,3) - svec(:,3) * vvec(:,2)
         vsvec(:,2)  =  svec(:,3) * vvec(:,1) - svec(:,1) * vvec(:,3)
         vsvec(:,3)  =  svec(:,1) * vvec(:,2) - svec(:,2) * vvec(:,1)

         wsvec(:,1)  =  svec(:,2) * wvec(:,3) - svec(:,3) * wvec(:,2)
         wsvec(:,2)  =  svec(:,3) * wvec(:,1) - svec(:,1) * wvec(:,3)
         wsvec(:,3)  =  svec(:,1) * wvec(:,2) - svec(:,2) * wvec(:,1)

         ursizvec = sqrt(urvec(:,1)**2 + urvec(:,2)**2 + urvec(:,3)**2)
         ussizvec = sqrt(usvec(:,1)**2 + usvec(:,2)**2 + usvec(:,3)**2)
         vssizvec = sqrt(vsvec(:,1)**2 + vsvec(:,2)**2 + vsvec(:,3)**2)
         wssizvec = sqrt(wsvec(:,1)**2 + wsvec(:,2)**2 + wsvec(:,3)**2)
         urvec(:,1)  = urvec(:,1) / ursizvec
         usvec(:,1)  = usvec(:,1) / ussizvec
         vsvec(:,1)  = vsvec(:,1) / vssizvec
         wsvec(:,1)  = wsvec(:,1) / wssizvec
         urvec(:,2)  = urvec(:,2) / ursizvec
         usvec(:,2)  = usvec(:,2) / ussizvec
         vsvec(:,2)  = vsvec(:,2) / vssizvec
         wsvec(:,2)  = wsvec(:,2) / wssizvec
         urvec(:,3)  = urvec(:,3) / ursizvec
         usvec(:,3)  = usvec(:,3) / ussizvec
         vsvec(:,3)  = vsvec(:,3) / vssizvec
         wsvec(:,3)  = wsvec(:,3) / wssizvec
      endwhere
c
c     get sine and cosine of angles between the rotation axes
c
      uvcosvec(1:nploc1) =  uvec(1:nploc1,1) * vvec(1:nploc1,1)
     &                    + uvec(1:nploc1,2) * vvec(1:nploc1,2)
     &                    + uvec(1:nploc1,3) * vvec(1:nploc1,3)
      uvsinvec(1:nploc1) =  sqrt(1.0_ti_p - uvcosvec(1:nploc1)**2)
      uwcosvec(1:nploc1) =  uvec(1:nploc1,1) * wvec(1:nploc1,1)
     &                    + uvec(1:nploc1,2) * wvec(1:nploc1,2)
     &                    + uvec(1:nploc1,3) * wvec(1:nploc1,3)
      uwsinvec(1:nploc1) =  sqrt(1.0_ti_p - uwcosvec(1:nploc1)**2)
      vwcosvec(1:nploc1) =  vvec(1:nploc1,1) * wvec(1:nploc1,1)
     &                    + vvec(1:nploc1,2) * wvec(1:nploc1,2)
     &                    + vvec(1:nploc1,3) * wvec(1:nploc1,3)
      vwsinvec(1:nploc1) =  sqrt(1.0_ti_p - vwcosvec(1:nploc1)**2)
c
      where (axetypvec(1:nploc1) == 'Z-Bisect') 
         urcosvec =  uvec(:,1) * rvec(:,1) + uvec(:,2) * rvec(:,2)
     &             + uvec(:,3) * rvec(:,3)
         ursinvec =  sqrt(1.0_ti_p - urcosvec(1:nploc1)**2)
         uscosvec =  uvec(:,1) * svec(:,1) + uvec(:,2) * svec(:,2)
     &             + uvec(:,3) * svec(:,3)
         ussinvec =  sqrt(1.0_ti_p - uscosvec(1:nploc1)**2)
         vscosvec =  vvec(:,1) * svec(:,1) + vvec(:,2) * svec(:,2)
     &             + vvec(:,3) * svec(:,3)
         vssinvec =  sqrt(1.0_ti_p - vscosvec(1:nploc1)**2)
         wscosvec =  wvec(:,1) * svec(:,1) + wvec(:,2) * svec(:,2)
     &             + wvec(:,3) * svec(:,3)
         wssinvec =  sqrt(1.0_ti_p - wscosvec(1:nploc1)**2)

c
c     compute the projection of v and w onto the ru-plane
c
         t1vec(:,1)  =  vvec(:,1) - svec(:,1) * vscosvec
         t2vec(:,1)  =  wvec(:,1) - svec(:,1) * wscosvec
         t1vec(:,2)  =  vvec(:,2) - svec(:,2) * vscosvec
         t2vec(:,2)  =  wvec(:,2) - svec(:,2) * wscosvec
         t1vec(:,3)  =  vvec(:,3) - svec(:,3) * vscosvec
         t2vec(:,3)  =  wvec(:,3) - svec(:,3) * wscosvec

         t1sizvec = sqrt(t1vec(:,1)**2 + t1vec(:,2)**2 + t1vec(:,3)**2)
         t2sizvec = sqrt(t2vec(:,1)**2 + t2vec(:,2)**2 + t2vec(:,3)**2)
         t1vec(:,1)  = t1vec(:,1) / t1sizvec
         t1vec(:,2)  = t1vec(:,2) / t1sizvec
         t1vec(:,3)  = t1vec(:,3) / t1sizvec
         t2vec(:,1)  = t2vec(:,1) / t2sizvec
         t2vec(:,2)  = t2vec(:,2) / t2sizvec
         t2vec(:,3)  = t2vec(:,3) / t2sizvec

         ut1cosvec =  uvec(:,1) * t1vec(:,1) + uvec(:,2) * t1vec(:,2)
     &              + uvec(:,3) * t1vec(:,3)
         ut1sinvec =  sqrt(1.0_ti_p - ut1cosvec(1:nploc1)**2)
         ut2cosvec =  uvec(:,1) * t2vec(:,1) + uvec(:,2) * t2vec(:,2)
     &              + uvec(:,3) * t2vec(:,3)
         ut2sinvec =  sqrt(1.0_ti_p - ut2cosvec(1:nploc1)**2)
      endwhere
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphiduvec(1:nploc1) = - trqvec(1:nploc1,1) * uvec(1:nploc1,1)
     &                      - trqvec(1:nploc1,2) * uvec(1:nploc1,2)
     &                      - trqvec(1:nploc1,3) * uvec(1:nploc1,3)
      dphidvvec(1:nploc1) = - trqvec(1:nploc1,1) * vvec(1:nploc1,1)
     &                      - trqvec(1:nploc1,2) * vvec(1:nploc1,2)
     &                      - trqvec(1:nploc1,3) * vvec(1:nploc1,3)
      dphidwvec(1:nploc1) = - trqvec(1:nploc1,1) * wvec(1:nploc1,1)
     &                      - trqvec(1:nploc1,2) * wvec(1:nploc1,2)
     &                      - trqvec(1:nploc1,3) * wvec(1:nploc1,3)

      where (axetypvec(1:nploc1) == 'Z-Bisect')
        dphidrvec = - trqvec(:,1) * rvec(:,1)
     &              - trqvec(:,2) * rvec(:,2)
     &              - trqvec(:,3) * rvec(:,3)
        dphidsvec = - trqvec(:,1) * svec(:,1)
     &              - trqvec(:,2) * svec(:,2)
     &              - trqvec(:,3) * svec(:,3)
      endwhere
c
c     force distribution for the Z-Only local coordinate method
c
      where (axetypvec(1:nploc1) == 'Z-Only')
         duvec(:,1) =  uvvec (:,1)  * dphidvvec / (usizvec * uvsinvec )
     &               + uwvec (:,1)  * dphidwvec /  usizvec
         duvec(:,2) =  uvvec (:,2)  * dphidvvec / (usizvec * uvsinvec )
     &               + uwvec (:,2)  * dphidwvec /  usizvec
         duvec(:,3) =  uvvec (:,3)  * dphidvvec / (usizvec * uvsinvec )
     &               + uwvec (:,3)  * dphidwvec /  usizvec

         devec(1,ialocvec) =  devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) =  devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) =  devec(3,ialocvec) + duvec(:,3)

         devec(1,iblocvec) =  devec(1,iblocvec) - duvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - duvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - duvec(:,3)

         frczvec1(:,1) =  frczvec1(:,1) + duvec(:,1)
         frczvec1(:,2) =  frczvec1(:,2) + duvec(:,2)
         frczvec1(:,3) =  frczvec1(:,3) + duvec(:,3)
          
c
c     force distribution for the Z-then-X local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == 'Z-then-X')
         duvec(:,1) =   uvvec(:,1) * dphidvvec / (usizvec * uvsinvec)
     &                + uwvec(:,1) * dphidwvec /  usizvec
         duvec(:,2) =   uvvec(:,2) * dphidvvec / (usizvec * uvsinvec)
     &                + uwvec(:,2) * dphidwvec /  usizvec
         duvec(:,3) =   uvvec(:,3) * dphidvvec / (usizvec * uvsinvec)
     &                + uwvec(:,3) * dphidwvec /  usizvec

         dvvec(:,1) = - uvvec(:,1) * dphiduvec / (vsizvec* uvsinvec)
         dvvec(:,2) = - uvvec(:,2) * dphiduvec / (vsizvec* uvsinvec)
         dvvec(:,3) = - uvvec(:,3) * dphiduvec / (vsizvec* uvsinvec)

         devec(1,ialocvec) = devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) = devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) = devec(3,ialocvec) + duvec(:,3)

         devec(1,iclocvec) = devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) = devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) = devec(3,iclocvec) + dvvec(:,3)

         devec(1,iblocvec) = devec(1,iblocvec) - duvec(:,1) - dvvec(:,1)
         devec(2,iblocvec) = devec(2,iblocvec) - duvec(:,2) - dvvec(:,2)
         devec(3,iblocvec) = devec(3,iblocvec) - duvec(:,3) - dvvec(:,3)

         frczvec1(:,1)     =  frczvec1(:,1) + duvec(:,1)
         frczvec1(:,2)     =  frczvec1(:,2) + duvec(:,2)
         frczvec1(:,3)     =  frczvec1(:,3) + duvec(:,3)
         frcxvec1(:,1)     =  frcxvec1(:,1) + dvvec(:,1)
         frcxvec1(:,2)     =  frcxvec1(:,2) + dvvec(:,2)
         frcxvec1(:,3)     =  frcxvec1(:,3) + dvvec(:,3)


c
c     force distribution for the Bisector local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == 'Bisector')
         duvec(:,1) =      uvvec(:,1) * dphidvvec / (usizvec * uvsinvec)
     &          +  0.5_ti_p * uwvec(:,1) * dphidwvec /  usizvec
         duvec(:,2) =      uvvec(:,2) * dphidvvec / (usizvec * uvsinvec)
     &          +  0.5_ti_p * uwvec(:,2) * dphidwvec /  usizvec
         duvec(:,3) =      uvvec(:,3) * dphidvvec / (usizvec * uvsinvec)
     &          +  0.5_ti_p * uwvec(:,3) * dphidwvec /  usizvec
         dvvec(:,1) =   -  uvvec(:,1) * dphiduvec / (vsizvec * uvsinvec)
     &          +  0.5_ti_p * vwvec(:,1) * dphidwvec /  vsizvec
         dvvec(:,2) =   -  uvvec(:,2) * dphiduvec / (vsizvec * uvsinvec)
     &          +  0.5_ti_p * vwvec(:,2) * dphidwvec /  vsizvec
         dvvec(:,3) =   -  uvvec(:,3) * dphiduvec / (vsizvec * uvsinvec)
     &          +  0.5_ti_p * vwvec(:,3) * dphidwvec /  vsizvec

         devec(1,ialocvec) = devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) = devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) = devec(3,ialocvec) + duvec(:,3)

         devec(1,iclocvec) = devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) = devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) = devec(3,iclocvec) + dvvec(:,3)

         devec(1,iblocvec) = devec(1,iblocvec) - duvec(:,1) - dvvec(:,1)
         devec(2,iblocvec) = devec(2,iblocvec) - duvec(:,2) - dvvec(:,2)
         devec(3,iblocvec) = devec(3,iblocvec) - duvec(:,3) - dvvec(:,3)

         frczvec1(:,1)     = frczvec1(:,1)     + duvec(:,1)
         frczvec1(:,2)     = frczvec1(:,2)     + duvec(:,2)
         frczvec1(:,3)     = frczvec1(:,3)     + duvec(:,3)
         frcxvec1(:,1)     = frcxvec1(:,1)     + dvvec(:,1)
         frcxvec1(:,2)     = frcxvec1(:,2)     + dvvec(:,2)
         frcxvec1(:,3)     = frcxvec1(:,3)     + dvvec(:,3)



c
c     force distribution for the Z-Bisect local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == 'Z-Bisect')
         duvec(:,1) =  urvec(:,1) * dphidrvec / (usizvec * ursinvec)
     &               + usvec(:,1) * dphidsvec /  usizvec
         duvec(:,2) =  urvec(:,2) * dphidrvec / (usizvec * uvsinvec)
     &               + usvec(:,2) * dphidsvec /  usizvec
         duvec(:,3) =  urvec(:,3) * dphidrvec / (usizvec * uvsinvec)
     &               + usvec(:,3) * dphidsvec /  usizvec
         dvvec(:,1) =  (vssinvec * svec (:,1) - vscosvec   * t1vec(:,1))
     &               * dphiduvec / (vsizvec   * (ut1sinvec + ut2sinvec))
         dvvec(:,2) =  (vssinvec * svec (:,2) - vscosvec   * t1vec(:,2))
     &               * dphiduvec / (vsizvec   * (ut1sinvec + ut2sinvec))
         dvvec(:,3) =  (vssinvec * svec (:,3) - vscosvec   * t1vec(:,3))
     &               * dphiduvec / (vsizvec   * (ut1sinvec + ut2sinvec))
         dwvec(:,1) =  (wssinvec * svec (:,1) - wscosvec   * t2vec(:,1))
     &               * dphiduvec / (wsizvec   * (ut1sinvec + ut2sinvec))
         dwvec(:,2) =  (wssinvec * svec (:,2) - wscosvec   * t2vec(:,2))
     &               * dphiduvec / (wsizvec   * (ut1sinvec + ut2sinvec))
         dwvec(:,3) =  (wssinvec * svec (:,3) - wscosvec   * t2vec(:,3))
     &               * dphiduvec / (wsizvec   * (ut1sinvec + ut2sinvec))

         devec(1,ialocvec) = devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) = devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) = devec(3,ialocvec) + duvec(:,3)

         devec(1,iclocvec) = devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) = devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) = devec(3,iclocvec) + dvvec(:,3)

         devec(1,idlocvec) = devec(1,idlocvec) + dwvec(:,1)
         devec(2,idlocvec) = devec(2,idlocvec) + dwvec(:,2)
         devec(3,idlocvec) = devec(3,idlocvec) + dwvec(:,3)

         devec(1,iblocvec) = devec(1,iblocvec) - duvec(:,1) - dvvec(:,1)
     &                                         - dwvec(:,1)
         devec(2,iblocvec) = devec(2,iblocvec) - duvec(:,2) - dvvec(:,2)
     &                                         - dwvec(:,2)
         devec(3,iblocvec) = devec(3,iblocvec) - duvec(:,3) - dvvec(:,3)
     &                                         - dwvec(:,3)
         frczvec1(:,1)     =  frczvec1(:,1)    + duvec(:,1)
         frczvec1(:,2)     =  frczvec1(:,2)    + duvec(:,2)
         frczvec1(:,3)     =  frczvec1(:,3)    + duvec(:,3)

         frcxvec1(:,1)     =  frcxvec1(:,1)    + dvvec(:,1)
         frcxvec1(:,2)     =  frcxvec1(:,2)    + dvvec(:,2)
         frcxvec1(:,3)     =  frcxvec1(:,3)    + dvvec(:,3)

         frcyvec1(:,1)     =  frcyvec1(:,1)    + dwvec(:,1)
         frcyvec1(:,2)     =  frcyvec1(:,2)    + dwvec(:,2)
         frcyvec1(:,3)     =  frcyvec1(:,3)    + dwvec(:,3)
c
c     force distribution for the 3-Fold local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == '3-Fold')
         pvec(:,1) = uvec(:,1) + vvec(:,1) + wvec(:,1)
         pvec(:,2) = uvec(:,2) + vvec(:,2) + wvec(:,2)
         pvec(:,3) = uvec(:,3) + vvec(:,3) + wvec(:,3)
         psizvec   = sqrt(pvec(:,1)**2 + pvec(:,2)**2 + pvec(:,3)**2)

         pvec(:,1) = pvec(:,1) / psizvec
         pvec(:,2) = pvec(:,2) / psizvec
         pvec(:,3) = pvec(:,3) / psizvec
         upcosvec  =  uvec(:,1) * pvec(:,1) + uvec(:,2) * pvec(:,2)
     &              + uvec(:,3) * pvec(:,3)
         vpcosvec  =  vvec(:,1) * pvec(:,1) + vvec(:,2) * pvec(:,2)
     &              + vvec(:,3) * pvec(:,3)
         wpcosvec  =  wvec(:,1) * pvec(:,1) + wvec(:,2) * pvec(:,2)
     &              + wvec(:,3) * pvec(:,3)
         upsinvec  =  sqrt(1.0_ti_p - upcosvec**2)
         vpsinvec  =  sqrt(1.0_ti_p - vpcosvec**2)
         wpsinvec  =  sqrt(1.0_ti_p - wpcosvec**2)

         rvec(:,1) = uvec(:,1) + vvec(:,1)
         rvec(:,2) = uvec(:,2) + vvec(:,2)
         rvec(:,3) = uvec(:,3) + vvec(:,3)
         rsizvec   = sqrt(rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2)

         rvec(:,1) = rvec(:,1) / rsizvec
         rvec(:,2) = rvec(:,2) / rsizvec
         rvec(:,3) = rvec(:,3) / rsizvec
         rwcosvec  =  rvec(:,1) * wvec(:,1) + rvec(:,2) * wvec(:,2)
     &              + rvec(:,3) * wvec(:,3)
         rwsinvec  = sqrt(1.0_ti_p - rwcosvec**2)
         dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &               - trqvec(:,3) * rvec(:,3)

         delvec(:,1)  =  rvec(:,2) * wvec(:,3) - rvec(:,3) * wvec(:,2)
         delvec(:,2)  =  rvec(:,3) * wvec(:,1) - rvec(:,1) * wvec(:,3)
         delvec(:,3)  =  rvec(:,1) * wvec(:,2) - rvec(:,2) * wvec(:,1)
         delsizvec    =  sqrt(  delvec(:,1)**2 + delvec(:,2)**2
     &                        + delvec(:,3)**2 )
         delvec(:,1)  =  delvec(:,1) / delsizvec
         delvec(:,2)  =  delvec(:,2) / delsizvec
         delvec(:,3)  =  delvec(:,3) / delsizvec
         dphiddelvec  = - trqvec(:,1) * delvec(:,1)
     &                  - trqvec(:,2) * delvec(:,2)
     &                  - trqvec(:,3) * delvec(:,3)

         epsvec(:,1) = delvec(:,2) * wvec(:,3) - delvec(:,3) * wvec(:,2)
         epsvec(:,2) = delvec(:,3) * wvec(:,1) - delvec(:,1) * wvec(:,3)
         epsvec(:,3) = delvec(:,1) * wvec(:,2) - delvec(:,2) * wvec(:,1)
         dwvec(:,1)  = delvec(:,1) * dphidrvec   / (wsizvec * rwsinvec)
     &              +  epsvec(:,1) * dphiddelvec * wpcosvec
     &                                           / (wsizvec * psizvec )
         dwvec(:,2) =  delvec(:,2) * dphidrvec   / (wsizvec * rwsinvec)
     &              +  epsvec(:,2) * dphiddelvec * wpcosvec
     &                                           / (wsizvec * psizvec )
         dwvec(:,3) =  delvec(:,3) * dphidrvec   / (wsizvec * rwsinvec)
     &              +  epsvec(:,3) * dphiddelvec * wpcosvec
     &                                           / (wsizvec * psizvec )
         devec(1,idlocvec) =  devec(1,idlocvec) + dwvec(:,1)
         devec(2,idlocvec) =  devec(2,idlocvec) + dwvec(:,2)
         devec(3,idlocvec) =  devec(3,idlocvec) + dwvec(:,3)
         devec(1,iblocvec) =  devec(1,iblocvec) - dwvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - dwvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - dwvec(:,3)
         frcyvec1(:,1)     =  frcyvec1  (:,1)   + dwvec(:,1)
         frcyvec1(:,2)     =  frcyvec1  (:,2)   + dwvec(:,2)
         frcyvec1(:,3)     =  frcyvec1  (:,3)   + dwvec(:,3)

         rvec(:,1)  =  vvec(:,1) + wvec(:,1)
         rvec(:,2)  =  vvec(:,2) + wvec(:,2)
         rvec(:,3)  =  vvec(:,3) + wvec(:,3)
         rsizvec = sqrt(rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2)

         rvec(:,1) =  rvec(:,1) / rsizvec
         rvec(:,2) =  rvec(:,2) / rsizvec
         rvec(:,3) =  rvec(:,3) / rsizvec
         rucosvec  =  rvec(:,1) * uvec(:,1) + rvec(:,2) * uvec(:,2)
     &              + rvec(:,3) * uvec(:,3)
         rusinvec  = sqrt(1.0_ti_p - rucosvec**2)

         dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &               - trqvec(:,3) * rvec(:,3)

         delvec(:,1) =  rvec(:,2) * uvec(:,3) - rvec(:,3) * uvec(:,2)
         delvec(:,2) =  rvec(:,3) * uvec(:,1) - rvec(:,1) * uvec(:,3)
         delvec(:,3) =  rvec(:,1) * uvec(:,2) - rvec(:,2) * uvec(:,1)

         delsizvec = sqrt(  delvec(:,1)**2 + delvec(:,2)**2
     &                    + delvec(:,3)**2 )
         delvec(:,1) =  delvec(:,1) / delsizvec
         delvec(:,2) =  delvec(:,2) / delsizvec
         delvec(:,3) =  delvec(:,3) / delsizvec
         dphiddelvec = - trqvec(:,1) * delvec(:,1)
     &                 - trqvec(:,2) * delvec(:,2)
     &                 - trqvec(:,3) * delvec(:,3)
         epsvec(:,1) = delvec(:,2) * uvec(:,3) - delvec(:,3) * uvec(:,2)
         epsvec(:,2) = delvec(:,3) * uvec(:,1) - delvec(:,1) * uvec(:,3)
         epsvec(:,3) = delvec(:,1) * uvec(:,2) - delvec(:,2) * uvec(:,1)

         duvec(:,1) =  delvec(:,1) * dphidrvec   / (usizvec * rusinvec)
     &               + epsvec(:,1) * dphiddelvec * upcosvec
     &                                           / (usizvec * psizvec )
         duvec(:,2) =  delvec(:,2) * dphidrvec   / (usizvec * rusinvec)
     &               + epsvec(:,2) * dphiddelvec * upcosvec
     &                                           / (usizvec * psizvec )
         duvec(:,3) =  delvec(:,3) * dphidrvec   / (usizvec * rusinvec)
     &               + epsvec(:,3) * dphiddelvec * upcosvec
     &                                           / (usizvec * psizvec )
         devec(1,ialocvec) =  devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) =  devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) =  devec(3,ialocvec) + duvec(:,3)

         devec(1,iblocvec) =  devec(1,iblocvec) - duvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - duvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - duvec(:,3)

         frczvec1(:,1)     =  frczvec1(:,1)     + duvec(:,1)
         frczvec1(:,2)     =  frczvec1(:,2)     + duvec(:,2)
         frczvec1(:,3)     =  frczvec1(:,3)     + duvec(:,3)

         rvec(:,1)  =  uvec(:,1) + wvec(:,1)
         rvec(:,2)  =  uvec(:,2) + wvec(:,2)
         rvec(:,3)  =  uvec(:,3) + wvec(:,3)
         rsizvec = sqrt( rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2 )

         rvec(:,1) =  rvec(:,1) / rsizvec
         rvec(:,2) =  rvec(:,2) / rsizvec
         rvec(:,3) =  rvec(:,3) / rsizvec
         rvcosvec  =  rvec(:,1) * vvec(:,1) + rvec(:,2) * vvec(:,2)
     &              + rvec(:,3) * vvec(:,3)
         rvsinvec = sqrt(1.0_ti_p - rvcosvec**2)
         dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &               - trqvec(:,3) * rvec(:,3)

         delvec(:,1) =  rvec(:,2) * vvec(:,3) - rvec(:,3) * vvec(:,2)
         delvec(:,2) =  rvec(:,3) * vvec(:,1) - rvec(:,1) * vvec(:,3)
         delvec(:,3) =  rvec(:,1) * vvec(:,2) - rvec(:,2) * vvec(:,1)

         delsizvec = sqrt(  delvec(:,1)**2 + delvec(:,2)**2
     &                    + delvec(:,3)**2 )
        delvec(:,1) =  delvec(:,1) / delsizvec
        delvec(:,2) =  delvec(:,2) / delsizvec
        delvec(:,3) =  delvec(:,3) / delsizvec

        dphiddelvec = - trqvec(:,1) * delvec(:,1)
     &                - trqvec(:,2) * delvec(:,2)
     &                - trqvec(:,3) * delvec(:,3)
         epsvec(:,1) = delvec(:,2) * vvec(:,3) - delvec(:,3) * vvec(:,2)
         epsvec(:,2) = delvec(:,3) * vvec(:,1) - delvec(:,1) * vvec(:,3)
         epsvec(:,3) = delvec(:,1) * vvec(:,2) - delvec(:,2) * vvec(:,1)

         dvvec(:,1)  =  delvec(:,1) * dphidrvec   / (vsizvec * rvsinvec)
     &                + epsvec(:,1) * dphiddelvec * vpcosvec
     &                                            / (vsizvec * psizvec)
         dvvec(:,2)  =  delvec(:,2) * dphidrvec   / (vsizvec * rvsinvec)
     &                + epsvec(:,2) * dphiddelvec * vpcosvec
     &                                            / (vsizvec * psizvec)
         dvvec(:,3)  =  delvec(:,3) * dphidrvec   / (vsizvec * rvsinvec)
     &                + epsvec(:,3) * dphiddelvec * vpcosvec
     &                                            / (vsizvec * psizvec)
         devec(1,iclocvec) =  devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) =  devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) =  devec(3,iclocvec) + dvvec(:,3)

         devec(1,iblocvec) =  devec(1,iblocvec) - dvvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - dvvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - dvvec(:,3)
         frcxvec1(:,1)     =  frcxvec1(:,1)     + dvvec(:,1)
         frcxvec1(:,2)     =  frcxvec1(:,2)     + dvvec(:,2)
         frcxvec1(:,3)     =  frcxvec1(:,3)     + dvvec(:,3)
      endwhere
         
cnew
cnew  send the results back in the original arrays
      do i=1,3
         frcxvec(1:nploc,i) =
     &                unpack(frcxvec1(1:nploc1,i), mask(1:nploc), 0.0_ti_p)
         frcyvec(1:nploc,i) =
     &                unpack(frcyvec1(1:nploc1,i), mask(1:nploc), 0.0_ti_p)
         frczvec(1:nploc,i) =
     &                unpack(frczvec1(1:nploc1,i), mask(1:nploc), 0.0_ti_p)
      enddo
cdeb  if(rank.eq.0.and.tinkerdebug) write(*,*) 'torquevec'
cdeb     if(rank.eq.0.and.tinkerdebug) then
cdeb        write(*,*) 'TORQUEVEC OUT nploc, trqvec'
cdeb        write(*,*) nploc,nploc1,nploop
cdeb        write(*,*)
cdeb &      'fixvec1 ',frcxvec(1:nploc,1),
cdeb &      'fixvec2 ',frcxvec(1:nploc,2),
cdeb &      'fixvec3 ',frcxvec(1:nploc,3),
cdeb &      'fiyvec1 ',frcyvec(1:nploc,1),
cdeb &      'fiyvec2 ',frcyvec(1:nploc,2),
cdeb &      'fiyvec3 ',frcyvec(1:nploc,3),
cdeb &      'fizvec1 ',frczvec(1:nploc,1),
cdeb &      'fizvec2 ',frczvec(1:nploc,2),
cdeb &      'fizvec3 ',frczvec(1:nploc,3)
cdeb     endif


      deallocate (usizvec)
      deallocate (uvec)
      deallocate (ssizvec)
      deallocate (svec)
      deallocate (rsizvec)
      deallocate (rvec)
      deallocate (idlocvec)
      deallocate (iclocvec)
      deallocate (iblocvec)
      deallocate (ialocvec)
      deallocate (idvec)
      deallocate (icvec)
      deallocate (ibvec)
      deallocate (iavec)
      deallocate (mask)
      deallocate (ivec1)
      deallocate (uvvec)
      deallocate (uvsizvec)
      deallocate (vvec)
      deallocate (vsizvec)
      deallocate (wvec)
      deallocate (wsizvec)
      deallocate (axetypvec)
      deallocate (uwvec,uwsizvec)
      deallocate (vwvec,vwsizvec)
      deallocate (usvec,ussizvec)
      deallocate (vsvec,vssizvec)
      deallocate (wsvec,wssizvec)
      deallocate (urvec,ursizvec)
      deallocate (t1vec,t1sizvec)
      deallocate (t2vec,t2sizvec)
      deallocate (rucosvec,rusinvec)
      deallocate (rvcosvec,rvsinvec)
      deallocate (rwcosvec,rwsinvec)
      deallocate (uvcosvec,uvsinvec)
      deallocate (uwcosvec,uwsinvec)
      deallocate (vwcosvec,vwsinvec)
      deallocate (upcosvec,upsinvec)
      deallocate (urcosvec,ursinvec)
      deallocate (uscosvec,ussinvec)
      deallocate (vpcosvec,vpsinvec)
      deallocate (vscosvec,vssinvec)
      deallocate (wpcosvec,wpsinvec)
      deallocate (wscosvec,wssinvec)
      deallocate (ut1cosvec,ut1sinvec)
      deallocate (ut2cosvec,ut2sinvec)
      deallocate (dphiduvec,dphidvvec)
      deallocate (dphidwvec,dphidrvec)
      deallocate (dphiddelvec)
      deallocate (dphidsvec)
      deallocate (duvec,dvvec,dwvec)

      return
      end
      subroutine torquevec_rec (nploc,ivec,trqvec,
     &                          frcxvec,frcyvec,frczvec,devec)
      use atoms
      use deriv
      use domdec
      use mpole
      use sizes

      implicit none
      integer i,j,nploc,nploc1
      integer ivec(*)
      integer ialoc,ibloc,icloc,idloc
      real(t_p) trqvec(nploc,3)
      real(t_p) frcxvec(nploc,3),frcyvec(nploc,3), frczvec(nploc,3)
      real(t_p) devec(3,*)
      !DIR$ ATTRIBUTES ALIGN:64::mask
      logical,allocatable :: mask(:)
      !DIR$ ATTRIBUTES ALIGN:64::ivec1
      integer,allocatable :: ivec1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: frcxvec1,frcyvec1,frczvec1
      real(t_p),allocatable :: frcxvec1(:,:),frcyvec1(:,:),frczvec1(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: devec1
      real(t_p),allocatable :: devec1(:,:)
      !DIR$ ATTRIBUTES ALIGN:64::iavec,ibvec,icvec,idvec
      integer,allocatable :: iavec(:),ibvec(:),icvec(:),idvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ialocvec,iblocvec
      integer,allocatable :: ialocvec(:),iblocvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::iclocvec,idlocvec
      integer,allocatable :: iclocvec(:),idlocvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::delvec,delsizvec
      real(t_p),allocatable :: delvec(:,:),delsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::epsvec,epssizvec
      real(t_p),allocatable :: epsvec(:,:),epssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::pvec,psizvec
      real(t_p),allocatable :: pvec(:,:),psizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rvec,rsizvec
      real(t_p),allocatable :: rvec(:,:),rsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::svec,ssizvec
      real(t_p),allocatable :: svec(:,:),ssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uvec,usizvec
      real(t_p),allocatable :: uvec(:,:),usizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vvec,vsizvec
      real(t_p),allocatable :: vvec(:,:),vsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wvec,wsizvec
      real(t_p),allocatable :: wvec(:,:),wsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uvvec,uvsizvec
      real(t_p),allocatable :: uvvec(:,:),uvsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uwvec,uwsizvec
      real(t_p),allocatable :: uwvec(:,:),uwsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vwvec,vwsizvec
      real(t_p),allocatable :: vwvec(:,:),vwsizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::usvec,ussizvec
      real(t_p),allocatable :: usvec(:,:),ussizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vsvec,vssizvec
      real(t_p),allocatable :: vsvec(:,:),vssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wsvec,wssizvec
      real(t_p),allocatable :: wsvec(:,:),wssizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::urvec,ursizvec
      real(t_p),allocatable :: urvec(:,:),ursizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::t1vec,t1sizvec
      real(t_p),allocatable :: t1vec(:,:),t1sizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::t2vec,t2sizvec
      real(t_p),allocatable :: t2vec(:,:),t2sizvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rucosvec,rusinvec
      real(t_p),allocatable :: rucosvec(:),rusinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rvcosvec,rvsinvec
      real(t_p),allocatable :: rvcosvec(:),rvsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::rwcosvec,rwsinvec
      real(t_p),allocatable :: rwcosvec(:),rwsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uvcosvec,uvsinvec
      real(t_p),allocatable :: uvcosvec(:),uvsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uwcosvec,uwsinvec
      real(t_p),allocatable :: uwcosvec(:),uwsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vwcosvec,vwsinvec
      real(t_p),allocatable :: vwcosvec(:),vwsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::upcosvec,upsinvec
      real(t_p),allocatable :: upcosvec(:),upsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::urcosvec,ursinvec
      real(t_p),allocatable :: urcosvec(:),ursinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::uscosvec,ussinvec
      real(t_p),allocatable :: uscosvec(:),ussinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vpcosvec,vpsinvec
      real(t_p),allocatable :: vpcosvec(:),vpsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::vscosvec,vssinvec
      real(t_p),allocatable :: vscosvec(:),vssinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wpcosvec,wpsinvec
      real(t_p),allocatable :: wpcosvec(:),wpsinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::wscosvec,wssinvec
      real(t_p),allocatable :: wscosvec(:),wssinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ut1cosvec,ut1sinvec
      real(t_p),allocatable :: ut1cosvec(:),ut1sinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ut2cosvec,ut2sinvec
      real(t_p),allocatable :: ut2cosvec(:),ut2sinvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphiduvec,dphidvvec
      real(t_p),allocatable :: dphiduvec(:),dphidvvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphidwvec,dphidrvec
      real(t_p),allocatable :: dphidwvec(:),dphidrvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphiddelvec
      real(t_p),allocatable :: dphiddelvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dphidsvec
      real(t_p),allocatable :: dphidsvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::duvec,dvvec,dwvec
      real(t_p),allocatable :: duvec(:,:),dvvec(:,:),dwvec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64::axetypvec
      character*8,allocatable :: axetypvec(:)

      allocate (ivec1(nploc))
      allocate (mask(nploc))
      allocate (iavec(nploc))
      allocate (ibvec(nploc))
      allocate (icvec(nploc))
      allocate (idvec(nploc))
      allocate (ialocvec(nploc))
      allocate (iblocvec(nploc))
      allocate (iclocvec(nploc))
      allocate (idlocvec(nploc))
      allocate (frcxvec1(nploc,3))
      allocate (frcyvec1(nploc,3))
      allocate (frczvec1(nploc,3))
      allocate (rvec(nploc,3))
      allocate (rsizvec(nploc))
      allocate (svec(nploc,3))
      allocate (ssizvec(nploc))
      allocate (uvec(nploc,3))
      allocate (usizvec(nploc))
      allocate (uvvec(nploc,3))
      allocate (uvsizvec(nploc))
      allocate (vvec(nploc,3))
      allocate (vsizvec(nploc))
      allocate (wvec(nploc,3))
      allocate (wsizvec(nploc))
      allocate (devec1(3,nploc))
      allocate (axetypvec(nploc))
      allocate (uwvec(nploc,3),uwsizvec(nploc))
      allocate (vwvec(nploc,3),vwsizvec(nploc))
      allocate (usvec(nploc,3),ussizvec(nploc))
      allocate (vsvec(nploc,3),vssizvec(nploc))
      allocate (wsvec(nploc,3),wssizvec(nploc))
      allocate (urvec(nploc,3),ursizvec(nploc))
      allocate (t1vec(nploc,3),t1sizvec(nploc))
      allocate (t2vec(nploc,3),t2sizvec(nploc))
      allocate (rucosvec(nploc),rusinvec(nploc))
      allocate (rvcosvec(nploc),rvsinvec(nploc))
      allocate (rwcosvec(nploc),rwsinvec(nploc))
      allocate (uvcosvec(nploc),uvsinvec(nploc))
      allocate (uwcosvec(nploc),uwsinvec(nploc))
      allocate (vwcosvec(nploc),vwsinvec(nploc))
      allocate (upcosvec(nploc),upsinvec(nploc))
      allocate (urcosvec(nploc),ursinvec(nploc))
      allocate (uscosvec(nploc),ussinvec(nploc))
      allocate (vpcosvec(nploc),vpsinvec(nploc))
      allocate (vscosvec(nploc),vssinvec(nploc))
      allocate (wpcosvec(nploc),wpsinvec(nploc))
      allocate (wscosvec(nploc),wssinvec(nploc))
      allocate (ut1cosvec(nploc),ut1sinvec(nploc))
      allocate (ut2cosvec(nploc),ut2sinvec(nploc))
      allocate (dphiduvec(nploc),dphidvvec(nploc))
      allocate (dphidwvec(nploc),dphidrvec(nploc))
      allocate (dphiddelvec(nploc))
      allocate (dphidsvec(nploc))
      allocate (duvec(nploc,3),dvvec(nploc,3),dwvec(nploc,3))
c
c
c     zero out force components on local frame-defining atoms
c
      frcxvec (1:nploc,:)  = 0.0_ti_p
      frcyvec (1:nploc,:)  = 0.0_ti_p
      frczvec (1:nploc,:)  = 0.0_ti_p
      frcxvec1(1:nploc,:)  = 0.0_ti_p
      frcyvec1(1:nploc,:)  = 0.0_ti_p
      frczvec1(1:nploc,:)  = 0.0_ti_p
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'torquevec_rec'

c
c     get the local frame type and the frame-defining atoms
c
      axetypvec(1:nploc)  = polaxe(ivec(1:nploc))
      mask     (1:nploc) = (axetypvec  (1:nploc) /= 'None') !mask for axetyp
      nploc1        = count(mask       (1:nploc))

      if (nploc1.eq.0) return ! All axetyp are 'None', so return
      ivec1 (1:nploc1) = pack (ivec (1:nploc)  , mask(1:nploc))
      iavec (1:nploc1) = zaxis(ivec1(1:nploc1))
      where (iavec(1:nploc1) > 0)
         ialocvec = locrec1(iavec)
      endwhere
      ibvec   (1:nploc1) = ipole  (ivec1(1:nploc1))
      iblocvec(1:nploc1) = locrec1(ibvec(1:nploc1))
      icvec   (1:nploc1) = xaxis  (ivec1(1:nploc1))
      where (icvec(1:nploc1) > 0)
            iclocvec = locrec1(icvec)
      endwhere
      idvec   (1:nploc1) = yaxis(ivec1(1:nploc1))
      where (idvec(1:nploc1) > 0)
            idlocvec = locrec1(idvec)
      endwhere
c
c     construct the three rotation axes for the local frame
c
      uvec(1:nploc1,1) = x(iavec(1:nploc1)) - x(ibvec(1:nploc1))
      uvec(1:nploc1,2) = y(iavec(1:nploc1)) - y(ibvec(1:nploc1))
      uvec(1:nploc1,3) = z(iavec(1:nploc1)) - z(ibvec(1:nploc1))
      usizvec(1:nploc1) = sqrt(  uvec(1:nploc1,1)**2
     &                         + uvec(1:nploc1,2)**2
     &                         + uvec(1:nploc1,3)**2
     &                        )

      where (axetypvec(1:nploc1) /= 'Z-Only')
           vvec(:,1) = x(icvec(:)) - x(ibvec(:))
           vvec(:,2) = y(icvec(:)) - y(ibvec(:))
           vvec(:,3) = z(icvec(:)) - z(ibvec(:))
           vsizvec   = sqrt(vvec(:,1)**2 + vvec(:,2)**2 + vvec(:,3)**2)
      elsewhere
           vvec(:,1) = merge (0.0_ti_p, 1.0_ti_p,
     &           (abs(uvec(:,1) / usizvec)) > 0.866_ti_p)
           vvec(:,2) = merge (1.0_ti_p, 0.0_ti_p,
     &           (abs(uvec(:,1) / usizvec)) > 0.866_ti_p)
           vvec(:,3) = 0.0_ti_p
           vsizvec   = 1.0_ti_p
      endwhere


      where (    axetypvec(1:nploc1) == 'Z-Bisect'
     &       .or.axetypvec(1:nploc1) == '3-Fold')
         wvec(:,1) = x(idvec) - x(ibvec)
         wvec(:,2) = y(idvec) - y(ibvec)
         wvec(:,3) = z(idvec) - z(ibvec)
      elsewhere 
         wvec(:,1) =  uvec(:,2) * vvec(:,3) - uvec(:,3) * vvec(:,2)
         wvec(:,2) =  uvec(:,3) * vvec(:,1) - uvec(:,1) * vvec(:,3)
         wvec(:,3) =  uvec(:,1) * vvec(:,2) - uvec(:,2) * vvec(:,1)
      endwhere
      wsizvec(1:nploc1) = sqrt(  wvec(1:nploc1,1)**2
     &                         + wvec(1:nploc1,2)**2
     &                         + wvec(1:nploc1,3)**2
     &                        )
      do j = 1, 3
         uvec(1:nploc1,j) = uvec(1:nploc1,j) / usizvec(1:nploc1)
         vvec(1:nploc1,j) = vvec(1:nploc1,j) / vsizvec(1:nploc1)
         wvec(1:nploc1,j) = wvec(1:nploc1,j) / wsizvec(1:nploc1)
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      where (axetypvec(1:nploc1) == 'Z-Bisect')
         rvec(:,1)  = vvec(:,1) + wvec(:,1)
         rvec(:,2)  = vvec(:,2) + wvec(:,2)
         rvec(:,3)  = vvec(:,3) + wvec(:,3)
         rsizvec = sqrt( rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2 )

         svec(:,1) =  uvec(:,2) * rvec(:,3) - uvec(:,3) * rvec(:,2)
         svec(:,2) =  uvec(:,3) * rvec(:,1) - uvec(:,1) * rvec(:,3)
         svec(:,3) =  uvec(:,1) * rvec(:,2) - uvec(:,2) * rvec(:,1)
         ssizvec = sqrt( svec(:,1)**2 + svec(:,2)**2 + svec(:,3)**2 )
         rvec(:,1) = rvec(:,1) / rsizvec
         svec(:,1) = svec(:,1) / ssizvec
         rvec(:,2) = rvec(:,2) / rsizvec
         svec(:,2) = svec(:,2) / ssizvec
         rvec(:,3) = rvec(:,3) / rsizvec
         svec(:,3) = svec(:,3) / ssizvec
      endwhere

c
c     find the perpendicular and angle for each pair of axes
c
         uvvec(1:nploc1,1)  =  vvec(1:nploc1,2) * uvec(1:nploc1,3)
     &                       - vvec(1:nploc1,3) * uvec(1:nploc1,2)
         uvvec(1:nploc1,2)  =  vvec(1:nploc1,3) * uvec(1:nploc1,1)
     &                       - vvec(1:nploc1,1) * uvec(1:nploc1,3)
         uvvec(1:nploc1,3)  =  vvec(1:nploc1,1) * uvec(1:nploc1,2)
     &                       - vvec(1:nploc1,2) * uvec(1:nploc1,1)

         uwvec(1:nploc1,1)  =  wvec(1:nploc1,2) * uvec(1:nploc1,3)
     &                       - wvec(1:nploc1,3) * uvec(1:nploc1,2)
         uwvec(1:nploc1,2)  =  wvec(1:nploc1,3) * uvec(1:nploc1,1)
     &                       - wvec(1:nploc1,1) * uvec(1:nploc1,3)
         uwvec(1:nploc1,3)  =  wvec(1:nploc1,1) * uvec(1:nploc1,2)
     &                       - wvec(1:nploc1,2) * uvec(1:nploc1,1)

         vwvec(1:nploc1,1)  =  wvec(1:nploc1,2) * vvec(1:nploc1,3)
     &                       - wvec(1:nploc1,3) * vvec(1:nploc1,2)
         vwvec(1:nploc1,2)  =  wvec(1:nploc1,3) * vvec(1:nploc1,1)
     &                       - wvec(1:nploc1,1) * vvec(1:nploc1,3)
         vwvec(1:nploc1,3)  =  wvec(1:nploc1,1) * vvec(1:nploc1,2)
     &                       - wvec(1:nploc1,2) * vvec(1:nploc1,1)

         uvsizvec(1:nploc1) = sqrt(  uvvec(1:nploc1,1)**2
     &                             + uvvec(1:nploc1,2)**2
     &                             + uvvec(1:nploc1,3)**2
     &                            )
         uwsizvec(1:nploc1) = sqrt(  uwvec(1:nploc1,1)**2
     &                             + uwvec(1:nploc1,2)**2
     &                             + uwvec(1:nploc1,3)**2
     &                            )
         vwsizvec(1:nploc1) = sqrt(  vwvec(1:nploc1,1)**2
     &                             + vwvec(1:nploc1,2)**2
     &                             + vwvec(1:nploc1,3)**2
     &                            )
      do j = 1, 3
         uvvec(1:nploc1,j) = uvvec(1:nploc1,j) / uvsizvec(1:nploc1)
         uwvec(1:nploc1,j) = uwvec(1:nploc1,j) / uwsizvec(1:nploc1)
         vwvec(1:nploc1,j) = vwvec(1:nploc1,j) / vwsizvec(1:nploc1)
      end do
      where (axetypvec(1:nploc1) == 'Z-Bisect')
         urvec(:,1)  =  rvec(:,2) * uvec(:,3) - rvec(:,3) * uvec(:,2)
         urvec(:,2)  =  rvec(:,3) * uvec(:,1) - rvec(:,1) * uvec(:,3)
         urvec(:,3)  =  rvec(:,1) * uvec(:,2) - rvec(:,2) * uvec(:,1)

         usvec(:,1)  =  svec(:,2) * uvec(:,3) - svec(:,3) * uvec(:,2)
         usvec(:,2)  =  svec(:,3) * uvec(:,1) - svec(:,1) * uvec(:,3)
         usvec(:,3)  =  svec(:,1) * uvec(:,2) - svec(:,2) * uvec(:,1)

         vsvec(:,1)  =  svec(:,2) * vvec(:,3) - svec(:,3) * vvec(:,2)
         vsvec(:,2)  =  svec(:,3) * vvec(:,1) - svec(:,1) * vvec(:,3)
         vsvec(:,3)  =  svec(:,1) * vvec(:,2) - svec(:,2) * vvec(:,1)

         wsvec(:,1)  =  svec(:,2) * wvec(:,3) - svec(:,3) * wvec(:,2)
         wsvec(:,2)  =  svec(:,3) * wvec(:,1) - svec(:,1) * wvec(:,3)
         wsvec(:,3)  =  svec(:,1) * wvec(:,2) - svec(:,2) * wvec(:,1)

         ursizvec = sqrt(urvec(:,1)**2 + urvec(:,2)**2 + urvec(:,3)**2)
         ussizvec = sqrt(usvec(:,1)**2 + usvec(:,2)**2 + usvec(:,3)**2)
         vssizvec = sqrt(vsvec(:,1)**2 + vsvec(:,2)**2 + vsvec(:,3)**2)
         wssizvec = sqrt(wsvec(:,1)**2 + wsvec(:,2)**2 + wsvec(:,3)**2)
         urvec(:,1)  = urvec(:,1) / ursizvec
         usvec(:,1)  = usvec(:,1) / ussizvec
         vsvec(:,1)  = vsvec(:,1) / vssizvec
         wsvec(:,1)  = wsvec(:,1) / wssizvec
         urvec(:,2)  = urvec(:,2) / ursizvec
         usvec(:,2)  = usvec(:,2) / ussizvec
         vsvec(:,2)  = vsvec(:,2) / vssizvec
         wsvec(:,2)  = wsvec(:,2) / wssizvec
         urvec(:,3)  = urvec(:,3) / ursizvec
         usvec(:,3)  = usvec(:,3) / ussizvec
         vsvec(:,3)  = vsvec(:,3) / vssizvec
         wsvec(:,3)  = wsvec(:,3) / wssizvec
      endwhere
c
c     get sine and cosine of angles between the rotation axes
c
      uvcosvec(1:nploc1) =  uvec(1:nploc1,1) * vvec(1:nploc1,1)
     &                    + uvec(1:nploc1,2) * vvec(1:nploc1,2)
     &                    + uvec(1:nploc1,3) * vvec(1:nploc1,3)
      uvsinvec(1:nploc1) =  sqrt(1.0_ti_p - uvcosvec(1:nploc1)**2)
      uwcosvec(1:nploc1) =  uvec(1:nploc1,1) * wvec(1:nploc1,1)
     &                    + uvec(1:nploc1,2) * wvec(1:nploc1,2)
     &                    + uvec(1:nploc1,3) * wvec(1:nploc1,3)
      uwsinvec(1:nploc1) =  sqrt(1.0_ti_p - uwcosvec(1:nploc1)**2)
      vwcosvec(1:nploc1) =  vvec(1:nploc1,1) * wvec(1:nploc1,1)
     &                    + vvec(1:nploc1,2) * wvec(1:nploc1,2)
     &                    + vvec(1:nploc1,3) * wvec(1:nploc1,3)
      vwsinvec(1:nploc1) =  sqrt(1.0_ti_p - vwcosvec(1:nploc1)**2)
c
      where (axetypvec(1:nploc1) == 'Z-Bisect')
         urcosvec =  uvec(:,1) * rvec(:,1) + uvec(:,2) * rvec(:,2)
     &             + uvec(:,3) * rvec(:,3)
         ursinvec =  sqrt(1.0_ti_p - urcosvec(1:nploc1)**2)
         uscosvec =  uvec(:,1) * svec(:,1) + uvec(:,2) * svec(:,2)
     &             + uvec(:,3) * svec(:,3)
         ussinvec =  sqrt(1.0_ti_p - uscosvec(1:nploc1)**2)
         vscosvec =  vvec(:,1) * svec(:,1) + vvec(:,2) * svec(:,2)
     &             + vvec(:,3) * svec(:,3)
         vssinvec =  sqrt(1.0_ti_p - vscosvec(1:nploc1)**2)
         wscosvec =  wvec(:,1) * svec(:,1) + wvec(:,2) * svec(:,2)
     &             + wvec(:,3) * svec(:,3)
         wssinvec =  sqrt(1.0_ti_p - wscosvec(1:nploc1)**2)
c
c     compute the projection of v and w onto the ru-plane
c
         t1vec(:,1)  =  vvec(:,1) - svec(:,1) * vscosvec
         t2vec(:,1)  =  wvec(:,1) - svec(:,1) * wscosvec
         t1vec(:,2)  =  vvec(:,2) - svec(:,2) * vscosvec
         t2vec(:,2)  =  wvec(:,2) - svec(:,2) * wscosvec
         t1vec(:,3)  =  vvec(:,3) - svec(:,3) * vscosvec
         t2vec(:,3)  =  wvec(:,3) - svec(:,3) * wscosvec

         t1sizvec = sqrt(t1vec(:,1)**2 + t1vec(:,2)**2 + t1vec(:,3)**2)
         t2sizvec = sqrt(t2vec(:,1)**2 + t2vec(:,2)**2 + t2vec(:,3)**2)
         t1vec(:,1)  = t1vec(:,1) / t1sizvec
         t1vec(:,2)  = t1vec(:,2) / t1sizvec
         t1vec(:,3)  = t1vec(:,3) / t1sizvec
         t2vec(:,1)  = t2vec(:,1) / t2sizvec
         t2vec(:,2)  = t2vec(:,2) / t2sizvec
         t2vec(:,3)  = t2vec(:,3) / t2sizvec

         ut1cosvec =  uvec(:,1) * t1vec(:,1) + uvec(:,2) * t1vec(:,2)
     &              + uvec(:,3) * t1vec(:,3)
         ut1sinvec =  sqrt(1.0_ti_p - ut1cosvec(1:nploc1)**2)
         ut2cosvec =  uvec(:,1) * t2vec(:,1) + uvec(:,2) * t2vec(:,2)
     &              + uvec(:,3) * t2vec(:,3)
         ut2sinvec =  sqrt(1.0_ti_p - ut2cosvec(1:nploc1)**2)
      endwhere

c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphiduvec(1:nploc1) = - trqvec(1:nploc1,1) * uvec(1:nploc1,1)
     &                      - trqvec(1:nploc1,2) * uvec(1:nploc1,2)
     &                      - trqvec(1:nploc1,3) * uvec(1:nploc1,3)
      dphidvvec(1:nploc1) = - trqvec(1:nploc1,1) * vvec(1:nploc1,1)
     &                      - trqvec(1:nploc1,2) * vvec(1:nploc1,2)
     &                      - trqvec(1:nploc1,3) * vvec(1:nploc1,3)
      dphidwvec(1:nploc1) = - trqvec(1:nploc1,1) * wvec(1:nploc1,1)
     &                      - trqvec(1:nploc1,2) * wvec(1:nploc1,2)
     &                      - trqvec(1:nploc1,3) * wvec(1:nploc1,3)

      where (axetypvec(1:nploc1) == 'Z-Bisect')
        dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &              - trqvec(:,3) * rvec(:,3)
        dphidsvec = - trqvec(:,1) * svec(:,1) - trqvec(:,2) * svec(:,2)
     &              - trqvec(:,3) * svec(:,3)
      endwhere

c
c     force distribution for the Z-Only local coordinate method
c
      where (axetypvec(1:nploc1) == 'Z-Only')
         duvec(:,1) =  uvvec (:,1)  * dphidvvec / (usizvec * uvsinvec )
     &               + uwvec (:,1)  * dphidwvec /  usizvec
         duvec(:,2) =  uvvec (:,2)  * dphidvvec / (usizvec * uvsinvec )
     &               + uwvec (:,2)  * dphidwvec /  usizvec
         duvec(:,3) =  uvvec (:,3)  * dphidvvec / (usizvec * uvsinvec )
     &               + uwvec (:,3)  * dphidwvec /  usizvec

         devec(1,ialocvec) =  devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) =  devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) =  devec(3,ialocvec) + duvec(:,3)

         devec(1,iblocvec) =  devec(1,iblocvec) - duvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - duvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - duvec(:,3)

         frczvec1(:,1) =  frczvec1(:,1) + duvec(:,1)
         frczvec1(:,2) =  frczvec1(:,2) + duvec(:,2)
         frczvec1(:,3) =  frczvec1(:,3) + duvec(:,3)
c
c     force distribution for the Z-then-X local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == 'Z-then-X')
         duvec(:,1) =   uvvec(:,1) * dphidvvec / (usizvec * uvsinvec)
     &                + uwvec(:,1) * dphidwvec /  usizvec
         duvec(:,2) =   uvvec(:,2) * dphidvvec / (usizvec * uvsinvec)
     &                + uwvec(:,2) * dphidwvec /  usizvec
         duvec(:,3) =   uvvec(:,3) * dphidvvec / (usizvec * uvsinvec)
     &                + uwvec(:,3) * dphidwvec /  usizvec

         dvvec(:,1) = - uvvec(:,1) * dphiduvec / (vsizvec* uvsinvec)
         dvvec(:,2) = - uvvec(:,2) * dphiduvec / (vsizvec* uvsinvec)
         dvvec(:,3) = - uvvec(:,3) * dphiduvec / (vsizvec* uvsinvec)

         devec(1,ialocvec) = devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) = devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) = devec(3,ialocvec) + duvec(:,3)

         devec(1,iclocvec) = devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) = devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) = devec(3,iclocvec) + dvvec(:,3)

         devec(1,iblocvec) = devec(1,iblocvec) - duvec(:,1) - dvvec(:,1)
         devec(2,iblocvec) = devec(2,iblocvec) - duvec(:,2) - dvvec(:,2)
         devec(3,iblocvec) = devec(3,iblocvec) - duvec(:,3) - dvvec(:,3)

         frczvec1(:,1)     =  frczvec1(:,1) + duvec(:,1)
         frczvec1(:,2)     =  frczvec1(:,2) + duvec(:,2)
         frczvec1(:,3)     =  frczvec1(:,3) + duvec(:,3)
         frcxvec1(:,1)     =  frcxvec1(:,1) + dvvec(:,1)
         frcxvec1(:,2)     =  frcxvec1(:,2) + dvvec(:,2)
         frcxvec1(:,3)     =  frcxvec1(:,3) + dvvec(:,3)
c
c     force distribution for the Bisector local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == 'Bisector')
         duvec(:,1) =      uvvec(:,1) * dphidvvec / (usizvec * uvsinvec)
     &          +  0.5_ti_p * uwvec(:,1) * dphidwvec /  usizvec
         duvec(:,2) =      uvvec(:,2) * dphidvvec / (usizvec * uvsinvec)
     &          +  0.5_ti_p * uwvec(:,2) * dphidwvec /  usizvec
         duvec(:,3) =      uvvec(:,3) * dphidvvec / (usizvec * uvsinvec)
     &          +  0.5_ti_p * uwvec(:,3) * dphidwvec /  usizvec
         dvvec(:,1) =   -  uvvec(:,1) * dphiduvec / (vsizvec * uvsinvec)
     &          +  0.5_ti_p * vwvec(:,1) * dphidwvec /  vsizvec
         dvvec(:,2) =   -  uvvec(:,2) * dphiduvec / (vsizvec * uvsinvec)
     &          +  0.5_ti_p * vwvec(:,2) * dphidwvec /  vsizvec
         dvvec(:,3) =   -  uvvec(:,3) * dphiduvec / (vsizvec * uvsinvec)
     &          +  0.5_ti_p * vwvec(:,3) * dphidwvec /  vsizvec

         devec(1,ialocvec) = devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) = devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) = devec(3,ialocvec) + duvec(:,3)

         devec(1,iclocvec) = devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) = devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) = devec(3,iclocvec) + dvvec(:,3)

         devec(1,iblocvec) = devec(1,iblocvec) - duvec(:,1) - dvvec(:,1)
         devec(2,iblocvec) = devec(2,iblocvec) - duvec(:,2) - dvvec(:,2)
         devec(3,iblocvec) = devec(3,iblocvec) - duvec(:,3) - dvvec(:,3)

         frczvec1(:,1)     = frczvec1(:,1)     + duvec(:,1)
         frczvec1(:,2)     = frczvec1(:,2)     + duvec(:,2)
         frczvec1(:,3)     = frczvec1(:,3)     + duvec(:,3)
         frcxvec1(:,1)     = frcxvec1(:,1)     + dvvec(:,1)
         frcxvec1(:,2)     = frcxvec1(:,2)     + dvvec(:,2)
         frcxvec1(:,3)     = frcxvec1(:,3)     + dvvec(:,3)
c
c     force distribution for the Z-Bisect local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == 'Z-Bisect')
         duvec(:,1) =  urvec(:,1) * dphidrvec / (usizvec * ursinvec)
     &               + usvec(:,1) * dphidsvec /  usizvec
         duvec(:,2) =  urvec(:,2) * dphidrvec / (usizvec * uvsinvec)
     &               + usvec(:,2) * dphidsvec /  usizvec
         duvec(:,3) =  urvec(:,3) * dphidrvec / (usizvec * uvsinvec)
     &               + usvec(:,3) * dphidsvec /  usizvec
         dvvec(:,1) =  (vssinvec * svec (:,1) - vscosvec   * t1vec(:,1))
     &               * dphiduvec / (vsizvec   * (ut1sinvec + ut2sinvec))
         dvvec(:,2) =  (vssinvec * svec (:,2) - vscosvec   * t1vec(:,2))
     &               * dphiduvec / (vsizvec   * (ut1sinvec + ut2sinvec))
         dvvec(:,3) =  (vssinvec * svec (:,3) - vscosvec   * t1vec(:,3))
     &               * dphiduvec / (vsizvec   * (ut1sinvec + ut2sinvec))
         dwvec(:,1) =  (wssinvec * svec (:,1) - wscosvec   * t2vec(:,1))
     &               * dphiduvec / (wsizvec   * (ut1sinvec + ut2sinvec))
         dwvec(:,2) =  (wssinvec * svec (:,2) - wscosvec   * t2vec(:,2))
     &               * dphiduvec / (wsizvec   * (ut1sinvec + ut2sinvec))
         dwvec(:,3) =  (wssinvec * svec (:,3) - wscosvec   * t2vec(:,3))
     &               * dphiduvec / (wsizvec   * (ut1sinvec + ut2sinvec))

         devec(1,ialocvec) = devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) = devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) = devec(3,ialocvec) + duvec(:,3)

         devec(1,iclocvec) = devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) = devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) = devec(3,iclocvec) + dvvec(:,3)

         devec(1,idlocvec) = devec(1,idlocvec) + dwvec(:,1)
         devec(2,idlocvec) = devec(2,idlocvec) + dwvec(:,2)
         devec(3,idlocvec) = devec(3,idlocvec) + dwvec(:,3)

         devec(1,iblocvec) = devec(1,iblocvec) - duvec(:,1) - dvvec(:,1)
     &                                         - dwvec(:,1)
         devec(2,iblocvec) = devec(2,iblocvec) - duvec(:,2) - dvvec(:,2)
     &                                         - dwvec(:,2)
         devec(3,iblocvec) = devec(3,iblocvec) - duvec(:,3) - dvvec(:,3)
     &                                         - dwvec(:,3)
         frczvec1(:,1)     =  frczvec1(:,1)    + duvec(:,1)
         frczvec1(:,2)     =  frczvec1(:,2)    + duvec(:,2)
         frczvec1(:,3)     =  frczvec1(:,3)    + duvec(:,3)

         frcxvec1(:,1)     =  frcxvec1(:,1)    + dvvec(:,1)
         frcxvec1(:,2)     =  frcxvec1(:,2)    + dvvec(:,2)
         frcxvec1(:,3)     =  frcxvec1(:,3)    + dvvec(:,3)

         frcyvec1(:,1)     =  frcyvec1(:,1)    + dwvec(:,1)
         frcyvec1(:,2)     =  frcyvec1(:,2)    + dwvec(:,2)
         frcyvec1(:,3)     =  frcyvec1(:,3)    + dwvec(:,3)
c
c     force distribution for the 3-Fold local coordinate method
c
      elsewhere (axetypvec(1:nploc1) == '3-Fold')
         pvec(:,1) = uvec(:,1) + vvec(:,1) + wvec(:,1)
         pvec(:,2) = uvec(:,2) + vvec(:,2) + wvec(:,2)
         pvec(:,3) = uvec(:,3) + vvec(:,3) + wvec(:,3)
         psizvec   = sqrt(pvec(:,1)**2 + pvec(:,2)**2 + pvec(:,3)**2)

         pvec(:,1) = pvec(:,1) / psizvec
         pvec(:,2) = pvec(:,2) / psizvec
         pvec(:,3) = pvec(:,3) / psizvec
         upcosvec  =  uvec(:,1) * pvec(:,1) + uvec(:,2) * pvec(:,2)
     &              + uvec(:,3) * pvec(:,3)
         vpcosvec  =  vvec(:,1) * pvec(:,1) + vvec(:,2) * pvec(:,2)
     &              + vvec(:,3) * pvec(:,3)
         wpcosvec  =  wvec(:,1) * pvec(:,1) + wvec(:,2) * pvec(:,2)
     &              + wvec(:,3) * pvec(:,3)
         upsinvec  =  sqrt(1.0_ti_p - upcosvec**2)
         vpsinvec  =  sqrt(1.0_ti_p - vpcosvec**2)
         wpsinvec  =  sqrt(1.0_ti_p - wpcosvec**2)

         rvec(:,1) = uvec(:,1) + vvec(:,1)
         rvec(:,2) = uvec(:,2) + vvec(:,2)
         rvec(:,3) = uvec(:,3) + vvec(:,3)
         rsizvec   = sqrt(rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2)

         rvec(:,1) = rvec(:,1) / rsizvec
         rvec(:,2) = rvec(:,2) / rsizvec
         rvec(:,3) = rvec(:,3) / rsizvec
         rwcosvec  =  rvec(:,1) * wvec(:,1) + rvec(:,2) * wvec(:,2)
     &              + rvec(:,3) * wvec(:,3)
         rwsinvec  = sqrt(1.0_ti_p - rwcosvec**2)
         dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &               - trqvec(:,3) * rvec(:,3)

         delvec(:,1)  =  rvec(:,2) * wvec(:,3) - rvec(:,3) * wvec(:,2)
         delvec(:,2)  =  rvec(:,3) * wvec(:,1) - rvec(:,1) * wvec(:,3)
         delvec(:,3)  =  rvec(:,1) * wvec(:,2) - rvec(:,2) * wvec(:,1)
         delsizvec    =  sqrt(  delvec(:,1)**2 + delvec(:,2)**2
     &                        + delvec(:,3)**2 )
         delvec(:,1)  =  delvec(:,1) / delsizvec
         delvec(:,2)  =  delvec(:,2) / delsizvec
         delvec(:,3)  =  delvec(:,3) / delsizvec
         dphiddelvec  = - trqvec(:,1) * delvec(:,1)
     &                  - trqvec(:,2) * delvec(:,2)
     &                  - trqvec(:,3) * delvec(:,3)

         epsvec(:,1) = delvec(:,2) * wvec(:,3) - delvec(:,3) * wvec(:,2)
         epsvec(:,2) = delvec(:,3) * wvec(:,1) - delvec(:,1) * wvec(:,3)
         epsvec(:,3) = delvec(:,1) * wvec(:,2) - delvec(:,2) * wvec(:,1)
         dwvec(:,1)  = delvec(:,1) * dphidrvec   / (wsizvec * rwsinvec)
     &              +  epsvec(:,1) * dphiddelvec * wpcosvec
     &                                           / (wsizvec * psizvec )
         dwvec(:,2) =  delvec(:,2) * dphidrvec   / (wsizvec * rwsinvec)
     &              +  epsvec(:,2) * dphiddelvec * wpcosvec
     &                                           / (wsizvec * psizvec )
         dwvec(:,3) =  delvec(:,3) * dphidrvec   / (wsizvec * rwsinvec)
     &              +  epsvec(:,3) * dphiddelvec * wpcosvec
     &                                           / (wsizvec * psizvec )
         devec(1,idlocvec) =  devec(1,idlocvec) + dwvec(:,1)
         devec(2,idlocvec) =  devec(2,idlocvec) + dwvec(:,2)
         devec(3,idlocvec) =  devec(3,idlocvec) + dwvec(:,3)
         devec(1,iblocvec) =  devec(1,iblocvec) - dwvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - dwvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - dwvec(:,3)
         frcyvec1(:,1)     =  frcyvec1  (:,1)   + dwvec(:,1)
         frcyvec1(:,2)     =  frcyvec1  (:,2)   + dwvec(:,2)
         frcyvec1(:,3)     =  frcyvec1  (:,3)   + dwvec(:,3)
         rvec(:,1)  =  vvec(:,1) + wvec(:,1)
         rvec(:,2)  =  vvec(:,2) + wvec(:,2)
         rvec(:,3)  =  vvec(:,3) + wvec(:,3)
         rsizvec = sqrt(rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2)

         rvec(:,1) =  rvec(:,1) / rsizvec
         rvec(:,2) =  rvec(:,2) / rsizvec
         rvec(:,3) =  rvec(:,3) / rsizvec
         rucosvec  =  rvec(:,1) * uvec(:,1) + rvec(:,2) * uvec(:,2)
     &              + rvec(:,3) * uvec(:,3)
         rusinvec  = sqrt(1.0_ti_p - rucosvec**2)

         dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &               - trqvec(:,3) * rvec(:,3)

         delvec(:,1) =  rvec(:,2) * uvec(:,3) - rvec(:,3) * uvec(:,2)
         delvec(:,2) =  rvec(:,3) * uvec(:,1) - rvec(:,1) * uvec(:,3)
         delvec(:,3) =  rvec(:,1) * uvec(:,2) - rvec(:,2) * uvec(:,1)

         delsizvec = sqrt(  delvec(:,1)**2 + delvec(:,2)**2
     &                    + delvec(:,3)**2 )
         delvec(:,1) =  delvec(:,1) / delsizvec
         delvec(:,2) =  delvec(:,2) / delsizvec
         delvec(:,3) =  delvec(:,3) / delsizvec
         dphiddelvec = - trqvec(:,1) * delvec(:,1)
     &                 - trqvec(:,2) * delvec(:,2)
     &                 - trqvec(:,3) * delvec(:,3)
         epsvec(:,1) = delvec(:,2) * uvec(:,3) - delvec(:,3) * uvec(:,2)
         epsvec(:,2) = delvec(:,3) * uvec(:,1) - delvec(:,1) * uvec(:,3)
         epsvec(:,3) = delvec(:,1) * uvec(:,2) - delvec(:,2) * uvec(:,1)

         duvec(:,1) =  delvec(:,1) * dphidrvec   / (usizvec * rusinvec)
     &               + epsvec(:,1) * dphiddelvec * upcosvec
     &                                           / (usizvec * psizvec )
         duvec(:,2) =  delvec(:,2) * dphidrvec   / (usizvec * rusinvec)
     &               + epsvec(:,2) * dphiddelvec * upcosvec
     &                                           / (usizvec * psizvec )
         duvec(:,3) =  delvec(:,3) * dphidrvec   / (usizvec * rusinvec)
     &               + epsvec(:,3) * dphiddelvec * upcosvec
     &                                           / (usizvec * psizvec )
         devec(1,ialocvec) =  devec(1,ialocvec) + duvec(:,1)
         devec(2,ialocvec) =  devec(2,ialocvec) + duvec(:,2)
         devec(3,ialocvec) =  devec(3,ialocvec) + duvec(:,3)

         devec(1,iblocvec) =  devec(1,iblocvec) - duvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - duvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - duvec(:,3)

         frczvec1(:,1)     =  frczvec1(:,1)     + duvec(:,1)
         frczvec1(:,2)     =  frczvec1(:,2)     + duvec(:,2)
         frczvec1(:,3)     =  frczvec1(:,3)     + duvec(:,3)
         rvec(:,1)  =  uvec(:,1) + wvec(:,1)
         rvec(:,2)  =  uvec(:,2) + wvec(:,2)
         rvec(:,3)  =  uvec(:,3) + wvec(:,3)
         rsizvec = sqrt( rvec(:,1)**2 + rvec(:,2)**2 + rvec(:,3)**2 )

         rvec(:,1) =  rvec(:,1) / rsizvec
         rvec(:,2) =  rvec(:,2) / rsizvec
         rvec(:,3) =  rvec(:,3) / rsizvec
         rvcosvec  =  rvec(:,1) * vvec(:,1) + rvec(:,2) * vvec(:,2)
     &              + rvec(:,3) * vvec(:,3)
         rvsinvec = sqrt(1.0_ti_p - rvcosvec**2)
         dphidrvec = - trqvec(:,1) * rvec(:,1) - trqvec(:,2) * rvec(:,2)
     &               - trqvec(:,3) * rvec(:,3)

         delvec(:,1) =  rvec(:,2) * vvec(:,3) - rvec(:,3) * vvec(:,2)
         delvec(:,2) =  rvec(:,3) * vvec(:,1) - rvec(:,1) * vvec(:,3)
         delvec(:,3) =  rvec(:,1) * vvec(:,2) - rvec(:,2) * vvec(:,1)

         delsizvec = sqrt(  delvec(:,1)**2 + delvec(:,2)**2
     &                    + delvec(:,3)**2 )
        delvec(:,1) =  delvec(:,1) / delsizvec
        delvec(:,2) =  delvec(:,2) / delsizvec
        delvec(:,3) =  delvec(:,3) / delsizvec

        dphiddelvec = - trqvec(:,1) * delvec(:,1)
     &                - trqvec(:,2) * delvec(:,2)
     &                - trqvec(:,3) * delvec(:,3)
         epsvec(:,1) = delvec(:,2) * vvec(:,3) - delvec(:,3) * vvec(:,2)
         epsvec(:,2) = delvec(:,3) * vvec(:,1) - delvec(:,1) * vvec(:,3)
         epsvec(:,3) = delvec(:,1) * vvec(:,2) - delvec(:,2) * vvec(:,1)

         dvvec(:,1)  =  delvec(:,1) * dphidrvec   / (vsizvec * rvsinvec)
     &                + epsvec(:,1) * dphiddelvec * vpcosvec
     &                                            / (vsizvec * psizvec)
         dvvec(:,2)  =  delvec(:,2) * dphidrvec   / (vsizvec * rvsinvec)
     &                + epsvec(:,2) * dphiddelvec * vpcosvec
     &                                            / (vsizvec * psizvec)
         dvvec(:,3)  =  delvec(:,3) * dphidrvec   / (vsizvec * rvsinvec)
     &                + epsvec(:,3) * dphiddelvec * vpcosvec
     &                                            / (vsizvec * psizvec)
         devec(1,iclocvec) =  devec(1,iclocvec) + dvvec(:,1)
         devec(2,iclocvec) =  devec(2,iclocvec) + dvvec(:,2)
         devec(3,iclocvec) =  devec(3,iclocvec) + dvvec(:,3)

         devec(1,iblocvec) =  devec(1,iblocvec) - dvvec(:,1)
         devec(2,iblocvec) =  devec(2,iblocvec) - dvvec(:,2)
         devec(3,iblocvec) =  devec(3,iblocvec) - dvvec(:,3)
         frcxvec1(:,1)     =  frcxvec1(:,1)     + dvvec(:,1)
         frcxvec1(:,2)     =  frcxvec1(:,2)     + dvvec(:,2)
         frcxvec1(:,3)     =  frcxvec1(:,3)     + dvvec(:,3)
      endwhere


cnew
cnew  send the results back in the original arrays
      do i=1,3
         frcxvec(1:nploc,i) =
     &                unpack(frcxvec1(1:nploc1,i), mask(1:nploc), 0.0_ti_p)
         frcyvec(1:nploc,i) =
     &                unpack(frcyvec1(1:nploc1,i), mask(1:nploc), 0.0_ti_p)
         frczvec(1:nploc,i) =
     &                unpack(frczvec1(1:nploc1,i), mask(1:nploc), 0.0_ti_p)
      enddo

      deallocate (usizvec)
      deallocate (uvec)
      deallocate (ssizvec)
      deallocate (svec)
      deallocate (rsizvec)
      deallocate (rvec)
      deallocate (idlocvec)
      deallocate (iclocvec)
      deallocate (iblocvec)
      deallocate (ialocvec)
      deallocate (idvec)
      deallocate (icvec)
      deallocate (ibvec)
      deallocate (iavec)
      deallocate (mask)
      deallocate (ivec1)
      deallocate (uvvec)
      deallocate (uvsizvec)
      deallocate (vvec)
      deallocate (vsizvec)
      deallocate (wvec)
      deallocate (wsizvec)
      deallocate (axetypvec)
      deallocate (uwvec,uwsizvec)
      deallocate (vwvec,vwsizvec)
      deallocate (usvec,ussizvec)
      deallocate (vsvec,vssizvec)
      deallocate (wsvec,wssizvec)
      deallocate (urvec,ursizvec)
      deallocate (t1vec,t1sizvec)
      deallocate (t2vec,t2sizvec)
      deallocate (rucosvec,rusinvec)
      deallocate (rvcosvec,rvsinvec)
      deallocate (rwcosvec,rwsinvec)
      deallocate (uvcosvec,uvsinvec)
      deallocate (uwcosvec,uwsinvec)
      deallocate (vwcosvec,vwsinvec)
      deallocate (upcosvec,upsinvec)
      deallocate (urcosvec,ursinvec)
      deallocate (uscosvec,ussinvec)
      deallocate (vpcosvec,vpsinvec)
      deallocate (vscosvec,vssinvec)
      deallocate (wpcosvec,wpsinvec)
      deallocate (wscosvec,wssinvec)
      deallocate (ut1cosvec,ut1sinvec)
      deallocate (ut2cosvec,ut2sinvec)
      deallocate (dphiduvec,dphidvvec)
      deallocate (dphidwvec,dphidrvec)
      deallocate (dphiddelvec)
      deallocate (dphidsvec)
      deallocate (duvec,dvvec,dwvec)

      return
      end
