c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c
c     "bspline_fill_site" finds B-spline coefficients and derivatives
c     for PME i-th atomic sites along the fractional coordinate axes
c
c
#include "tinker_precision.h"
      subroutine bspline_fill_site(i,impi)
!$acc routine(bsplgen) seq
      use atoms
      use boxes
      use pme
      use tinheader
      implicit none
      integer i,ifr,impi
      real(t_p) xi,yi,zi
      real(t_p) w,fr,eps
c
c     offset used to shift sites off exact lattice bounds
c
#if defined(SINGLE) || defined(MIXED)
      eps = 1.0d-4
#else
      eps = 1.0d-8
#endif
c
c     get the b-spline coefficients for the i-th atomic site
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
      w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
      fr = real(nfft1,t_p) * (w-anint(w)+0.5_ti_p)
      ifr = int(fr-eps)
      w = fr - real(ifr,t_p)
      igrid(1,i) = ifr - bsorder
      call bsplgen (w,impi,thetai1(1,1,impi))
      w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
      fr = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
      ifr = int(fr-eps)
      w = fr - real(ifr,t_p)
      igrid(2,i) = ifr - bsorder
      call bsplgen (w,impi,thetai2(1,1,impi))
      w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
      fr = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
      ifr = int(fr-eps)
      w = fr - real(ifr,t_p)
      igrid(3,i) = ifr - bsorder
      call bsplgen (w,impi,thetai3(1,1,impi))
      return
      end

      subroutine AssociateThetai_p
      use pme
      use potent
      use domdec
      implicit none
      integer tsize

      tsize = size(thetai1,dim=3)
      if (use_charge) then
         thetai1_p(1:2,1:bsorder,1:2*tsize)=> thetai1(:,:,:)
         thetai2_p(1:2,1:bsorder,1:2*tsize)=> thetai2(:,:,:)
         thetai3_p(1:2,1:bsorder,1:2*tsize)=> thetai3(:,:,:)
      else
         thetai1_p(1:4,1:bsorder,1:1*tsize)=> thetai1(:,:,:)
         thetai2_p(1:4,1:bsorder,1:1*tsize)=> thetai2(:,:,:)
         thetai3_p(1:4,1:bsorder,1:1*tsize)=> thetai3(:,:,:)
      end if
      end subroutine
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bsplgen  --  B-spline coefficients for an atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bsplgen" gets B-spline coefficients and derivatives for
c     a single PME atomic site along a particular direction
c
c
      subroutine bsplgen (w,isite,thetai)
!$acc routine
      use tinheader
      use pme,only:bsorder,maxorder
      use potent,only:use_mpole,use_polar
      implicit none
      integer i,j,k
      integer isite
      integer level
      real(t_p) w,denom
      real(t_p) thetai(4,bsorder)
      real(t_p) temp(maxorder,maxorder)
c
c     set B-spline depth for partial charges or multipoles
c
      level = 2
      if (use_mpole .or. use_polar)  level = 4
c
c     initialization to get to 2nd order recursion
c
      temp(2,2) = w
      temp(2,1) = 1.0_ti_p - w
c
c     perform one pass to get to 3rd order recursion
c
      temp(3,3) = 0.5_ti_p * w * temp(2,2)
      temp(3,2) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2,1)+
     &                        (2.0_ti_p-w)*temp(2,2))
      temp(3,1) = 0.5_ti_p * (1.0_ti_p-w) * temp(2,1)
c
c     compute standard B-spline recursion to desired order
c
!$acc loop seq
      do i = 4, bsorder
         k = i - 1
         denom = 1.0_ti_p / real(k,t_p)
         temp(i,i) = denom * w * temp(k,k)
         do j = 1, i-2
            temp(i,i-j) = denom * ((w+real(j,t_p))*temp(k,i-j-1)
     &                           +(real(i-j,t_p)-w)*temp(k,i-j))
         end do
         temp(i,1) = denom * (1.0_ti_p-w) * temp(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      temp(k,bsorder) = temp(k,bsorder-1)
!$acc loop seq
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
c
c     get coefficients for the B-spline second derivative
c
      if (level .eq. 4) then
         k = bsorder - 2
         temp(k,bsorder-1) = temp(k,bsorder-2)
         do i = bsorder-2, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
         temp(k,bsorder) = temp(k,bsorder-1)
         do i = bsorder-1, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
c
c     get coefficients for the B-spline third derivative
c
         k = bsorder - 3
         temp(k,bsorder-2) = temp(k,bsorder-3)
         do i = bsorder-3, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
         temp(k,bsorder-1) = temp(k,bsorder-2)
         do i = bsorder-2, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
         temp(k,bsorder) = temp(k,bsorder-1)
         do i = bsorder-1, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
      end if
c
c     copy coefficients from temporary to permanent storage
c
!$acc loop seq
      do i = 1, bsorder
!$acc loop seq
          do j = 1, level
            thetai(j,i) = temp(bsorder-j+1,i)
          end do
      end do
      return
      end

      ! bslgen special version 
      ! To be used only with point charge forcefield (use_charge)
      subroutine bsplgen_chg (w,isite,thetai)
!$acc routine
      use tinheader,only: ti_p
      use pme      ,only: bsorder,maxorder
      use potent   ,only: use_charge
      implicit none
      integer i,j,k
      integer isite
      integer,parameter:: level=2
      real(t_p) w,denom
      real(t_p) thetai(level,bsorder,*)
      real(t_p) temp(maxorder,maxorder)

#ifdef TINKER_DEBUG
!$acc routine(fatal_acc)
      if (.not.use_charge) then
         print*,"FATAL ERROR !! bsplgen routine is specific to point",
     &   "charge"
         call fatal_acc
      end if
#endif
c
c     initialization to get to 2nd order recursion
c
      temp(2,2) = w
      temp(2,1) = 1.0_ti_p - w
c
c     perform one pass to get to 3rd order recursion
c
      temp(3,3) = 0.5_ti_p * w * temp(2,2)
      temp(3,2) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2,1)+
     &                        (2.0_ti_p-w)*temp(2,2))
      temp(3,1) = 0.5_ti_p * (1.0_ti_p-w) * temp(2,1)
c
c     compute standard B-spline recursion to desired order
c
!$acc loop seq
      do i = 4, bsorder
         k = i - 1
         denom = 1.0_ti_p / real(k,t_p)
         temp(i,i) = denom * w * temp(k,k)
         do j = 1, i-2
            temp(i,i-j) = denom * ((w+real(j,t_p))*temp(k,i-j-1)
     &                           +(real(i-j,t_p)-w)*temp(k,i-j))
         end do
         temp(i,1) = denom * (1.0_ti_p-w) * temp(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      temp(k,bsorder) = temp(k,bsorder-1)
!$acc loop seq
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
c
c     copy coefficients from temporary to permanent storage
c
!$acc loop seq
      do i = 1, bsorder
!$acc loop seq
         do j = 1, level
            thetai(j,i,isite) = temp(bsorder-j+1,i)
         end do
      end do
      end
c
c
c       "grid_mpole_site" places the i-th fractional atomic multipole onto
c       the particle mesh Ewald grid
c
c
      subroutine grid_mpole_site(isite,impi,fmp)
      use chunks
      use domdec
      use fft
      use pme
      use potent
      use sizes
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,m,impi,rankloc
      integer iproc,proc
      integer ii,jj,kk
      integer isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer nearpt(3)
      integer abound(6)
      real(t_p) v0,u0,t0
      real(t_p) v1,u1,t1
      real(t_p) v2,u2,t2
      real(t_p) term,term0,term1,term2
      real(t_p) fmp(10)
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
c
c       put the permanent multipole moments onto the grid
c
      iatm = isite
      nearpt(1) = igrid(1,iatm) + grdoff
      nearpt(2) = igrid(2,iatm) + grdoff
      nearpt(3) = igrid(3,iatm) + grdoff
      abound(1) = nearpt(1) - nlpts
      abound(2) = nearpt(1) + nrpts
      abound(3) = nearpt(2) - nlpts
      abound(4) = nearpt(2) + nrpts
      abound(5) = nearpt(3) - nlpts
      abound(6) = nearpt(3) + nrpts
      call adjust (offsetx,nfft1,1,abound(1),
     &               abound(2),1,nfft1)
      call adjust (offsety,nfft2,1,abound(3),
     &               abound(4),1,nfft2)
      call adjust (offsetz,nfft3,1,abound(5),
     &               abound(6),1,nfft3)
      do kk = abound(5), abound(6)
         k = kk
         m = k + offsetz
         if (k .lt. 1)  k = k + nfft3
         v0 = thetai3(1,m,impi)
         v1 = thetai3(2,m,impi)
         v2 = thetai3(3,m,impi)
         do jj = abound(3), abound(4)
            j = jj
            m = j + offsety
            if (j .lt. 1)  j = j + nfft2
            u0 = thetai2(1,m,impi)
            u1 = thetai2(2,m,impi)
            u2 = thetai2(3,m,impi)
            term0 = fmp(1)*u0*v0 + fmp(3)*u1*v0
     &            + fmp(4)*u0*v1 + fmp(6)*u2*v0
     &            + fmp(7)*u0*v2 + fmp(10)*u1*v1
            term1 = fmp(2)*u0*v0 + fmp(8)*u1*v0
     &                 + fmp(9)*u0*v1
            term2 = fmp(5) * u0 * v0
            do ii = abound(1), abound(2)
               i = ii
               m = i + offsetx
               if (i .lt. 1)  i = i + nfft1
               t0 = thetai1(1,m,impi)
               t1 = thetai1(2,m,impi)
               t2 = thetai1(3,m,impi)
c
               kstart = kstart1(rankloc+1)
               kend = kend1(rankloc+1)
               jstart = jstart1(rankloc+1)
               jend = jend1(rankloc+1)
               istart = istart1(rankloc+1)
               iend = iend1(rankloc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 term = term0*t0 + term1*t1 + term2*t2
!$acc atomic update
                 qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1) =
     $             qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)
     $             + term
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   term = term0*t0 + term1*t1 + term2*t2
!$acc atomic update
                   qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $               iproc+1) = qgridin_2d(1,i-istart+1,j-jstart+1,
     $               k-kstart+1,iproc+1) + term
                   goto 10
                 end if
               end do
 10            continue
            end do
         end do
      end do
      return
      end
c
c
c
c       "grid_uind_site" places the fractional induced dipoles onto the
c       particle mesh Ewald grid
c
c
      subroutine grid_uind_site(isite,impi,fuind,fuinp,qgrid2loc)
      use chunks
      use domdec
      use fft
      use pme
      use potent
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,m,impi
      integer iproc,proc
      integer ii,jj,kk
      integer isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer nearpt(3)
      integer abound(6)
      real(t_p) v0,u0,t0
      real(t_p) v1,u1,t1
      real(t_p) term
      real(t_p) term01,term11
      real(t_p) term02,term12
      real(t_p) fuind(3)
      real(t_p) fuinp(3)
      real(t_p) qgrid2loc(2,n1mpimax,n2mpimax,n3mpimax,
     &                      nrec_send+1)
c
c       put the induced dipole moments onto the grid
c
      iatm = isite
      nearpt(1) = igrid(1,iatm) + grdoff
      nearpt(2) = igrid(2,iatm) + grdoff
      nearpt(3) = igrid(3,iatm) + grdoff
      abound(1) = nearpt(1) - nlpts
      abound(2) = nearpt(1) + nrpts
      abound(3) = nearpt(2) - nlpts
      abound(4) = nearpt(2) + nrpts
      abound(5) = nearpt(3) - nlpts
      abound(6) = nearpt(3) + nrpts
      call adjust (offsetx,nfft1,1,abound(1),
     &               abound(2),1,nfft1)
      call adjust (offsety,nfft2,1,abound(3),
     &               abound(4),1,nfft2)
      call adjust (offsetz,nfft3,1,abound(5),
     &               abound(6),1,nfft3)
      do kk = abound(5), abound(6)
         k  = kk
         m  = k + offsetz
         if (k .lt. 1)  k = k + nfft3
         v0 = thetai3(1,m,impi)
         v1 = thetai3(2,m,impi)
         do jj = abound(3), abound(4)
            j  = jj
            m  = j + offsety
            if (j .lt. 1)  j = j + nfft2
            u0 = thetai2(1,m,impi)
            u1 = thetai2(2,m,impi)
            term01 = fuind(2)*u1*v0
     &                  + fuind(3)*u0*v1
            term11 = fuind(1)*u0*v0
            term02 = fuinp(2)*u1*v0
     &                  + fuinp(3)*u0*v1
            term12 = fuinp(1)*u0*v0
            do ii = abound(1), abound(2)
               i  = ii
               m  = i + offsetx
               if (i .lt. 1)  i = i + nfft1
               t0 = thetai1(1,m,impi)
               t1 = thetai1(2,m,impi)
c
               if (use_pmecore) then
                 kstart = kstart1(rank_bis+1)
                 kend   = kend1(rank_bis+1)
                 jstart = jstart1(rank_bis+1)
                 jend   = jend1(rank_bis+1)
                 istart = istart1(rank_bis+1)
                 iend   = iend1(rank_bis+1)
               else
                 kstart = kstart1(rank+1)
                 kend   = kend1(rank+1)
                 jstart = jstart1(rank+1)
                 jend   = jend1(rank+1)
                 istart = istart1(rank+1)
                 iend   = iend1(rank+1)
               end if
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 term = term01*t0 + term11*t1
!$acc atomic update
                 qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,1)
     $         = qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,1)
     $         + term
                 term = term02*t0 + term12*t1 
!$acc atomic update
                 qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,1)
     $         = qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,1)
     $         + term
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc   = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend   = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend   = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend   = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $               ((j.ge.jstart).and.(j.le.jend)).and.
     $               ((i.ge.istart).and.(i.le.iend))) then
                   term = term01*t0 + term11*t1
!$acc atomic update
                   qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     $           = qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     $           + term
                   term = term02*t0 + term12*t1
!$acc atomic update
                   qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     $           = qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     $           + term
                   goto 10
                 end if
               end do
 10            continue
c
            end do
         end do
      end do
      return
      end
c
c       #################################################################
c       ##                                                             ##
c       ##  subroutine adjust  --  alter site bounds for the PME grid  ##
c       ##                                                             ##
c       #################################################################
c
c
c       "adjust" modifies site bounds on the PME grid and returns
c       an offset into the B-spline coefficient arrays
c
c
      subroutine adjust (offset,nfft,nchk,amin,amax,cmin,cmax)
!$acc routine seq
      implicit none
      integer offset
      integer nfft,nchk
      integer amin,amax
      integer cmin,cmax
c
c
c     modify grid offset and bounds for site at edge of chunk
c
      offset = 0
      if (nchk .ne. 1) then
         if (amin.lt.cmin .or. amax.gt.cmax) then
            if (amin.lt.1 .or. amax.gt.nfft) then
               if (cmin .eq. 1) then
                  offset = 1 - amin
                  amin = 1
               else if (cmax .eq. nfft) then
                  amax = nfft
                  amin = amin + nfft
               end if
            else
               if (cmin .gt. amin) then
                  offset = cmin - amin
                  amin = cmin
               else
                  amax = cmax
               end if
            end if
         end if
      end if
      offset = offset + 1 - amin
      return
      end
c
c
c       "fphi_mpole_site" extracts the permanent multipole potential on the i-th site from
c       the particle mesh Ewald grid
c
c
      subroutine fphi_mpole_site(isite,impi)
      use domdec
      use fft
      use pme
      use potent
      use tinheader ,only:ti_p,re_p
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,impi,rankloc
      integer iproc,proc
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real(t_p) v0,v1,v2,v3
      real(t_p) u0,u1,u2,u3
      real(t_p) t0,t1,t2,t3,tq
      real(t_p) tu00,tu10,tu01,tu20,tu11
      real(t_p) tu02,tu21,tu12,tu30,tu03
      real(t_p) tuv000,tuv100,tuv010,tuv001
      real(t_p) tuv200,tuv020,tuv002,tuv110
      real(t_p) tuv101,tuv011,tuv300,tuv030
      real(t_p) tuv003,tuv210,tuv201,tuv120
      real(t_p) tuv021,tuv102,tuv012,tuv111
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
c
      iatm = isite
      igrd0 = igrid(1,iatm)
      jgrd0 = igrid(2,iatm)
      kgrd0 = igrid(3,iatm)
      tuv000 = 0.0_ti_p
      tuv001 = 0.0_ti_p
      tuv010 = 0.0_ti_p
      tuv100 = 0.0_ti_p
      tuv200 = 0.0_ti_p
      tuv020 = 0.0_ti_p
      tuv002 = 0.0_ti_p
      tuv110 = 0.0_ti_p
      tuv101 = 0.0_ti_p
      tuv011 = 0.0_ti_p
      tuv300 = 0.0_ti_p
      tuv030 = 0.0_ti_p
      tuv003 = 0.0_ti_p
      tuv210 = 0.0_ti_p
      tuv201 = 0.0_ti_p
      tuv120 = 0.0_ti_p
      tuv021 = 0.0_ti_p
      tuv102 = 0.0_ti_p
      tuv012 = 0.0_ti_p
      tuv111 = 0.0_ti_p
      k0 = kgrd0
      do it3 = 1, bsorder
         k0 = k0 + 1
c        k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         k = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
         v0 = thetai3(1,it3,impi)
         v1 = thetai3(2,it3,impi)
         v2 = thetai3(3,it3,impi)
         v3 = thetai3(4,it3,impi)
         tu00 = 0.0_ti_p
         tu10 = 0.0_ti_p
         tu01 = 0.0_ti_p
         tu20 = 0.0_ti_p
         tu11 = 0.0_ti_p
         tu02 = 0.0_ti_p
         tu30 = 0.0_ti_p
         tu21 = 0.0_ti_p
         tu12 = 0.0_ti_p
         tu03 = 0.0_ti_p
         j0 = jgrd0
         do it2 = 1, bsorder
            j0 = j0 + 1
c           j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            j = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
            u0 = thetai2(1,it2,impi)
            u1 = thetai2(2,it2,impi)
            u2 = thetai2(3,it2,impi)
            u3 = thetai2(4,it2,impi)
            t0 = 0.0_ti_p
            t1 = 0.0_ti_p
            t2 = 0.0_ti_p
            t3 = 0.0_ti_p
            i0 = igrd0
            do it1 = 1, bsorder
               i0 = i0 + 1
c              i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
               i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
c
                tq = 0.0_ti_p
                kstart = kstart1(rankloc+1)
                kend = kend1(rankloc+1)
                jstart = jstart1(rankloc+1)
                jend = jend1(rankloc+1)
                istart = istart1(rankloc+1)
                iend = iend1(rankloc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $              iproc+1)
                   goto 10
                 end if
               end do
 10            continue
c
               t0 = t0 + tq*thetai1(1,it1,impi)
               t1 = t1 + tq*thetai1(2,it1,impi)
               t2 = t2 + tq*thetai1(3,it1,impi)
               t3 = t3 + tq*thetai1(4,it1,impi)
            end do
            tu00 = tu00 + t0*u0
            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
            tu20 = tu20 + t2*u0
            tu11 = tu11 + t1*u1
            tu02 = tu02 + t0*u2
            tu30 = tu30 + t3*u0
            tu21 = tu21 + t2*u1
            tu12 = tu12 + t1*u2
            tu03 = tu03 + t0*u3
         end do
         tuv000 = tuv000 + tu00*v0
         tuv100 = tuv100 + tu10*v0
         tuv010 = tuv010 + tu01*v0
         tuv001 = tuv001 + tu00*v1
         tuv200 = tuv200 + tu20*v0
         tuv020 = tuv020 + tu02*v0
         tuv002 = tuv002 + tu00*v2
         tuv110 = tuv110 + tu11*v0
         tuv101 = tuv101 + tu10*v1
         tuv011 = tuv011 + tu01*v1
         tuv300 = tuv300 + tu30*v0
         tuv030 = tuv030 + tu03*v0
         tuv003 = tuv003 + tu00*v3
         tuv210 = tuv210 + tu21*v0
         tuv201 = tuv201 + tu20*v1
         tuv120 = tuv120 + tu12*v0
         tuv021 = tuv021 + tu02*v1
         tuv102 = tuv102 + tu10*v2
         tuv012 = tuv012 + tu01*v2
         tuv111 = tuv111 + tu11*v1
      end do
      fphirec(1,impi) = tuv000
      fphirec(2,impi) = tuv100
      fphirec(3,impi) = tuv010
      fphirec(4,impi) = tuv001
      fphirec(5,impi) = tuv200
      fphirec(6,impi) = tuv020
      fphirec(7,impi) = tuv002
      fphirec(8,impi) = tuv110
      fphirec(9,impi) = tuv101
      fphirec(10,impi) = tuv011
      fphirec(11,impi) = tuv300
      fphirec(12,impi) = tuv030
      fphirec(13,impi) = tuv003
      fphirec(14,impi) = tuv210
      fphirec(15,impi) = tuv201
      fphirec(16,impi) = tuv120
      fphirec(17,impi) = tuv021
      fphirec(18,impi) = tuv102
      fphirec(19,impi) = tuv012
      fphirec(20,impi) = tuv111
      return
      end
c
c     "fphi_uind_site" extracts the induced dipole potential at the i-th site from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uind_site(isite,impi,fdip_phi1,fdip_phi2,
     $  fdip_sum_phi)
      use domdec
      use fft
      use pme
      use potent
      use tinheader ,only:ti_p,re_p
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,impi
      integer iproc,proc
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real(t_p) v0,v1,v2,v3
      real(t_p) u0,u1,u2,u3
      real(t_p) t0,t1,t2,t3
      real(t_p) t0_1,t0_2,t1_1,t1_2
      real(t_p) t2_1,t2_2,tq_1,tq_2
      real(t_p) tu00,tu10,tu01,tu20,tu11
      real(t_p) tu02,tu30,tu21,tu12,tu03
      real(t_p) tu00_1,tu01_1,tu10_1
      real(t_p) tu00_2,tu01_2,tu10_2
      real(t_p) tu20_1,tu11_1,tu02_1
      real(t_p) tu20_2,tu11_2,tu02_2
      real(t_p) tuv100_1,tuv010_1,tuv001_1
      real(t_p) tuv100_2,tuv010_2,tuv001_2
      real(t_p) tuv200_1,tuv020_1,tuv002_1
      real(t_p) tuv110_1,tuv101_1,tuv011_1
      real(t_p) tuv200_2,tuv020_2,tuv002_2
      real(t_p) tuv110_2,tuv101_2,tuv011_2
      real(t_p) tuv000,tuv100,tuv010,tuv001
      real(t_p) tuv200,tuv020,tuv002,tuv110
      real(t_p) tuv101,tuv011,tuv300,tuv030
      real(t_p) tuv003,tuv210,tuv201,tuv120
      real(t_p) tuv021,tuv102,tuv012,tuv111
      real(t_p) fdip_phi1(10)
      real(t_p) fdip_phi2(10)
      real(t_p) fdip_sum_phi(20)
c
      iatm  = isite
      igrd0 = igrid(1,iatm)
      jgrd0 = igrid(2,iatm)
      kgrd0 = igrid(3,iatm)
      tuv100_1 = 0.0_ti_p
      tuv010_1 = 0.0_ti_p
      tuv001_1 = 0.0_ti_p
      tuv200_1 = 0.0_ti_p
      tuv020_1 = 0.0_ti_p
      tuv002_1 = 0.0_ti_p
      tuv110_1 = 0.0_ti_p
      tuv101_1 = 0.0_ti_p
      tuv011_1 = 0.0_ti_p
      tuv100_2 = 0.0_ti_p
      tuv010_2 = 0.0_ti_p
      tuv001_2 = 0.0_ti_p
      tuv200_2 = 0.0_ti_p
      tuv020_2 = 0.0_ti_p
      tuv002_2 = 0.0_ti_p
      tuv110_2 = 0.0_ti_p
      tuv101_2 = 0.0_ti_p
      tuv011_2 = 0.0_ti_p
      tuv000 = 0.0_ti_p
      tuv001 = 0.0_ti_p
      tuv010 = 0.0_ti_p
      tuv100 = 0.0_ti_p
      tuv200 = 0.0_ti_p
      tuv020 = 0.0_ti_p
      tuv002 = 0.0_ti_p
      tuv110 = 0.0_ti_p
      tuv101 = 0.0_ti_p
      tuv011 = 0.0_ti_p
      tuv300 = 0.0_ti_p
      tuv030 = 0.0_ti_p
      tuv003 = 0.0_ti_p
      tuv210 = 0.0_ti_p
      tuv201 = 0.0_ti_p
      tuv120 = 0.0_ti_p
      tuv021 = 0.0_ti_p
      tuv102 = 0.0_ti_p
      tuv012 = 0.0_ti_p
      tuv111 = 0.0_ti_p
      k0 = kgrd0
      do it3 = 1, bsorder
         k0 = k0 + 1
c        k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         k = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
         v0 = thetai3(1,it3,impi)
         v1 = thetai3(2,it3,impi)
         v2 = thetai3(3,it3,impi)
         v3 = thetai3(4,it3,impi)
         tu00_1 = 0.0_ti_p
         tu01_1 = 0.0_ti_p
         tu10_1 = 0.0_ti_p
         tu20_1 = 0.0_ti_p
         tu11_1 = 0.0_ti_p
         tu02_1 = 0.0_ti_p
         tu00_2 = 0.0_ti_p
         tu01_2 = 0.0_ti_p
         tu10_2 = 0.0_ti_p
         tu20_2 = 0.0_ti_p
         tu11_2 = 0.0_ti_p
         tu02_2 = 0.0_ti_p
         tu00 = 0.0_ti_p
         tu10 = 0.0_ti_p
         tu01 = 0.0_ti_p
         tu20 =0.0_ti_p
         tu11 =0.0_ti_p
         tu02 = 0.0_ti_p
         tu30 = 0.0_ti_p
         tu21 = 0.0_ti_p
         tu12 = 0.0_ti_p
         tu03 = 0.0_ti_p
         j0 = jgrd0
         do it2 = 1, bsorder
            j0 = j0 + 1
c           j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            j = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
            u0 = thetai2(1,it2,impi)
            u1 = thetai2(2,it2,impi)
            u2 = thetai2(3,it2,impi)
            u3 = thetai2(4,it2,impi)
            t0_1 = 0.0_ti_p
            t1_1 = 0.0_ti_p
            t2_1 = 0.0_ti_p
            t0_2 = 0.0_ti_p
            t1_2 = 0.0_ti_p
            t2_2 = 0.0_ti_p
            t3 = 0.0_ti_p
            i0 = igrd0
            do it1 = 1, bsorder
               i0 = i0 + 1
c              i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
               i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
c
               tq_1 = 0_ti_p
               tq_2 = 0_ti_p
               if (use_pmecore) then
                 kstart = kstart1(rank_bis+1)
                 kend = kend1(rank_bis+1)
                 jstart = jstart1(rank_bis+1)
                 jend = jend1(rank_bis+1)
                 istart = istart1(rank_bis+1)
                 iend = iend1(rank_bis+1)
               else
                 kstart = kstart1(rank+1)
                 kend = kend1(rank+1)
                 jstart = jstart1(rank+1)
                 jend = jend1(rank+1)
                 istart = istart1(rank+1)
                 iend = iend1(rank+1)
               end if
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 tq_1  = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $             1)
                 tq_2  = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,
     $             1)
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   tq_1  = qgrid2in_2d(1,i-istart+1,j-jstart+1,
     $               k-kstart+1,iproc+1)
                   tq_2  = qgrid2in_2d(2,i-istart+1,j-jstart+1,
     $               k-kstart+1,iproc+1)
                   goto 10
                 end if
               end do
 10            continue
c
               t0_1 = t0_1 + tq_1*thetai1(1,it1,impi)
               t1_1 = t1_1 + tq_1*thetai1(2,it1,impi)
               t2_1 = t2_1 + tq_1*thetai1(3,it1,impi)
               t0_2 = t0_2 + tq_2*thetai1(1,it1,impi)
               t1_2 = t1_2 + tq_2*thetai1(2,it1,impi)
               t2_2 = t2_2 + tq_2*thetai1(3,it1,impi)
               t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,impi)
            end do
            tu00_1 = tu00_1 + t0_1*u0
            tu10_1 = tu10_1 + t1_1*u0
            tu01_1 = tu01_1 + t0_1*u1
            tu20_1 = tu20_1 + t2_1*u0
            tu11_1 = tu11_1 + t1_1*u1
            tu02_1 = tu02_1 + t0_1*u2
            tu00_2 = tu00_2 + t0_2*u0
            tu10_2 = tu10_2 + t1_2*u0
            tu01_2 = tu01_2 + t0_2*u1
            tu20_2 = tu20_2 + t2_2*u0
            tu11_2 = tu11_2 + t1_2*u1
            tu02_2 = tu02_2 + t0_2*u2
            t0 = t0_1 + t0_2
            t1 = t1_1 + t1_2
            t2 = t2_1 + t2_2
            tu00 = tu00 + t0*u0
            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
            tu20 = tu20 + t2*u0
            tu11 = tu11 + t1*u1
            tu02 = tu02 + t0*u2
            tu30 = tu30 + t3*u0
            tu21 = tu21 + t2*u1
            tu12 = tu12 + t1*u2
            tu03 = tu03 + t0*u3
         end do
         tuv100_1 = tuv100_1 + tu10_1*v0
         tuv010_1 = tuv010_1 + tu01_1*v0
         tuv001_1 = tuv001_1 + tu00_1*v1
         tuv200_1 = tuv200_1 + tu20_1*v0
         tuv020_1 = tuv020_1 + tu02_1*v0
         tuv002_1 = tuv002_1 + tu00_1*v2
         tuv110_1 = tuv110_1 + tu11_1*v0
         tuv101_1 = tuv101_1 + tu10_1*v1
         tuv011_1 = tuv011_1 + tu01_1*v1
         tuv100_2 = tuv100_2 + tu10_2*v0
         tuv010_2 = tuv010_2 + tu01_2*v0
         tuv001_2 = tuv001_2 + tu00_2*v1
         tuv200_2 = tuv200_2 + tu20_2*v0
         tuv020_2 = tuv020_2 + tu02_2*v0
         tuv002_2 = tuv002_2 + tu00_2*v2
         tuv110_2 = tuv110_2 + tu11_2*v0
         tuv101_2 = tuv101_2 + tu10_2*v1
         tuv011_2 = tuv011_2 + tu01_2*v1
         tuv000 = tuv000 + tu00*v0
         tuv100 = tuv100 + tu10*v0
         tuv010 = tuv010 + tu01*v0
         tuv001 = tuv001 + tu00*v1
         tuv200 = tuv200 + tu20*v0
         tuv020 = tuv020 + tu02*v0
         tuv002 = tuv002 + tu00*v2
         tuv110 = tuv110 + tu11*v0
         tuv101 = tuv101 + tu10*v1
         tuv011 = tuv011 + tu01*v1
         tuv300 = tuv300 + tu30*v0
         tuv030 = tuv030 + tu03*v0
         tuv003 = tuv003 + tu00*v3
         tuv210 = tuv210 + tu21*v0
         tuv201 = tuv201 + tu20*v1
         tuv120 = tuv120 + tu12*v0
         tuv021 = tuv021 + tu02*v1
         tuv102 = tuv102 + tu10*v2
         tuv012 = tuv012 + tu01*v2
         tuv111 = tuv111 + tu11*v1
      end do
      fdip_phi1(2) = tuv100_1
      fdip_phi1(3) = tuv010_1
      fdip_phi1(4) = tuv001_1
      fdip_phi1(5) = tuv200_1
      fdip_phi1(6) = tuv020_1
      fdip_phi1(7) = tuv002_1
      fdip_phi1(8) = tuv110_1
      fdip_phi1(9) = tuv101_1
      fdip_phi1(10) = tuv011_1
      fdip_phi2(2) = tuv100_2
      fdip_phi2(3) = tuv010_2
      fdip_phi2(4) = tuv001_2
      fdip_phi2(5) = tuv200_2
      fdip_phi2(6) = tuv020_2
      fdip_phi2(7) = tuv002_2
      fdip_phi2(8) = tuv110_2
      fdip_phi2(9) = tuv101_2
      fdip_phi2(10) = tuv011_2
      fdip_sum_phi(1) = tuv000
      fdip_sum_phi(2) = tuv100
      fdip_sum_phi(3) = tuv010
      fdip_sum_phi(4) = tuv001
      fdip_sum_phi(5) = tuv200
      fdip_sum_phi(6) = tuv020
      fdip_sum_phi(7) = tuv002
      fdip_sum_phi(8) = tuv110
      fdip_sum_phi(9) = tuv101
      fdip_sum_phi(10) = tuv011
      fdip_sum_phi(11) = tuv300
      fdip_sum_phi(12) = tuv030
      fdip_sum_phi(13) = tuv003
      fdip_sum_phi(14) = tuv210
      fdip_sum_phi(15) = tuv201
      fdip_sum_phi(16) = tuv120
      fdip_sum_phi(17) = tuv021
      fdip_sum_phi(18) = tuv102
      fdip_sum_phi(19) = tuv012
      fdip_sum_phi(20) = tuv111
c
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine cmp_to_fmp  --  transformation of multipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "cmp_to_fmp" transforms the atomic multipoles from Cartesian
c     to fractional coordinates
c
c
      subroutine cmp_to_fmp(cmp,fmp)
      use mpole
      implicit none
      integer i,j,k
      real(t_p) ctf(10,10)
      real(t_p) cmp(10,*)
      real(t_p) fmp(10,*)
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac (ctf)
c
c     apply the transformation to get the fractional multipoles
c
      do i = 1, npole
         fmp(1,i) = ctf(1,1) * cmp(1,i)
         do j = 2, 4
            fmp(j,i) = 0.0
            do k = 2, 4
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
         do j = 5, 10
            fmp(j,i) = 0.0
            do k = 5, 10
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
      end do
      return
      end
c
c     "cmp_to_fmp_site" transforms the ith atomic multipole from Cartesian
c     to fractional coordinates
c
c
      subroutine cmp_to_fmp_site(cmp,fmp)
!$acc routine seq
!$acc routine(cart_to_frac) seq
      use sizes
      implicit none
      integer j,k
      real(t_p) ctf(10,10)
      real(t_p) cmp(*)
      real(t_p) fmp(*)
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac (ctf)
c
c     apply the transformation to get the fractional multipoles
c
      fmp(1) = ctf(1,1) * cmp(1)
!$acc loop seq
      do j = 2, 4
         fmp(j) = 0.0
         do k = 2, 4
            fmp(j) = fmp(j) + ctf(j,k)*cmp(k)
         end do
      end do
!$acc loop seq
      do j = 5, 10
         fmp(j) = 0.0
         do k = 5, 10
            fmp(j) = fmp(j) + ctf(j,k)*cmp(k)
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine cart_to_frac  --  Cartesian to fractional  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "cart_to_frac" computes a transformation matrix to convert
c     a multipole object in Cartesian coordinates to fractional
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine cart_to_frac (ctf)
!$acc routine
      use boxes
      use pme
      use tinheader
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real(t_p) a(3,3)
      real(t_p) ctf(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do
c
c     get the Cartesian to fractional conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ctf(j,i) = 0.0_ti_p
         end do
      end do
      ctf(1,1) = 1.0_ti_p
      do i = 2, 4
         do j = 2, 4
            ctf(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i)*a(m,j) + a(k,j)*a(m,i)
         end do
      end do
      return
      end
c
      subroutine cart_to_fracgpu
      use boxes
      use pme
      use tinheader
      use utilgpu,only:rec_queue,ctf,qi1,qi2
      implicit none
      integer i,j,k,m
      integer i1,i2
      real(t_p) a(3,3)
c
c     set the reciprocal vector transformation matrix
c

!$acc kernels async(rec_queue) present(ctf)
!$acc loop collapse(2)
      do i = 1, 10
         do j = 1, 10
            ctf(j,i) = 0.0_ti_p
         end do
      end do
      ctf(1,1) = 1.0_ti_p
!$acc end kernels

!$acc kernels async(rec_queue)
!$acc&        create(a) present(ctf,recip)
!$acc loop
      do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do
c
c     get the Cartesian to fractional conversion matrix
c
!$acc loop collapse(2)
      do i = 2, 4
         do j = 2, 4
            ctf(i,j) = a(i-1,j-1)
         end do
      end do
!$acc loop collapse(2)
      do i2 = 1, 6
         do i1 = 1, 3
            k = qi1(i1)
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i) * a(k,j)
         end do
      end do
!$acc loop seq collapse(2)
      do i2 = 1, 6
         do i1 = 4, 6
            k = qi1(i1)
            m = qi2(i1)
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i)*a(m,j) + a(k,j)*a(m,i)
         end do
      end do
!$acc end kernels
      end
c
      subroutine fphi_to_cphi_site(fphi,cphi)
!$acc routine
      use sizes
      use tinheader
      implicit none
      integer j,k
      real(t_p) ftc(10,10)
      real(t_p) cphi(10)
      real(t_p) fphi(20)
!$acc routine(frac_to_cart) seq
c
c     find the matrix to convert fractional to Cartesian
c
      call frac_to_cart (ftc)
c
c     apply the transformation to get the Cartesian potential
c
      cphi(1) = ftc(1,1) * fphi(1)
      do j = 2, 4
         cphi(j) = 0.0_ti_p
         do k = 2, 4
            cphi(j) = cphi(j) + ftc(j,k)*fphi(k)
         end do
      end do
      do j = 5, 10
         cphi(j) = 0.0_ti_p
         do k = 5, 10
            cphi(j) = cphi(j) + ftc(j,k)*fphi(k)
         end do
      end do
      return
      end
c
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine frac_to_cart  --  fractional to Cartesian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "frac_to_cart" computes a transformation matrix to convert
c     a multipole object in fraction coordinates to Cartesian
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine frac_to_cart (ftc)
!$acc routine seq
c      use sizes
      use boxes
      use pme
      use tinheader
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real(t_p) a(3,3)
      real(t_p) ftc(10,10)
      parameter (
     &       qi1=(/ 1, 2, 3, 1, 1, 2 /),
     &       qi2=(/ 1, 2, 3, 2, 3, 3 /))
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(i,1) = real(nfft1,t_p) * recip(i,1)
         a(i,2) = real(nfft2,t_p) * recip(i,2)
         a(i,3) = real(nfft3,t_p) * recip(i,3)
      end do
c
c     get the fractional to Cartesian conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ftc(j,i) = 0.0_ti_p
         end do
      end do
      ftc(1,1) = 1.0_ti_p
      do i = 2, 4
         do j = 2, 4
            ftc(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(k,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = 2.0_ti_p * a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(m,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = a(k,i)*a(m,j) + a(m,i)*a(k,j)
         end do
      end do
      return
      end
c
      subroutine frac_to_cartgpu
c     use sizes
      use boxes
      use pme
      use tinheader ,only: ti_p
      use utilgpu,only:rec_queue,ftc,qi1,qi2
      implicit none
      integer i,j,k,m
      integer i1,i2
      real(t_p) a(3,3)
c
c     set the reciprocal vector transformation matrix
c
!$acc data create(a)
!$acc&     present(recip,ftc) 
!$acc&     async(rec_queue)

!$acc kernels async(rec_queue)
!$acc loop collapse(2)
      do i = 1, 10
         do j = 1, 10
            ftc(j,i) = 0.0_ti_p
         end do
      end do
      ftc(1,1) = 1.0_ti_p
!$acc end kernels
!$acc kernels async(rec_queue)
!$acc loop
      do i = 1, 3
         a(i,1) = real(nfft1,t_p) * recip(i,1)
         a(i,2) = real(nfft2,t_p) * recip(i,2)
         a(i,3) = real(nfft3,t_p) * recip(i,3)
      end do
c
c     get the fractional to Cartesian conversion matrix
c
!$acc loop collapse(2)
      do i = 2, 4
         do j = 2, 4
            ftc(i,j) = a(i-1,j-1)
         end do
      end do
!$acc loop collapse(2)
      do i2 = 1, 6
         do i1 = 1, 3
            k = qi1(i1)
            i = qi1(i2)
            if (i2.lt.4) then
            ftc(i1+4,i2+4) = a(k,i)*a(k,i)
            else
            j = qi2(i2)
            ftc(i1+4,i2+4) = 2.0_ti_p * a(k,i)*a(k,j)
            end if
         end do
      end do
!$acc loop collapse(2)
      do i2 = 1, 6
         do i1 = 4, 6
            k = qi1(i1)
            m = qi2(i1)
            i = qi1(i2)
            if (i2.lt.4) then
               ftc(i1+4,i2+4) = a(k,i)*a(m,i)
            else
               j = qi2(i2)
               ftc(i1+4,i2+4) = a(k,i)*a(m,j) + a(m,i)*a(k,j)
            end if
         end do
      end do
!$acc end kernels

!$acc end data
      return
      end
c
c     subroutine reassignpme: initialize arrays to compute reciprocal part of 
c     electrostatic and polarization interactions, also reassign atoms on 'reciprocal cores'
c     between time steps
c
      subroutine reassignpme(init)
      use atomsMirror
      use boxes
      use chunks
      use domdec
      use fft
      use iso_c_binding
      use inform    ,only: deb_Path
      use math
      use pme
      use pme1
      use potent
      use mpi
      use timestat  ,only: timer_enter,timer_exit,timer_eneig
     &              ,quiet_timers
      use tinheader ,only:ti_p,re_p
      use tinMemory
      use utilcomm  ,only: repartrectemp=>buffMpi_i1
      use utils     ,only: set_to_zero1
      use utilgpu   ,only: rec_queue
      implicit none
      real(t_p) disttemp
      real(t_p) jtemp(4),ktemp(4),jtempbis(4),ktempbis(4)
      integer ierr,iglob,j
      integer i,k
      integer iproc,iproc1,temp
      integer proc1,proc2
      integer iprec,idomlen,ibeg
      integer nprocloc,rankloc,nprocloc1
      integer gridin_size,gridout_size
      integer, allocatable :: req(:)
      integer, allocatable :: counttemp(:)
      integer, allocatable :: domlentemp(:)
      integer, allocatable :: nrecdirrecep(:),precdirrecep(:,:)
      integer, allocatable :: nrecdirreceptemp(:)
      integer, allocatable :: nrecrecep(:),precrecep(:,:)
      integer ifr,kbeg,kend,kmin,kmax,sizek,sizej
      integer jbeg,jend,jmin,jmax
      real(t_p) wbeg,wend
      integer tag, status(MPI_STATUS_SIZE)
      integer nlocrecold
      real(t_p) fr,xi,yi,zi,eps,w
      real(t_p) dist,xr,yr
      type(c_ptr) grid_ptr
      logical init
c
#if defined(SINGLE) || defined(MIXED)
      eps = 1.0d-4
#else
      eps = 1.0d-8
#endif
      allocate (req(2*nproc))
      allocate (counttemp(nproc))
      call prmem_request(repartrectemp,n)
      allocate (domlentemp(nproc))
      counttemp = 0
c
c     get geometrical correspondance between processes
c
      if (init) then
        if (deb_Path) print*, 'reassignpme init'
        call timer_enter( timer_eneig )
        if (use_pmecore) then
          nprocloc = ndir
          rankloc  = rank_bis
        else
          nprocloc = nproc
          rankloc  = rank
        end if
        call set_pme_eps
        allocate (nrecdirrecep(nproc))
        allocate (precdirrecep(nproc,nproc))
        allocate (nrecdirreceptemp(nproc))
        allocate (nrecrecep(nproc))
        allocate (precrecep(nproc,nproc))
        nrecdirrecep = 0
        nrecdirreceptemp = 0
        precdirrecep = 0
        nrecrecep = 0
        precrecep = 0
c
        if (.not.allocated(prec_send))
     $    allocate (prec_send(nproc))
        if (.not.allocated(prec_recep))
     $    allocate (prec_recep(nproc))
        !Remove declare directive before
        !call prmem_request(prec_send,nproc)
        !call prmem_request(prec_recep,nproc)
        prec_recep = 0
        prec_send  = 0
        call prmem_request(repartrec,n)
        call prmem_request(domlenrec,nproc)
c
        precdir_recep  = 0
        precdir_send   = 0
        precdir_recep1 = 0
        precdir_send1  = 0
c
        if (use_pmecore) then
          nprocloc1 = nrec
        else
          nprocloc1 = nproc
        end if
        do iproc = 0, nprocloc1 - 1
          kmin = kstart1(iproc+1)
          kmax = kend1(iproc+1)
          jmin = jstart1(iproc+1)
          jmax = jend1(iproc+1)
          if (kmin.lt.1) kmin = kmin + nfft3
          if (kmax.gt.nfft3) kmax = kmax - nfft3
          if (jmin.lt.1) jmin = jmin + nfft2
          if (jmax.gt.nfft2) jmax = jmax - nfft2
c
          do iproc1 = 0, nprocloc - 1
c
            wbeg = (zbegproc(iproc1+1))*recip(3,3)
            fr = real(nfft3,t_p) * (wbeg+0.5_ti_p)
            ifr = int(fr-pme_eps)
            !if (fr-ifr.eq.0.0.and.ifr>0) ifr = ifr-1
            kbeg = ifr - bsorder + grdoff
            wend = (zendproc(iproc1+1))*recip(3,3)
            fr = real(nfft3,t_p) * (wend+0.5_ti_p)
            ifr = int(fr-pme_eps)
            !if (fr-ifr.eq.0.0.and.ifr>0) ifr = ifr-1
            kend = ifr - bsorder + grdoff
            sizek = abs(kend-kbeg)
            if (kbeg.lt.1) kbeg = kbeg + nfft3
            if (kbeg+sizek-1.le.nfft3) then
              kend = kbeg + sizek + 1
            else
              temp = nfft3 - kbeg + 1
              kend = sizek - temp + 1
            end if
c
            wbeg = (ybegproc(iproc1+1))*recip(2,2)
            fr = real(nfft2,t_p) * (wbeg+0.5_ti_p)
            ifr = int(fr-pme_eps)
            !if (fr-ifr.eq.0.0.and.ifr>0) ifr = ifr-1
            jbeg = ifr - bsorder + grdoff
            wend = (yendproc(iproc1+1))*recip(2,2)
            fr = real(nfft2,t_p) * (wend+0.5_ti_p)
            ifr = int(fr-pme_eps)
            !if (fr-ifr.eq.0.0.and.ifr>0) ifr = ifr-1
            jend = ifr - bsorder + grdoff
            sizej = abs(jend-jbeg)
            if (jbeg.lt.1) jbeg = jbeg + nfft2
            if (jbeg+sizej-1.le.nfft2) then
              jend = jbeg + sizej + 1
            else
              temp = nfft2 - jbeg + 1
              jend = sizej - temp + 1
            end if
c
            if ((kbeg.le.kend).and.(jbeg.le.jend)) then
              if ((((kbeg.le.kmin).and.(kend.ge.kmin)).or.
     $          ((kbeg.ge.kmin).and.(kbeg.le.kmax))).and.
     $          (((jbeg.le.jmin).and.(jend.ge.jmin)).or.
     $          ((jbeg.ge.jmin).and.(jbeg.le.jmax))))
     $         then
                nrecdirrecep(iproc+1) = nrecdirrecep(iproc+1)+1
                precdirrecep(nrecdirrecep(iproc+1),iproc+1) = iproc1
              end if
            else if ((kbeg.le.kend).and.(jbeg.gt.jend)) then
              if ((((kbeg.le.kmin).and.(kend.ge.kmin)).or.
     $          ((kbeg.ge.kmin).and.(kbeg.le.kmax))).and.
     $          (((jbeg.le.jmax)).or.((jend.ge.jmin))))
     $         then
                nrecdirrecep(iproc+1) = nrecdirrecep(iproc+1)+1
                precdirrecep(nrecdirrecep(iproc+1),iproc+1) = iproc1
              end if
            else if ((kbeg.gt.kend).and.(jbeg.le.jend)) then
              if ((((jbeg.le.jmin).and.(jend.ge.jmin)).or.
     $          ((jbeg.ge.jmin).and.(jbeg.le.jmax))).and.
     $          (((kbeg.le.kmax)).or.((kend.ge.kmin))))
     $         then
                nrecdirrecep(iproc+1) = nrecdirrecep(iproc+1)+1
                precdirrecep(nrecdirrecep(iproc+1),iproc+1) = iproc1
              end if
            else
              if (((jbeg.le.jmax).or.(jend.ge.jmin)).and.
     $          ((kbeg.le.kmax).or.(kend.ge.kmin))) then
                nrecdirrecep(iproc+1) = nrecdirrecep(iproc+1)+1
                precdirrecep(nrecdirrecep(iproc+1),iproc+1) = iproc1
              end if
            end if
          end do
        end do
        if (use_pmecore) then
          if (rank.gt.ndir-1) then
            nrecdir_recep = nrecdirrecep(rankloc+1)
            precdir_recep = precdirrecep(:,rankloc+1)
          end if
        else
          nrecdir_recep = nrecdirrecep(rankloc+1)
          precdir_recep = precdirrecep(:,rankloc+1)
        end if
 345    format(A5,10I5)
c
c   Receive also the neighbors within a cutoff to rotate reciprocal multipoles
c
        if (use_pmecore) then
          nprocloc = ndir
          rankloc  = rank
        else
          nprocloc = nproc
          rankloc  = rank
        end if
c
        nrecdirreceptemp = nrecdirrecep
c
c      store this list (without the processes to rotate the reciprocal multipoles)
c
        nrecdir_recep1 = nrecdir_recep
        precdir_recep1 = precdir_recep
c
c     get the processes to send data to
c
        nrecdir_send1 = 0
c
        if (use_pmecore) then
          if (rank.le.ndir-1) then
            do iproc = 0, nrec-1
              do i = 1, nrecdirrecep(iproc+1)
                if (precdirrecep(i,iproc+1).eq.rankloc) then
                  nrecdir_send1 = nrecdir_send1 + 1
                  precdir_send1(nrecdir_send1) = iproc+ndir
                end if
              end do
            end do
          end if
        else
          do iproc = 0, nprocloc-1
            do i = 1, nrecdirrecep(iproc+1)
              if (precdirrecep(i,iproc+1).eq.rankloc) then
                nrecdir_send1 = nrecdir_send1 + 1
                precdir_send1(nrecdir_send1) = iproc
              end if
            end do
          end do
        end if
        if (use_pmecore) then
          if (rank.gt.ndir-1) then
            nrecdir_recep = nrecdirreceptemp(rankloc+1-ndir)
            nrecdirrecep  = nrecdirreceptemp
            precdir_recep = precdirrecep(:,rankloc+1-ndir)
          end if
        else
          nrecdir_recep = nrecdirreceptemp(rankloc+1)
          precdir_recep = precdirrecep(:,rankloc+1)
        end if
c
c     get the processes to send data to
c
        nrecdir_send = 0
c
        if (use_pmecore) then
          if (rank.le.ndir-1) then
            do iproc = 0, nrec-1
              do i = 1, nrecdirreceptemp(iproc+1)
                if (precdirrecep(i,iproc+1).eq.rankloc) then
                  nrecdir_send = nrecdir_send + 1
                  precdir_send(nrecdir_send) = iproc+ndir
                end if
              end do
            end do
          end if
        else
          do iproc = 0, nprocloc-1
            do i = 1, nrecdirreceptemp(iproc+1)
              if (precdirrecep(i,iproc+1).eq.rankloc) then
                nrecdir_send = nrecdir_send + 1
                precdir_send(nrecdir_send) = iproc
              end if
            end do
          end do
        end if
c
c      Remove procs already in real space neighbor list and put them in another list
c
        nrecdir_recep2 = 0
        precdir_recep2 = 0
        do i = 1, nrecdir_recep
          proc1 = precdir_recep(i)
          do j = 1, nbig_recep
            proc2 = pbig_recep(j)
            if (proc1.eq.proc2) goto 20
          end do
          nrecdir_recep2 = nrecdir_recep2 + 1
          precdir_recep2(nrecdir_recep2) = proc1
 20     continue
        end do
cc
        nrecdir_send2 = 0
        precdir_send2 = 0
        do i = 1, nrecdir_send
          proc1 = precdir_send(i)
          do j = 1, nbig_send
            proc2 = pbig_send(j)
            if (proc1.eq.proc2) goto 30
          end do
          nrecdir_send2 = nrecdir_send2 + 1
          precdir_send2(nrecdir_send2) = proc1
 30     continue
        end do
c
c      Remove procs already in recdir space neighbor list and put them in another list
c
c       nrecdir_recep3 = 0
c       precdir_recep3 = 0
c       do i = 1, nrecdir_recep
c         proc1 = precdir_recep(i)
c         do j = 1, nrecdir_recep2
c           proc2 = precdir_recep2(j)
c           if (proc1.eq.proc2) goto 40
c         end do
c         nrecdir_recep3 = nrecdir_recep3 + 1
c         precdir_recep3(nrecdir_recep3) = proc1
c40     continue
c       end do
cc
c       nrecdir_send3 = 0
c       precdir_send3 = 0
c       do i = 1, nrecdir_send
c         proc1 = precdir_send(i)
c         do j = 1, nrecdir_send3
c           proc2 = precdir_send3(j)
c           if (proc1.eq.proc2) goto 50
c         end do
c         nrecdir_send3 = nrecdir_send3 + 1
c         precdir_send3(nrecdir_send3) = proc1
c50     continue
c       end do
c
        if (use_pmecore) then
          nprocloc = nrec
          rankloc  = rank_bis
        else
          nprocloc = nproc
          rankloc  = rank
        end if
c
        if ((.not.(use_pmecore)).or.(rank.gt.ndir-1)) then
c
c     get the reciprocal processes to send and receive data from
c
c
          nrec_recep = 0
          do iproc = 0, nprocloc-1
            do iproc1 = 0, nprocloc-1
              dist = 1000.0_ti_p
              if (iproc.eq.iproc1) goto 100
              jtemp(1) = jstart1(iproc+1)
              ktemp(1) = kstart1(iproc+1)
              jtemp(2) = jend1(iproc+1)
              ktemp(2) = kstart1(iproc+1)
              jtemp(3) = jend1(iproc+1)
              ktemp(3) = kend1(iproc+1)
              jtemp(4) = jstart1(iproc+1)
              ktemp(4) = kend1(iproc+1)
              jtempbis(1) = jstart1(iproc1+1)
              ktempbis(1) = kstart1(iproc1+1)
              jtempbis(2) = jend1(iproc1+1)
              ktempbis(2) = kstart1(iproc1+1)
              jtempbis(3) = jend1(iproc1+1)
              ktempbis(3) = kend1(iproc1+1)
              jtempbis(4) = jstart1(iproc1+1)
              ktempbis(4) = kend1(iproc1+1)
              do i = 1, 4
                do j = 1, 4
                   xr = abs(jtemp(i) - jtempbis(j))
                   yr = abs(ktemp(i) - ktempbis(j))
                   if (xr.ge.nfft2-2) xr = nfft2-xr
                   if (yr.ge.nfft3-2) yr = nfft3-yr
                   disttemp = sqrt(xr*xr+yr*yr)
                   if (disttemp.le.dist) dist = disttemp
                 end do
               end do
              if ((iproc1.ne.iproc).and.(dist.le.1.5*max(nlpts,nrpts)))
     $          then
                nrecrecep(iproc+1) = nrecrecep(iproc+1)+1
                precrecep(nrecrecep(iproc+1),iproc+1) = iproc1
              end if
 100          continue
            end do
          end do
          nrec_recep = nrecrecep(rankloc+1)
          prec_recep = precrecep(:,rankloc+1)
c
c     get the processes to send data to
c
          nrec_send = 0
c
          do iproc = 0, nprocloc-1
            do i = 1, nrecrecep(iproc+1)
              if (precrecep(i,iproc+1).eq.rankloc) then
                nrec_send = nrec_send + 1
                prec_send(nrec_send) = iproc
              end if
            end do
          end do

          if ( rank.eq.0.and.nproc.gt.1.and.tinkerdebug.gt.0 ) then
  46         format(A20,I3,';',1x,$)
  47         format(I3,$)
              write(*,46)'nrec_send      = ',nrec_send
              write(*,47) (prec_send(i),i=1,nrec_send)
              write(*,*)
              write(*,46) 'nrec_recep     = ',nrec_recep
              write(*,47) (prec_recep(i),i=1,nrec_recep)
              write(*,*)
              write(*,46) 'nrecdir_send   = ',nrecdir_send
              write(*,47) (precdir_send(i),i=1,nrecdir_send)
              write(*,*)
              write(*,46) 'nrecdir_send1  = ',nrecdir_send1
              write(*,47) (precdir_send1(i),i=1,nrecdir_send1)
              write(*,*)
              write(*,46) 'nrecdir_send2  = ',nrecdir_send2
              write(*,47) (precdir_send2(i),i=1,nrecdir_send2)
              write(*,*)
              write(*,46) 'nrecdir_recep  = ',nrecdir_recep
              write(*,47) (precdir_recep(i),i=1,nrecdir_recep)
              write(*,*)
              write(*,46) 'nrecdir_recep1 = ',nrecdir_recep1
              write(*,47) (precdir_recep1(i),i=1,nrecdir_recep1)
              write(*,*)
              write(*,46) 'nrecdir_recep2 = ',nrecdir_recep2
              write(*,47) (precdir_recep2(i),i=1,nrecdir_recep2)
              write(*,*)
          end if
!$acc update device(prec_send,prec_recep) async
!$acc update device(nrec_send) async
c
          gridin_size  = 2*(nrec_send+1)*n1mpimax*n2mpimax*n3mpimax
          gridout_size = 2*isize2(rankloc+1)*jsize2(rankloc+1)*
     $                     ksize2(rankloc+1)
c
          if (allocated(qgridin_2d)) then  !Track memory consumption
          s_prmem = s_prmem -2*sizeof(qgridin_2d) -2*sizeof(qgridout_2d)
#ifdef _OPENACC
         sd_prmem =sd_prmem -2*sizeof(qgridin_2d) -2*sizeof(qgridout_2d)
!$acc exit data delete(qgridin_2d,qgrid2in_2d,qgridout_2d,qgrid2out_2d)
#endif
          end if
          ! Removed declare directive before using tinker allocator
          if (allocated(qgridin_2d)) deallocate (qgridin_2d)
          allocate (qgridin_2d(2,n1mpimax,n2mpimax,n3mpimax,
     $       nrec_send+1))
          if (allocated(qgrid2in_2d)) deallocate (qgrid2in_2d)
          allocate (qgrid2in_2d(2,n1mpimax,n2mpimax,n3mpimax,
     $       nrec_send+1))
          if (allocated(qgridout_2d)) deallocate (qgridout_2d)
          allocate (qgridout_2d(2,isize2(rankloc+1),jsize2(rankloc+1),
     $       ksize2(rankloc+1)))
          if (allocated(qgrid2out_2d)) deallocate (qgrid2out_2d)
          allocate (qgrid2out_2d(2,isize2(rankloc+1),jsize2(rankloc+1),
     $       ksize2(rankloc+1)))
          s_prmem = s_prmem +2*sizeof(qgridin_2d) +2*sizeof(qgridout_2d)
#ifdef _OPENACC
         sd_prmem =sd_prmem +2*sizeof(qgridin_2d) +2*sizeof(qgridout_2d)
!$acc enter data create(qgridin_2d,qgrid2in_2d,qgridout_2d,qgrid2out_2d)
#endif
c
c     Associate 1D grid pointer with his 5D grid target
c
          grid_ptr = c_loc(qgridin_2d)
          call c_f_pointer(grid_ptr,qgridin_p,(/gridin_size/))
          grid_ptr = c_loc(qgrid2in_2d)
          call c_f_pointer(grid_ptr,qgrid2in_p,(/gridin_size/))
          grid_ptr = c_loc(qgridout_2d)
          call c_f_pointer(grid_ptr,qgridout_p,(/gridout_size/))
          grid_ptr = c_loc(qgrid2out_2d)
          call c_f_pointer(grid_ptr,qgrid2out_p,(/gridout_size/))
c
!$acc enter data attach(qgridin_p,qgrid2in_p,
!$acc&      qgridout_p,qgrid2out_p)
!!$acc enter data create(qgridin_2d,qgrid2in_2d,
!!$acc&    qgridout_2d,qgrid2out_2d)

#ifdef _OPENACC
          ! Attach CUDA grid device pointers
          call attach_pmecu_pointer(0)
#endif

          if (nrec_recep.gt.0) then
             ! Allocate mpi grid
             if (allocated(qgridmpi)) then
!$acc exit data delete(qgridmpi)
                deallocate(qgridmpi)
                s_prmem = s_prmem - size(qgridmpi)*szoTp
#ifdef _OPENACC
               sd_prmem =sd_prmem - size(qgridmpi)*szoTp
#endif
             end if
             allocate(qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
!$acc enter data create(qgridmpi)
             s_prmem = s_prmem + size(qgridmpi)*szoTp
#ifdef _OPENACC
            sd_prmem =sd_prmem + size(qgridmpi)*szoTp
#endif
            call set_to_zero1(qgridmpi,size(qgridmpi),rec_queue)
c
             if (GridDecomp1d) then
                call set_MPI_GRID1D_type(bsorder/2)
             end if
          end if
        end if
c
        domlenrec = 0
        do i = 1, n
          xi = x(i)
          yi = y(i)
          zi = z(i)
          w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
          fr = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
          ifr = int(fr-pme_eps)
          !if (ifr-fr.eq.0.0.and.ifr>0) ifr=ifr-1
          k = ifr- bsorder + grdoff
          if (k .lt. 1)  k = k + nfft3
          w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
          fr = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
          ifr = int(fr-pme_eps)
          !if (ifr-fr.eq.0.0.and.ifr>0) ifr=ifr-1
          j = ifr- bsorder + grdoff
          if (j .lt. 1)  j = j + nfft2
          do iproc = 0, nprocloc-1
            if (((k.ge.kstart1(iproc+1)).and.(k.le.kend1(iproc+1))).and.
     $          ((j.ge.jstart1(iproc+1)).and.(j.le.jend1(iproc+1))))then
              repartrec(i) = iproc
            end if
          end do
        end do
!$acc update device(repartrec)

        domlentemp = domlen
        if(nproc.gt.1) call orderbuffer_gpu(.false.)
        call orderbufferrec
c
        deallocate (precdirrecep)
        deallocate (nrecdirrecep)
        deallocate (nrecdirreceptemp)
        deallocate (precrecep)
        deallocate (nrecrecep)
        call timer_exit( timer_eneig,quiet_timers )
      end if
c
c
      if (nproc.eq.1) return

      call timer_enter( timer_eneig )
      if (deb_Path.and..not.init) print*,'reassignpme'

      if (use_pmecore) then
        nprocloc = nrec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        rankloc  = rank
      end if
c
      nlocrecold = nlocrec
!$acc data present(glob,x,y,z,repartrec,repartrectemp
!$acc&            ,recip,kstart1,kend1,jstart1,jend1)

!$acc parallel loop async
      do i = 1, nloc
        iglob = glob(i)
        xi    = x(iglob)
        yi    = y(iglob)
        zi    = z(iglob)
        w     = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
        fr    = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
        ifr   = int(fr-pme_eps)
        !if (ifr-fr.eq.0.0.and.ifr>0) ifr=ifr-1
        k     = ifr- bsorder + grdoff
        if (k .lt. 1)  k = k + nfft3
        w     = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
        fr    = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
        ifr   = int(fr-pme_eps)
        !if (ifr-fr.eq.0.0) ifr=ifr-1
        j     = ifr- bsorder + grdoff
        if (j .lt. 1)  j = j + nfft2
!$acc loop seq
        do iproc = 0, nprocloc-1
          if (((k.ge.kstart1(iproc+1)).and.(k.le.kend1(iproc+1))).and.
     $        ((j.ge.jstart1(iproc+1)).and.(j.le.jend1(iproc+1)))) then
             repartrectemp(i) = iproc
          end if
        end do
      end do
c
c     send and receive repartrectemp
c
!$acc host_data use_device(repartrectemp)
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag     = nproc*rank + precdir_recep(i) + 1
          iprec   = precdir_recep1(i)
          idomlen = domlen(iprec+1)
          call MPI_IRECV(repartrectemp(bufbeg(iprec+1)),idomlen,MPI_INT
     $        ,iprec,tag,COMM_TINKER,req(i),ierr)
        end if
      end do
c
!$acc wait
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
          tag = nproc*precdir_send1(i) + rank + 1
          call MPI_ISEND(repartrectemp,nloc,MPI_INT
     $        ,precdir_send1(i),tag,COMM_TINKER,req(nproc+i),ierr)
        end if
      end do
!$acc end host_data
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank)
     &    call MPI_WAIT(req(i),status,ierr)
      end do
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank)
     &    call MPI_WAIT(req(nproc+i),status,ierr)
      end do

!$acc parallel loop async
      do i = 1, nloc
         iglob = glob(i)
         repartrec(iglob) =  repartrectemp(i)
      end do

      do iproc = 1, nrecdir_recep1
        if (precdir_recep1(iproc).ne.rank) then
           iprec   = precdir_recep1(iproc)
           idomlen = domlen(iprec+1)
           ibeg    = bufbeg(iprec+1)-1
!$acc parallel loop async
           do i = 1, idomlen
             iglob            = glob(ibeg+i)
             repartrec(iglob) = repartrectemp(ibeg+i)
           end do
        end if
      end do

      call orderbufferrec_gpu
c
!$acc end data

      deallocate (domlentemp)
      deallocate (counttemp)
      deallocate (req)
      call timer_exit( timer_eneig,quiet_timers )
      end
c
c     Set MPI type for fftGrid
c
      subroutine set_MPI_GRID1D_type(nfloor)
      use domdec,only: rank
      use fft
      use mpi
      use pme1 , pm_nfloor=>nfloor
      implicit none
      integer,intent(in)::nfloor
      integer i,ierr
      integer ncount,block,gridin_size

      if (Is_MPI_GRID1D) then
         call MPI_Type_free(MPI_GRID1D,ierr)
         Is_MPI_GRID1D=.false.
      end if
      pm_nfloor   = nfloor
      MPI_GRID1D  = 0
      gridin_size = 2*n1mpimax*n2mpimax*n3mpimax

      ! Check limits of 1d subdivision
      if (2*nfloor.ge.n3mpimax) then
         ncount = 1
         block  = 2*n1mpimax*n2mpimax*n3mpimax 
         stride = 1
         call MPI_type_contiguous(block,MPI_TPREC,MPI_GRID1D,ierr)
      else
         ncount = 2
         block  = 2*n1mpimax*n2mpimax*nfloor
         stride = 2*n1mpimax*n2mpimax*(n3mpimax-nfloor)
         ! Set MPI type for 1D fftGrid
         call MPI_Type_vector(ncount,block,stride,MPI_TPREC
     &                       ,MPI_GRID1D,ierr)
         if (ierr.ne.MPI_SUCCESS)
     &      write(0,*) 'MPI_Type_vector Issue !!!',ierr,rank
      end if

      call MPI_Type_commit(MPI_GRID1D,ierr)
      Is_MPI_GRID1D=.true.
      if (ierr.ne.MPI_SUCCESS)
     &   write(0,*) 'MPI_Type_commit Issue !!!',ierr,rank

c 43  format( "I rank",I2," nc(",I2,") bl(",I7,") stri("
c    &      ,I8,") GrSize (",I9,") MPI_TYPE:",I3)
c     write(*,43) rank,ncount,block,stride,gridin_size
c    &           ,MPI_GRID1D
      end subroutine
c
c     "grid_pchg_site" places the i-th fractional atomic charge onto
c     the particle mesh Ewald grid
c
      subroutine grid_pchg_site(isite,impi,q)
      use chunks
      use domdec
      use fft
      use pme
      use potent
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,m,impi,rankloc
      integer iproc,proc
      integer ii,jj,kk
      integer isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer nearpt(3)
      integer abound(6)
      real(t_p) q
      real(t_p) v0,u0,t0
      real(t_p) term,uterm
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
c
c     put the permanent multipole moments onto the grid
c
      iatm = isite
      nearpt(1) = igrid(1,iatm) + grdoff
      nearpt(2) = igrid(2,iatm) + grdoff
      nearpt(3) = igrid(3,iatm) + grdoff
      abound(1) = nearpt(1) - nlpts
      abound(2) = nearpt(1) + nrpts
      abound(3) = nearpt(2) - nlpts
      abound(4) = nearpt(2) + nrpts
      abound(5) = nearpt(3) - nlpts
      abound(6) = nearpt(3) + nrpts
      call adjust (offsetx,nfft1,1,abound(1),
     &               abound(2),1,nfft1)
      call adjust (offsety,nfft2,1,abound(3),
     &               abound(4),1,nfft2)
      call adjust (offsetz,nfft3,1,abound(5),
     &               abound(6),1,nfft3)
      do kk = abound(5), abound(6)
         k = kk
         m = k + offsetz
         if (k .lt. 1)  k = k + nfft3
         v0 = thetai3(1,m,impi)
         do jj = abound(3), abound(4)
            j = jj
            m = j + offsety
            if (j .lt. 1)  j = j + nfft2
            u0 = thetai2(1,m,impi)
            term = q*u0*v0
            do ii = abound(1), abound(2)
               i = ii
               m = i + offsetx
               if (i .lt. 1)  i = i + nfft1
               t0 = thetai1(1,m,impi)
               uterm = term*t0
               
c
               kstart = kstart1(rankloc+1)
               kend = kend1(rankloc+1)
               jstart = jstart1(rankloc+1)
               jend = jend1(rankloc+1)
               istart = istart1(rankloc+1)
               iend = iend1(rankloc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
!$acc atomic update
                 qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1) =
     $           qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1) +
     $              uterm
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   uterm = term*t0
!$acc atomic update
                   qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $               iproc+1) = qgridin_2d(1,i-istart+1,j-jstart+1,
     $               k-kstart+1,iproc+1) + uterm
                   goto 10
                 end if
               end do
 10            continue
            end do
         end do
      end do
      return
      end
