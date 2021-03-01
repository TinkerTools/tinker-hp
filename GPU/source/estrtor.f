c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine estrtor  --  stretch-torsion cross term energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrtor" calculates the stretch-torsion potential energy
c
c
#include "tinker_precision.h"
      module estrtor_inl
        contains
#include "image.f.inc"
      end module

      subroutine estrtor
      use atmlst
      use atomsMirror
      use bond
      use bound
      use energi
      use estrtor_inl
      use group
      use nvshmem
      use strtor
      use tinheader
      use torpot
      use tors
      use usage
      implicit none
      integer i,k,istrtor,iistrtor
      integer ia,ib,ic,id
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,dr
      real(t_p) rt2,ru2,rtru
      real(t_p) rba,rcb,rdc
      real(t_p) e1,e2,e3
      real(t_p) xt,yt,zt
      real(t_p) xu,yu,zu
      real(t_p) xtu,ytu,ztu
      real(t_p) v1,v2,v3
      real(t_p) c1,c2,c3
      real(t_p) s1,s2,s3
      real(t_p) sine,cosine
      real(t_p) sine2,cosine2
      real(t_p) sine3,cosine3
      real(t_p) phi1,phi2,phi3
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      real(t_p),parameter::rtiny=0.000001_ti_p
      logical proceed
c
c
c     zero out the stretch-torsion energy
c
      ebt = 0.0_ti_p
c
c     calculate the stretch-torsion interaction energy term
c
!$acc parallel loop default(present) present(ebt) async
      do istrtor = 1, nstrtorloc
         iistrtor = strtorglob(istrtor)
         i    = ist(1,iistrtor)
#ifdef USE_NVSHMEM_CUDA
         ipe  = (i-1)/ntors_pe
         ind  = mod((i-1),ntors_pe) +1
         ia   = d_itors(ipe)%pel(1,ind)
         ib   = d_itors(ipe)%pel(2,ind)
         ic   = d_itors(ipe)%pel(3,ind)
         id   = d_itors(ipe)%pel(4,ind)
#else
         ia   = itors(1,i)
         ib   = itors(2,i)
         ic   = itors(3,i)
         id   = itors(4,i)
#endif
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                            use(ic) .or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            if (use_polymer) then
               call image_inl (xba,yba,zba)
               call image_inl (xcb,ycb,zcb)
               call image_inl (xdc,ydc,zdc)
            end if
            rba = sqrt(xba*xba + yba*yba + zba*zba)
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
            if (min(rba,rcb,rdc) .ne. 0.0_ti_p) then
               xt     = yba*zcb - ycb*zba
               yt     = zba*xcb - zcb*xba
               zt     = xba*ycb - xcb*yba
               xu     = ycb*zdc - ydc*zcb
               yu     = zcb*xdc - zdc*xcb
               zu     = xcb*ydc - xdc*ycb
               xtu    = yt*zu - yu*zt
               ytu    = zt*xu - zu*xt
               ztu    = xt*yu - xu*yt
               rt2    = max( xt*xt + yt*yt + zt*zt,rtiny )
               ru2    = max( xu*xu + yu*yu + zu*zu,rtiny )
               rtru   = sqrt(rt2 * ru2)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine   = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
#ifdef USE_NVSHMEM_CUDA
               c1 = d_tors1(ipe)%pel(3,ind)
               s1 = d_tors1(ipe)%pel(4,ind)
               c2 = d_tors2(ipe)%pel(3,ind)
               s2 = d_tors2(ipe)%pel(4,ind)
               c3 = d_tors3(ipe)%pel(3,ind)
               s3 = d_tors3(ipe)%pel(4,ind)
#else
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
#endif
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2   = 2.0_ti_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3   = cosine*sine2 + sine*cosine2
               phi1    = 1.0_ti_p + (cosine*c1 + sine*s1)
               phi2    = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               phi3    = 1.0_ti_p + (cosine3*c3 + sine3*s3)
c
c     get the stretch-torsion values for the first bond
c
               v1 = kst(1,iistrtor)
               v2 = kst(2,iistrtor)
               v3 = kst(3,iistrtor)
               k  = ist(2,iistrtor)
#ifdef USE_NVSHMEM_CUDA
               ipe = (k-1)/nbond_pe
               ind = mod((k-1),nbond_pe) +1
               dr  = rba - d_bl(ipe)%pel(ind)
#else
               dr = rba - bl(k)
#endif
               e1 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the stretch-torsion values for the second bond
c
               v1 = kst(4,iistrtor)
               v2 = kst(5,iistrtor)
               v3 = kst(6,iistrtor)
               k  = ist(3,iistrtor)
#ifdef USE_NVSHMEM_CUDA
               ipe = (k-1)/nbond_pe
               ind = mod((k-1),nbond_pe) +1
               dr  = rcb - d_bl(ipe)%pel(ind)
#else
               dr = rcb - bl(k)
#endif
               e2 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the stretch-torsion values for the third bond
c
               v1 = kst(7,iistrtor)
               v2 = kst(8,iistrtor)
               v3 = kst(9,iistrtor)
               k  = ist(4,iistrtor)
#ifdef USE_NVSHMEM_CUDA
               ipe = (k-1)/nbond_pe
               ind = mod((k-1),nbond_pe) +1
               dr  = rdc - d_bl(ipe)%pel(ind)
#else
               dr = rdc - bl(k)
#endif
               e3 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     increment the total stretch-torsion energy
c
               ebt = ebt + e1 + e2 + e3
            end if
         end if
      end do
      end
