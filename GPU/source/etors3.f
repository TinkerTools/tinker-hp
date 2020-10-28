c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine etors3  --  torsional energy & analysis  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "etors3" calculates the torsional potential energy; also
c     partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      module etors3_inl
        contains
#include "image.f.inc"
      end module

      subroutine etors3
      implicit none
      call etors3a
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine etors3a  --  standard torsional analysis  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "etors3a" calculates the torsional potential energy using
c     a standard sum of Fourier terms and partitions the energy
c     among the atoms
c
c
      subroutine etors3a
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use domdec
      use energi
      use etors3_inl
      use group
      use inform
      use iounit
      use math
      use tinheader ,only:ti_p,re_p
      use torpot
      use tors
      use usage
      implicit none
      integer i,ia,ib,ic,id,ibloc,icloc
      integer itor
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,rcb,angle
      real(t_p) xt,yt,zt,rt2
      real(t_p) xu,yu,zu,ru2
      real(t_p) xtu,ytu,ztu,rtru
      real(t_p) v1,v2,v3,v4,v5,v6
      real(t_p) c1,c2,c3,c4,c5,c6
      real(t_p) s1,s2,s3,s4,s5,s6
      real(t_p) cosine,cosine2
      real(t_p) cosine3,cosine4
      real(t_p) cosine5,cosine6
      real(t_p) sine,sine2,sine3
      real(t_p) sine4,sine5,sine6
      real(t_p) phi1,phi2,phi3
      real(t_p) phi4,phi5,phi6
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xdc,ydc,zdc
      real(t_p) xcb,ycb,zcb
      logical proceed
      logical header,huge
c
c
c     zero out the torsional energy and partitioning terms
c
      if(deb_Path) write(*,*) 'etors3'
      net = 0
      et = 0.0_ti_p
      aet = 0.0_ti_p
      if (rank.eq.0) then
         header = .true.
      else
         header=.false.
      end if
c
c     calculate the torsional angle energy term
c
!$acc parallel loop default(present) present(et) async
      do itor = 1, ntorsloc
         i = torsglob(itor)
#ifdef USE_NVSHMEM_CUDA
         ipe = (i-1)/ntors_pe
         ind = mod((i-1),ntors_pe) +1
         ia  = d_itors(ipe)%pel(1,ind)
         ib  = d_itors(ipe)%pel(2,ind)
         ic  = d_itors(ipe)%pel(3,ind)
         id  = d_itors(ipe)%pel(4,ind)
#else
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
#endif
         ibloc = loc(ib)
         icloc = loc(ic)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
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
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0_ti_p) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0_ti_p)  angle = -angle
c
c     set the torsional parameters for this angle
c
#ifdef USE_NVSHMEM_CUDA
               v1 = d_tors1(ipe)%pel(1,ind)
               c1 = d_tors1(ipe)%pel(3,ind)
               s1 = d_tors1(ipe)%pel(4,ind)
               v2 = d_tors2(ipe)%pel(1,ind)
               c2 = d_tors2(ipe)%pel(3,ind)
               s2 = d_tors2(ipe)%pel(4,ind)
               v3 = d_tors3(ipe)%pel(1,ind)
               c3 = d_tors3(ipe)%pel(3,ind)
               s3 = d_tors3(ipe)%pel(4,ind)
               v4 = d_tors4(ipe)%pel(1,ind)
               c4 = d_tors4(ipe)%pel(3,ind)
               s4 = d_tors4(ipe)%pel(4,ind)
               v5 = d_tors5(ipe)%pel(1,ind)
               c5 = d_tors5(ipe)%pel(3,ind)
               s5 = d_tors5(ipe)%pel(4,ind)
               v6 = d_tors6(ipe)%pel(1,ind)
               c6 = d_tors6(ipe)%pel(3,ind)
               s6 = d_tors6(ipe)%pel(4,ind)
#else
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
#endif
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2   = 2.0_ti_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3   = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4   = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5   = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6   = cosine*sine5 + sine*cosine5
               phi1 = 1.0_ti_p + (cosine *c1 + sine *s1)
               phi2 = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               phi3 = 1.0_ti_p + (cosine3*c3 + sine3*s3)
               phi4 = 1.0_ti_p + (cosine4*c4 + sine4*s4)
               phi5 = 1.0_ti_p + (cosine5*c5 + sine5*s5)
               phi6 = 1.0_ti_p + (cosine6*c6 + sine6*s6)
c
c     calculate the torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                       + v4*phi4 + v5*phi5 + v6*phi6)
c
c     increment the total torsional angle energy
c
               net = net + 1
               et = et + e
!$acc atomic
               aet(ibloc) = aet(ibloc) + 0.5*e
!$acc atomic
               aet(icloc) = aet(icloc) + 0.5*e
#ifndef _OPENACC
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 3.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Torsional Angle',
     &                          ' Interactions :',
     &                       //,' Type',25x,'Atom Names',21x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   20             format (' Torsion',3x,4(i7,'-',a3),f11.4,f12.4)
               end if
#endif
            end if
         end if
      end do
!$acc update self(aet) async
      end
