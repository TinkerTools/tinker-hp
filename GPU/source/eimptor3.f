c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eimptor3  --  impr. torsion energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eimptor3" calculates the improper torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
#include "tinker_precision.h"
      subroutine eimptor3
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use domdec
      use energi
      use group
      use imptor
      use inform
      use iounit
      use math
      use tinheader ,only:ti_p,re_p
      use torpot
      use usage
      implicit none
      integer i,ia,ib,ic,id
      integer iimptor
      integer ibloc,icloc
      real(t_p) e,rcb
      real(t_p) angle
      real(t_p) xt,yt,zt
      real(t_p) xu,yu,zu
      real(t_p) xtu,ytu,ztu
      real(t_p) rt2,ru2,rtru
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
      logical proceed
      logical header,huge
c
c     zero out the torsional energy and partitioning terms
c
      neit = 0
      eit = 0.0_ti_p
      aeit = 0.0_ti_p
      header = .true.
c
c     calculate the improper torsional angle energy term
c
      do iimptor = 1, nitorsloc
         i = imptorglob(iimptor)
         ia = iitors(1,i)
         ib = iitors(2,i)
         ibloc = loc(ib)
         ic = iitors(3,i)
         icloc = loc(ic)
         id = iitors(4,i)
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
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
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
c     set the improper torsional parameters for this angle
c
               v1 = itors1(1,i)
               c1 = itors1(3,i)
               s1 = itors1(4,i)
               v2 = itors2(1,i)
               c2 = itors2(3,i)
               s2 = itors2(4,i)
               v3 = itors3(1,i)
               c3 = itors3(3,i)
               s3 = itors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0_ti_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0_ti_p + (cosine*c1 + sine*s1)
               phi2 = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               phi3 = 1.0_ti_p + (cosine3*c3 + sine3*s3)
c
c     calculate the improper torsional energy for this angle
c
               e = itorunit * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     increment the total torsional angle energy
c
               neit = neit + 1
               eit = eit + e
               aeit(ibloc) = aeit(ibloc) + 0.5_ti_p*e
               aeit(icloc) = aeit(icloc) + 0.5_ti_p*e
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 5.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Improper Torsion',
     &                          ' Interactions :',
     &                       //,' Type',25x,'Atom Names',21x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   20             format (' Improper',2x,4(i7,'-',a3),f11.4,f12.4)
               end if
            end if
         end if
      end do
      return
      end
