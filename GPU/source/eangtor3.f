c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eangtor3  --  angle-torsion energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eangtor3" calculates the angle-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
#include "tinker_precision.h"
      module eangtor3_inl
      contains
#include "image.f.inc"
      end module

      subroutine eangtor3
      use action
      use analyz
      use angle
      use angtor
      use atmtyp
      use atmlst
      use atoms
      use bound
      use domdec   ,only: loc
      use eangtor3_inl
      use energi
      use group
      use inform
      use iounit
      use math
      use tinheader,only: ti_p
      use torpot
      use tors
      use usage
      implicit none
      integer i,k,iangtor,iiangtor
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      real(t_p) e,e1,e2
      real(t_p) rcb,fgrp
      real(t_p) rt2,ru2,rtru
      real(t_p) rba2,rcb2,rdc2
      real(t_p) dot,dt,tangle
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
      real(t_p) angle1,cosang
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
      if (deb_Path) print*, "eangtor3"
#ifdef USE_NVSHMEM_CUDA
      ! TODO Remove this check
      ! Implement NVSHMEM Access to anat & itors, tors[123]
      print*, '  FATAL ERROR  '
      print*, 'NVSHMEM feature not implemented inside eangtor3'
      call fatal
#endif
c
c     zero out the energy due to extra potential terms
c
      neat = 0
      eat  = 0.0_re_p
      aet  = 0.0_ti_p
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nangtor.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Angle-Torsion Interactions :',
     &           //,' Type',25x,'Atom Names',21x,'Angle',
     &              6x,'Energy',/)
      end if
c
c
c     calculate the angle-torsion interaction energy term
c
!$acc parallel loop async reduction(+:eat,neat)
!$acc&         present(eat,aeat,angtorglob,iat,anat
!$acc&     ,itors,tors1,tors2,tors3,x,y,z,use,kant,loc)
      do iangtor = 1, nangtorloc
         iiangtor = angtorglob(iangtor)
         i = iat(1,iiangtor)
         ia = itors(1,i)
         ialoc = loc(ia)
         ib = itors(2,i)
         ibloc = loc(ib)
         ic = itors(3,i)
         icloc = loc(ic)
         id = itors(4,i)
         idloc = loc(id)
c
c     decide whether to compute the current interaction
c
c         if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
         proceed = (use(ia) .or. use(ib) .or.
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
            rba2 = xba*xba + yba*yba + zba*zba
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
            if (min(rba2,rcb2,rdc2) .ne. 0.0_ti_p) then
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
               rt2 = max(rt2,0.000001_ti_p)
               ru2 = xu*xu + yu*yu + zu*zu
               ru2 = max(ru2,0.000001_ti_p)
               rtru = sqrt(rt2*ru2)
               rcb = sqrt(rcb2)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               tangle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  tangle = -tangle
c
c     compute multiple angle trigonometry and phase terms
c
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0_ti_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0_ti_p + (cosine*c1 + sine*s1)
               phi2 = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               phi3 = 1.0_ti_p + (cosine3*c3 + sine3*s3)
c
c     get the angle-torsion values for the first angle
c
               v1 = kant(1,iiangtor)
               v2 = kant(2,iiangtor)
               v3 = kant(3,iiangtor)
               k = iat(2,iiangtor)
               dot = xba*xcb + yba*ycb + zba*zcb
               cosang = -dot / sqrt(rba2*rcb2)
               angle1 = radian * acos(cosang)
               dt = angle1 - anat(k)
               e1 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the angle-torsion values for the second angle
c
               v1 = kant(4,iiangtor)
               v2 = kant(5,iiangtor)
               v3 = kant(6,iiangtor)
               k = iat(3,iiangtor)
               dot = xcb*xdc + ycb*ydc + zcb*zdc
               cosang = acos(-dot / sqrt(rcb2*rdc2))
               angle1 = radian * (cosang)
               dt = angle1 - anat(k)
               e2 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
cc
cc     scale the interaction based on its group membership
cc
c               if (use_group) then
c                  e1 = e1 * fgrp
c                  e2 = e2 * fgrp
c               end if
c
c     increment the total angle-torsion energy
c
               neat = neat + 1
               e    = e1 + e2
               eat  = eat + e
!$acc atomic
               aeat(ialoc) = aeat(ialoc) + e1/3.0_ti_p
!$acc atomic
               aeat(ibloc) = aeat(ibloc) + e/3.0_ti_p
!$acc atomic
               aeat(icloc) = aeat(icloc) + e/3.0_ti_p
!$acc atomic
               aeat(idloc) = aeat(idloc) + e2/3.0_ti_p
#ifndef _OPENACC
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 3.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Angle-Torsion',
     &                          ' Interactions :',
     &                       //,' Type',25x,'Atom Names',21x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,30)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),tangle,e
   30             format (' AngTors',3x,4(i7,'-',a3),f11.4,f12.4)
               end if
#endif
            end if
         end if
      end do
!$acc update host(aeat) async
      end
