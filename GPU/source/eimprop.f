c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine eimprop  --  improper dihedral energy  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "eimprop" calculates the improper dihedral potential energy
c
c
#include "tinker_macro.h"
      module eimprop_inl
        contains
#include "image.f.inc"
#include "groups.inc.f"
      end module

      subroutine eimprop
      use atmlst
      use atmtyp
      use atoms
      use bound
      use eimprop_inl
      use energi
      use group
      use improp
      use inform
      use math
      use torpot
      use tinheader
      use usage
      implicit none
      integer i,ia,ib,ic,id
      integer iimprop
      real(t_p) e,dt,fgrp
      real(t_p) ideal,force
      real(t_p) cosine,sine
      real(t_p) rcb,angle
      real(t_p) xt,yt,zt
      real(t_p) xu,yu,zu
      real(t_p) xtu,ytu,ztu
      real(t_p) rt2,ru2,rtru
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      integer iga,igb,igc,igd,gmin,gmax
      logical proceed

      if (deb_Path) write(*,*) 'eimprop'
c
c     zero out improper dihedral energy
c
      eid = 0.0_ti_p
c
c     calculate the improper dihedral angle energy term
c
!$acc parallel loop async
!$acc&         default(present) present(eid)
      do iimprop = 1, niproploc
         i = impropglob(iimprop)
         ia = iiprop(1,i)
         ib = iiprop(2,i)
         ic = iiprop(3,i)
         id = iiprop(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
         if (use_group)
     &      call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
c
c     compute the value of the improper dihedral angle
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
c     set the improper dihedral parameters for this angle
c
               ideal = vprop(i)
               force = kprop(i)
               if (abs(angle+ideal) .lt. abs(angle-ideal))
     &            ideal = -ideal
               dt = angle - ideal
               do while (dt .gt. 180.0_ti_p)
                  dt = dt - 360.0_ti_p
               end do
               do while (dt .lt. -180.0_ti_p)
                  dt = dt + 360.0_ti_p
               end do
c
c     calculate the improper dihedral energy
c
               e = idihunit * force * dt**2
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total improper dihedral energy
c
               eid = eid + e
            end if
         end if
      end do
      end
