c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epitors  --  pi-orbit torsion potential energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epitors" calculates the pi-orbital torsion potential energy
c
c
#include "tinker_macro.h"
      module epitors_inl
        contains
#include "image.f.inc"
#include "groups.inc.f"
      end module

      subroutine epitors
      use atmlst
      use atoms
      use bound
      use energi
      use epitors_inl
      use group
      use pitors
      use torpot
      use tinheader
      use usage
      implicit none
      integer i,ia,ib,ic
      integer id,ie,ig,ipitors
      real(t_p) e,rdc
      real(t_p) xt,yt,zt,rt2
      real(t_p) xu,yu,zu,ru2
      real(t_p) xtu,ytu,ztu,rtru
      real(t_p) v2,c2,s2,phi2
      real(t_p) sine,cosine
      real(t_p) sine2,cosine2
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xie,yie,zie
      real(t_p) xig,yig,zig
      real(t_p) xip,yip,zip
      real(t_p) xiq,yiq,ziq
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xec,yec,zec
      real(t_p) xgc,ygc,zgc
      real(t_p) xcp,ycp,zcp
      real(t_p) xdc,ydc,zdc
      real(t_p) xqd,yqd,zqd
      real(t_p) fgrp
      integer iga,igb,igc,igd,ige,igg,gmin,gmax
      logical proceed
c
c
c     zero out the pi-orbital torsion potential energy
c
      ept = 0.0_ti_p
c
c     calculate the pi-orbital torsion angle energy term
c
!$acc parallel loop default(present) present(ept) async
      do ipitors = 1, npitorsloc
         i = pitorsglob(ipitors)
         ia = ipit(1,i)
         ib = ipit(2,i)
         ic = ipit(3,i)
         id = ipit(4,i)
         ie = ipit(5,i)
         ig = ipit(6,i)
c
c     decide whether to compute the current interaction
c
         if (use_group)
     &      call groups6_inl(fgrp,ia,ib,ic,id,ie,ig,ngrp,grplist,wgrp)
         proceed = (use(ia) .or. use(ib) .or. use(ic) .or.
     &              use(id) .or. use(ie) .or. use(ig))
c
c     compute the value of the pi-orbital torsion angle
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
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
            xig = x(ig)
            yig = y(ig)
            zig = z(ig)
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xbd = xib - xid
            ybd = yib - yid
            zbd = zib - zid
            xec = xie - xic
            yec = yie - yic
            zec = zie - zic
            xgc = xig - xic
            ygc = yig - yic
            zgc = zig - zic
            if (use_polymer) then
               call image_inl (xad,yad,zad)
               call image_inl (xbd,ybd,zbd)
               call image_inl (xec,yec,zec)
               call image_inl (xgc,ygc,zgc)
            end if
            xip = yad*zbd - ybd*zad + xic
            yip = zad*xbd - zbd*xad + yic
            zip = xad*ybd - xbd*yad + zic
            xiq = yec*zgc - ygc*zec + xid
            yiq = zec*xgc - zgc*xec + yid
            ziq = xec*ygc - xgc*yec + zid
            xcp = xic - xip
            ycp = yic - yip
            zcp = zic - zip
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xqd = xiq - xid
            yqd = yiq - yid
            zqd = ziq - zid
            if (use_polymer) then
               call image_inl (xcp,ycp,zcp)
               call image_inl (xdc,ydc,zdc)
               call image_inl (xqd,yqd,zqd)
            end if
            xt = ycp*zdc - ydc*zcp
            yt = zcp*xdc - zdc*xcp
            zt = xcp*ydc - xdc*ycp
            xu = ydc*zqd - yqd*zdc
            yu = zdc*xqd - zqd*xdc
            zu = xdc*yqd - xqd*ydc
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0_ti_p) then
               rdc    = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine   = (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru)
c
c     set the pi-orbital torsion parameters for this angle
c
               v2 = kpit(i)
               c2 = -1.0_ti_p
               s2 = 0.0_ti_p
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2   = 2.0_ti_p * cosine * sine
               phi2    = 1.0_ti_p + (cosine2*c2 + sine2*s2)
c
c     calculate the pi-orbital torsion energy for this angle
c
               e = ptorunit * v2 * phi2
               if (use_group) e = e * fgrp
c
c     increment the total pi-orbital torsion angle energy
c
               ept = ept + e
            end if
         end if
      end do
      end
