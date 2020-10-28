c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epitors1  --  pi-orbit torsion energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epitors1" calculates the pi-orbital torsion potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      subroutine epitors1
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use group
      use pitors
      use torpot
      use tinheader ,only:ti_p,re_p
      use usage
      use virial
      use potent,only:use_amd_dih
      use mamd
      implicit none
      integer i,ia,ib,ic
      integer id,ie,ig
      integer ialoc,ibloc,icloc,idloc,ieloc,igloc
      integer ipitors
      real(t_p) e,dedphi
      real(t_p) xt,yt,zt,rt2
      real(t_p) xu,yu,zu,ru2
      real(t_p) xtu,ytu,ztu
      real(t_p) rdc,rtru
      real(t_p) v2,c2,s2
      real(t_p) phi2,dphi2
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
      real(t_p) xdp,ydp,zdp
      real(t_p) xqc,yqc,zqc
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(t_p) dedxie,dedyie,dedzie
      real(t_p) dedxig,dedyig,dedzig
      real(t_p) dedxip,dedyip,dedzip
      real(t_p) dedxiq,dedyiq,dedziq
      real(t_p) vxterm,vyterm,vzterm
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed
c
!$acc update host(dept,vir)
c
c     zero out the pi-orbital torsion energy and first derivatives
c
      ept = 0.0_ti_p
c
c     calculate the pi-orbital torsion angle energy term
c
      do ipitors = 1, npitorsloc
         i = pitorsglob(ipitors)
         ia = ipit(1,i)
         ib = ipit(2,i)
         ic = ipit(3,i)
         id = ipit(4,i)
         ie = ipit(5,i)
         ig = ipit(6,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
         ieloc = loc(ie)
         igloc = loc(ig)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic) .or.
     &                              use(id) .or. use(ie) .or. use(ig))
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
               call image (xad,yad,zad)
               call image (xbd,ybd,zbd)
               call image (xec,yec,zec)
               call image (xgc,ygc,zgc)
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
               call image (xcp,ycp,zcp)
               call image (xdc,ydc,zdc)
               call image (xqd,yqd,zqd)
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
               rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru)
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
               sine2 = 2.0_ti_p * cosine * sine
               phi2 = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
c
c     calculate pi-orbital torsion energy and master chain rule term
c
               e = ptorunit * v2 * phi2
               dedphi = ptorunit * v2 * dphi2
c
c     chain rule terms for first derivative components
c
               xdp = xid - xip
               ydp = yid - yip
               zdp = zid - zip
               xqc = xiq - xic
               yqc = yiq - yic
               zqc = ziq - zic
               dedxt = dedphi * (yt*zdc - ydc*zt) / (rt2*rdc)
               dedyt = dedphi * (zt*xdc - zdc*xt) / (rt2*rdc)
               dedzt = dedphi * (xt*ydc - xdc*yt) / (rt2*rdc)
               dedxu = -dedphi * (yu*zdc - ydc*zu) / (ru2*rdc)
               dedyu = -dedphi * (zu*xdc - zdc*xu) / (ru2*rdc)
               dedzu = -dedphi * (xu*ydc - xdc*yu) / (ru2*rdc)
c
c     compute first derivative components for pi-orbital angle
c
               dedxip = zdc*dedyt - ydc*dedzt
               dedyip = xdc*dedzt - zdc*dedxt
               dedzip = ydc*dedxt - xdc*dedyt
               dedxic = ydp*dedzt - zdp*dedyt + zqd*dedyu - yqd*dedzu
               dedyic = zdp*dedxt - xdp*dedzt + xqd*dedzu - zqd*dedxu
               dedzic = xdp*dedyt - ydp*dedxt + yqd*dedxu - xqd*dedyu
               dedxid = zcp*dedyt - ycp*dedzt + yqc*dedzu - zqc*dedyu
               dedyid = xcp*dedzt - zcp*dedxt + zqc*dedxu - xqc*dedzu
               dedzid = ycp*dedxt - xcp*dedyt + xqc*dedyu - yqc*dedxu
               dedxiq = zdc*dedyu - ydc*dedzu
               dedyiq = xdc*dedzu - zdc*dedxu
               dedziq = ydc*dedxu - xdc*dedyu
c
c     compute first derivative components for individual atoms
c
               dedxia = ybd*dedzip - zbd*dedyip
               dedyia = zbd*dedxip - xbd*dedzip
               dedzia = xbd*dedyip - ybd*dedxip
               dedxib = zad*dedyip - yad*dedzip
               dedyib = xad*dedzip - zad*dedxip
               dedzib = yad*dedxip - xad*dedyip
               dedxie = ygc*dedziq - zgc*dedyiq
               dedyie = zgc*dedxiq - xgc*dedziq
               dedzie = xgc*dedyiq - ygc*dedxiq
               dedxig = zec*dedyiq - yec*dedziq
               dedyig = xec*dedziq - zec*dedxiq
               dedzig = yec*dedxiq - xec*dedyiq
               dedxic = dedxic + dedxip - dedxie - dedxig
               dedyic = dedyic + dedyip - dedyie - dedyig
               dedzic = dedzic + dedzip - dedzie - dedzig
               dedxid = dedxid + dedxiq - dedxia - dedxib
               dedyid = dedyid + dedyiq - dedyia - dedyib
               dedzid = dedzid + dedziq - dedzia - dedzib
c
c     increment the total pi-orbital torsion energy and gradient
c
               ept = ept + e
               dept(1,icloc) = dept(1,icloc) + dedxic
               dept(2,icloc) = dept(2,icloc) + dedyic
               dept(3,icloc) = dept(3,icloc) + dedzic
c
               dept(1,ialoc) = dept(1,ialoc) + dedxia
               dept(2,ialoc) = dept(2,ialoc) + dedyia
               dept(3,ialoc) = dept(3,ialoc) + dedzia
c
               dept(1,ibloc) = dept(1,ibloc) + dedxib
               dept(2,ibloc) = dept(2,ibloc) + dedyib
               dept(3,ibloc) = dept(3,ibloc) + dedzib
c
               dept(1,idloc) = dept(1,idloc) + dedxid
               dept(2,idloc) = dept(2,idloc) + dedyid
               dept(3,idloc) = dept(3,idloc) + dedzid
c
               dept(1,ieloc) = dept(1,ieloc) + dedxie
               dept(2,ieloc) = dept(2,ieloc) + dedyie
               dept(3,ieloc) = dept(3,ieloc) + dedzie
c
               dept(1,igloc) = dept(1,igloc) + dedxig
               dept(2,igloc) = dept(2,igloc) + dedyig
               dept(3,igloc) = dept(3,igloc) + dedzig
c
c     increment the internal virial tensor components
c
               vxterm = dedxid + dedxia + dedxib
               vyterm = dedyid + dedyia + dedyib
               vzterm = dedzid + dedzia + dedzib
               vxx = xdc*vxterm + xcp*dedxip - xqd*dedxiq
               vyx = ydc*vxterm + ycp*dedxip - yqd*dedxiq
               vzx = zdc*vxterm + zcp*dedxip - zqd*dedxiq
               vyy = ydc*vyterm + ycp*dedyip - yqd*dedyiq
               vzy = zdc*vyterm + zcp*dedyip - zqd*dedyiq
               vzz = zdc*vzterm + zcp*dedzip - zqd*dedziq
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end if
      end do
c
c     increment the aMD/GaMD derivs table
c
      if (use_amd_dih) deamdD(:,1:nloc) = deamdD(:,1:nloc) +
     $ dept(:,1:nloc)
      if (use_amd_dih) eDaMD = eDaMD + ept
!$acc update device(dept,vir)
      return
      end
