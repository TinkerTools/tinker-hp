c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrtor1  --  stretch-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrtor1" calculates the stretch-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      subroutine estrtor1
      use atmlst
      use atoms
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use group
      use strtor
      use torpot
      use tors
      use tinheader,only:ti_p,re_p
      use potent   ,only:use_amd_dih
      use usage
      use virial
      use mamd
      implicit none
      integer i,k,istrtor,iistrtor
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      real(t_p) e,rcb,dr
      real(t_p) ddr,dedphi
      real(t_p) rt2,ru2,rtru
      real(t_p) ddrdx,ddrdy,ddrdz
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
      real(t_p) dphi1,dphi2,dphi3
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed
c
!$acc wait
!$acc update host(debt,vir)
c
c     zero out the stretch-torsion energy and first derivatives
c
      ebt = 0.0_ti_p
c
c     calculate the stretch-torsion interaction energy term
c
      do istrtor = 1, nstrtorloc
         iistrtor = strtorglob(istrtor)
         i = ist(1,iistrtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
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
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,iistrtor)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = kst(2,iistrtor)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = kst(3,iistrtor)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
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
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
c
c     calculate the bond-stretch for the central bond
c
               k = ist(2,iistrtor)
               dr = rcb - bl(k)
c
c     calculate stretch-torsion energy and chain rule terms
c
               e = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddr = e / (dr * rcb)
c
c     first derivative components for the bond stretch
c
               ddrdx = xcb * ddr
               ddrdy = ycb * ddr
               ddrdz = zcb * ddr
c
c     chain rule terms for first derivative components
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               if (use_polymer) then
                  call image (xca,yca,zca)
                  call image (xdb,ydb,zdb)
               end if
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute derivative components for this interaction
c
               dedxia = zcb*dedyt - ycb*dedzt
               dedyia = xcb*dedzt - zcb*dedxt
               dedzia = ycb*dedxt - xcb*dedyt
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu
     &                     - ydc*dedzu - ddrdx
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu
     &                     - zdc*dedxu - ddrdy
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu
     &                     - xdc*dedyu - ddrdz
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu
     &                     - zdb*dedyu + ddrdx
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu
     &                     - xdb*dedzu + ddrdy
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu
     &                     - ydb*dedxu + ddrdz
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     increment the stretch-torsion energy and gradient
c
               ebt = ebt + e
               debt(1,ibloc) = debt(1,ibloc) + dedxib
               debt(2,ibloc) = debt(2,ibloc) + dedyib
               debt(3,ibloc) = debt(3,ibloc) + dedzib
c
               debt(1,ialoc) = debt(1,ialoc) + dedxia
               debt(2,ialoc) = debt(2,ialoc) + dedyia
               debt(3,ialoc) = debt(3,ialoc) + dedzia
c
               debt(1,icloc) = debt(1,icloc) + dedxic
               debt(2,icloc) = debt(2,icloc) + dedyic
               debt(3,icloc) = debt(3,icloc) + dedzic
c
               debt(1,idloc) = debt(1,idloc) + dedxid
               debt(2,idloc) = debt(2,idloc) + dedyid
               debt(3,idloc) = debt(3,idloc) + dedzid
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
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
     $ debt(:,1:nloc)
      if (use_amd_dih) eDaMD = eDaMD + ebt
!$acc update device(debt,vir)
!$acc data present(ebt)
!$acc update host(ebt)
!$acc end data
      return
      end
