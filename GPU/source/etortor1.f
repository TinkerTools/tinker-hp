c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine etortor1  --  torsion-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "etortor1" calculates the torsion-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine etortor1
      use sizes
      use atmlst
      use atmtyp ,only: atomic
      use atoms
      use bitor
      use bound
      use deriv
      use domdec
      use energi
      use group
      use ktrtor
      use mamd
      use math
      use potent   ,only:use_amd_dih
      use tinheader ,only:ti_p,re_p
      use torpot
      use tortor
      use usage
      use virial
      implicit none
      integer i,k,itortor,iitortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
      integer ialoc,ibloc,icloc,idloc,ieloc
      integer nlo,nhi,nt
      integer xlo,ylo
      real(t_p) e,sign
      real(t_p) angle1,angle2
      real(t_p) value1,value2
      real(t_p) cosine1,cosine2
      real(t_p) xt,yt,zt,rt2
      real(t_p) xu,yu,zu,ru2
      real(t_p) xv,yv,zv,rv2
      real(t_p) rtru
      real(t_p) rurv
      real(t_p) x1l,x1u
      real(t_p) y1l,y1u
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xie,yie,zie
      real(t_p) xba,yba,zba
      real(t_p) xdc,ydc,zdc
      real(t_p) xcb,ycb,zcb
      real(t_p) xed,yed,zed
      real(t_p) rcb,rdc
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(t_p) xec,yec,zec
      real(t_p) dedang1,dedang2
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(t_p) dedxu2,dedyu2,dedzu2
      real(t_p) dedxv2,dedyv2,dedzv2
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(t_p) dedxib2,dedyib2,dedzib2
      real(t_p) dedxic2,dedyic2,dedzic2
      real(t_p) dedxid2,dedyid2,dedzid2
      real(t_p) dedxie2,dedyie2,dedzie2
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      real(t_p) vxx2,vyy2,vzz2
      real(t_p) vyx2,vzx2,vzy2
      real(t_p) ftt(4),ft12(4)
      real(t_p) ft1(4),ft2(4)
      real(t_p) fgrp
      logical proceed
c
!$acc update host(dett,vir)
c
c     zero out the torsion-torsion energy and first derivatives
c
      ett = 0.0_ti_p
c
c     calculate the torsion-torsion interaction energy term
c
      do itortor = 1, ntortorloc
         iitortor = tortorglob(itortor)
         i = itt(1,iitortor)
         k = itt(2,iitortor)
         if (itt(3,iitortor) .eq. 1) then
            ia = ibitor(1,i)
            ib = ibitor(2,i)
            ic = ibitor(3,i)
            id = ibitor(4,i)
            ie = ibitor(5,i)
         else
            ia = ibitor(5,i)
            ib = ibitor(4,i)
            ic = ibitor(3,i)
            id = ibitor(2,i)
            ie = ibitor(1,i)
         end if
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
         ieloc = loc(ie)
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ib,ic,id,ie,0)
         proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                 .or. use(id) .or. use(ie))
c
c     compute the values of the torsional angles
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xed = xie - xid
            yed = yie - yid
            zed = zie - zid
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
               call image (xed,yed,zed)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            xv = ydc*zed - yed*zdc
            yv = zdc*xed - zed*xdc
            zv = xdc*yed - xed*ydc
            rv2 = xv*xv + yv*yv + zv*zv
            rurv = sqrt(ru2 * rv2)
            if (rtru.ne.0.0_ti_p .and. rurv.ne.0.0_ti_p) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
               cosine1 = min(1.0_ti_p,max(-1.0_ti_p,cosine1))
               angle1 = radian * acos(cosine1)
               sign = xba*xu + yba*yu + zba*zu
               if (sign .lt. 0.0_ti_p)  angle1 = -angle1
               value1 = angle1
               rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
               cosine2 = min(1.0_ti_p,max(-1.0_ti_p,cosine2))
               angle2 = radian * acos(cosine2)
               sign = xcb*xv + ycb*yv + zcb*zv
               if (sign .lt. 0.0_ti_p)  angle2 = -angle2
               value2 = angle2
c
c     check for inverted chirality at the central atom
c
               call chkttor (ib,ic,id,sign,value1,value2,x,y,z,
     &                       type,atomic)
c
c     use bicubic interpolation to compute spline values
c
               nlo = 1
               nhi = tnx(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi+nlo) / 2
                  if (ttx(nt,k) .gt. value1) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               xlo = nlo
               nlo = 1
               nhi = tny(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi + nlo)/2
                  if (tty(nt,k) .gt. value2) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               ylo = nlo
               x1l = ttx(xlo,k)
               x1u = ttx(xlo+1,k)
               y1l = tty(ylo,k)
               y1u = tty(ylo+1,k)
               pos2 = ylo*tnx(k) + xlo
               pos1 = pos2 - tnx(k)
               ftt(1) = tbf(pos1,k)
               ftt(2) = tbf(pos1+1,k)
               ftt(3) = tbf(pos2+1,k)
               ftt(4) = tbf(pos2,k)
               ft1(1) = tbx(pos1,k)
               ft1(2) = tbx(pos1+1,k)
               ft1(3) = tbx(pos2+1,k)
               ft1(4) = tbx(pos2,k)
               ft2(1) = tby(pos1,k)
               ft2(2) = tby(pos1+1,k)
               ft2(3) = tby(pos2+1,k)
               ft2(4) = tby(pos2,k)
               ft12(1) = tbxy(pos1,k)
               ft12(2) = tbxy(pos1+1,k)
               ft12(3) = tbxy(pos2+1,k)
               ft12(4) = tbxy(pos2,k)
               call bcuint1 (ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u,
     &                       value1,value2,e,dedang1,dedang2)
               e = ttorunit * e
               dedang1 = sign * ttorunit * radian * dedang1
               dedang2 = sign * ttorunit * radian * dedang2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedang1 = dedang1 * fgrp
                  dedang2 = dedang2 * fgrp
               endif
c
c     chain rule terms for first angle derivative components
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
               dedxt = dedang1 * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedang1 * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedang1 * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedang1 * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedang1 * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedang1 * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute first derivative components for first angle
c
               dedxia = zcb*dedyt - ycb*dedzt
               dedyia = xcb*dedzt - zcb*dedxt
               dedzia = ycb*dedxt - xcb*dedyt
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     chain rule terms for second angle derivative components
c
               xec = xie - xic
               yec = yie - yic
               zec = zie - zic
               if (use_polymer) then
                  call image (xdb,ydb,zdb)
                  call image (xec,yec,zec)
               end if
               dedxu2 = dedang2 * (yu*zdc - ydc*zu) / (ru2*rdc)
               dedyu2 = dedang2 * (zu*xdc - zdc*xu) / (ru2*rdc)
               dedzu2 = dedang2 * (xu*ydc - xdc*yu) / (ru2*rdc)
               dedxv2 = -dedang2 * (yv*zdc - ydc*zv) / (rv2*rdc)
               dedyv2 = -dedang2 * (zv*xdc - zdc*xv) / (rv2*rdc)
               dedzv2 = -dedang2 * (xv*ydc - xdc*yv) / (rv2*rdc)
c
c     compute first derivative components for second angle
c
               dedxib2 = zdc*dedyu2 - ydc*dedzu2
               dedyib2 = xdc*dedzu2 - zdc*dedxu2
               dedzib2 = ydc*dedxu2 - xdc*dedyu2
               dedxic2 = ydb*dedzu2 - zdb*dedyu2
     &                      + zed*dedyv2 - yed*dedzv2
               dedyic2 = zdb*dedxu2 - xdb*dedzu2
     &                      + xed*dedzv2 - zed*dedxv2
               dedzic2 = xdb*dedyu2 - ydb*dedxu2
     &                      + yed*dedxv2 - xed*dedyv2
               dedxid2 = zcb*dedyu2 - ycb*dedzu2
     &                      + yec*dedzv2 - zec*dedyv2
               dedyid2 = xcb*dedzu2 - zcb*dedxu2
     &                      + zec*dedxv2 - xec*dedzv2
               dedzid2 = ycb*dedxu2 - xcb*dedyu2
     &                      + xec*dedyv2 - yec*dedxv2
               dedxie2 = zdc*dedyv2 - ydc*dedzv2
               dedyie2 = xdc*dedzv2 - zdc*dedxv2
               dedzie2 = ydc*dedxv2 - xdc*dedyv2
c
c     increment the torsion-torsion energy and gradient
c
               ett = ett + e
               dett(1,ibloc) = dett(1,ibloc) + dedxib + dedxib2
               dett(2,ibloc) = dett(2,ibloc) + dedyib + dedyib2
               dett(3,ibloc) = dett(3,ibloc) + dedzib + dedzib2
c
               dett(1,ialoc) = dett(1,ialoc) + dedxia
               dett(2,ialoc) = dett(2,ialoc) + dedyia
               dett(3,ialoc) = dett(3,ialoc) + dedzia
c
               dett(1,icloc) = dett(1,icloc) + dedxic + dedxic2
               dett(2,icloc) = dett(2,icloc) + dedyic + dedyic2
               dett(3,icloc) = dett(3,icloc) + dedzic + dedzic2
c
               dett(1,idloc) = dett(1,idloc) + dedxid + dedxid2
               dett(2,idloc) = dett(2,idloc) + dedyid + dedyid2
               dett(3,idloc) = dett(3,idloc) + dedzid + dedzid2
c
               dett(1,ieloc) = dett(1,ieloc) + dedxie2
               dett(2,ieloc) = dett(2,ieloc) + dedyie2
               dett(3,ieloc) = dett(3,ieloc) + dedzie2
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
               vxx2 = xdc*(dedxid2+dedxie2) - xcb*dedxib2 + xed*dedxie2
               vyx2 = ydc*(dedxid2+dedxie2) - ycb*dedxib2 + yed*dedxie2
               vzx2 = zdc*(dedxid2+dedxie2) - zcb*dedxib2 + zed*dedxie2
               vyy2 = ydc*(dedyid2+dedyie2) - ycb*dedyib2 + yed*dedyie2
               vzy2 = zdc*(dedyid2+dedyie2) - zcb*dedyib2 + zed*dedyie2
               vzz2 = zdc*(dedzid2+dedzie2) - zcb*dedzib2 + zed*dedzie2
               vir(1,1) = vir(1,1) + vxx + vxx2
               vir(2,1) = vir(2,1) + vyx + vyx2
               vir(3,1) = vir(3,1) + vzx + vzx2
               vir(1,2) = vir(1,2) + vyx + vyx2
               vir(2,2) = vir(2,2) + vyy + vyy2
               vir(3,2) = vir(3,2) + vzy + vzy2
               vir(1,3) = vir(1,3) + vzx + vzx2
               vir(2,3) = vir(2,3) + vzy + vzy2
               vir(3,3) = vir(3,3) + vzz + vzz2
            end if
         end if
      end do
!$acc update device(dett,vir)
c
c     increment the aMD/GaMD derivs table
c
      if (use_amd_dih) then
         deamdD(:,1:nloc) = deamdD(:,1:nloc) + dett(:,1:nloc)
         eDaMD = eDaMD + ett
!$acc update device(deamdD,eDaMD)
      end if
      return
      end
