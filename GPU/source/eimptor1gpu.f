c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eimptor1  --  impr. torsion energy & gradient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eimptor1gpu" calculates on device improper torsion energy and its
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module eimptor1gpu_inl
#include "atomicOp.h.f"
      contains
#include "image.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"
      end module

      subroutine eimptor1gpu
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eimptor1gpu_inl
      use group
      use inform    ,only: deb_Path
      use imptor
      use torpot
      use tinheader ,only:ti_p,re_p
      use usage
      use virial
      implicit none
      integer i,ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer iimptor
      real(t_p) e,dedphi,fgrp
      real(t_p) rcb
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
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(t_p) dphi1,dphi2,dphi3
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(r_p) dedx,dedy,dedz
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      integer iga,igb,igc,igd,gmin,gmax
      logical proceed

      if (deb_Path) write(*,*) "eimptor1gpu",nitorsloc
c
c
c     zero out energy and first derivative components
c
c     eit = 0.0_ti_p
c
c     calculate the improper torsional angle energy term
c
!$acc parallel loop async default(present)
!$acc&             present(eit,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         reduction(+:eit,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do iimptor = 1, nitorsloc
         i  = imptorglob(iimptor)
         ia = iitors(1,i)
         ib = iitors(2,i)
         ic = iitors(3,i)
         id = iitors(4,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
c
c     decide whether to compute the current interaction
c
         if (use_group)
     &      call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
         proceed = (use(ia) .or. use(ib) .or.
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
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
c
c     calculate improper torsion energy and master chain rule term
c
               e = itorunit * (v1*phi1+v2*phi2+v3*phi3)
               dedphi = itorunit * (v1*dphi1+v2*dphi2+v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
               end if
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
                  call image_inl (xca,yca,zca)
                  call image_inl (xdb,ydb,zdb)
               end if
               dedxt =  dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt =  dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt =  dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute first derivative components for this angle
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
c     increment the improper torsion energy and gradient
c
               eit = eit + e

               deit(1,icloc) = deit(1,icloc) + dedxic
               deit(2,icloc) = deit(2,icloc) + dedyic
               deit(3,icloc) = deit(3,icloc) + dedzic
c
               deit(1,ialoc) = deit(1,ialoc) + dedxia
               deit(2,ialoc) = deit(2,ialoc) + dedyia
               deit(3,ialoc) = deit(3,ialoc) + dedzia
c
               deit(1,ibloc) = deit(1,ibloc) + dedxib
               deit(2,ibloc) = deit(2,ibloc) + dedyib
               deit(3,ibloc) = deit(3,ibloc) + dedzib
c
               deit(1,idloc) = deit(1,idloc) + dedxid
               deit(2,idloc) = deit(2,idloc) + dedyid
               deit(3,idloc) = deit(3,idloc) + dedzid
c
c     increment the internal virial tensor components
c
           g_vxx = g_vxx + xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
           g_vxy = g_vxy + ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
           g_vxz = g_vxz + zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
           g_vyy = g_vyy + ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
           g_vyz = g_vyz + zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
           g_vzz = g_vzz + zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
            end if
         end if
      end do
      end
