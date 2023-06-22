c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eimprop1  --  impr. dihedral energy & gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eimprop1" calculates improper dihedral energy and its
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module eimprop1gpu_inl
#include "atomicOp.h.f"
        contains
#include "image.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"
      end module

      subroutine eimprop1gpu
      use atmlst
      use atmtyp
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eimprop1gpu_inl
      use group
      use improp
      use inform
      use math
      use torpot
      use tinheader ,only:ti_p,re_p
      use tintypes
      use usage
      use virial
      implicit none
      integer i,ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer iimprop
      real(t_p) e,dedphi,fgrp
      real(t_p) dt
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
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(r_p) dedx,dedy,dedz
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      integer iga,igb,igc,igd,gmin,gmax
      logical proceed
c
      if (deb_Path) write(*,*) 'eimprop1gpu'
c
c     zero out energy and first derivative components
c
      eid = 0.0_ti_p
c
c     calculate the improper dihedral angle energy term
c
!$acc parallel loop async
!$acc&         present(impropglob,iiprop,kprop,vprop,grplist,wgrp,
!$acc&   loc,use,x,y,z,deid)
!$acc&         present(eid,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         reduction(+:eid,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do iimprop = 1, niproploc
         i     = impropglob(iimprop)
         ia    = iiprop(1,i)
         ib    = iiprop(2,i)
         ic    = iiprop(3,i)
         id    = iiprop(4,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
c
c     decide whether to compute the current interaction
c
         if (use_group)
     &      call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
         proceed = (use(ia) .or. use(ib)
     &                       .or. use(ic) .or. use(id))
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
c     calculate improper energy and master chain rule term
c
               e = idihunit * force * dt**2
               dedphi = 2.0_ti_p * radian * idihunit * force * dt
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
               dedxt =  dedphi * (yt*zcb - ycb*zt) /(rt2*rcb)
               dedyt =  dedphi * (zt*xcb - zcb*xt) /(rt2*rcb)
               dedzt =  dedphi * (xt*ycb - xcb*yt) /(rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) /(ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) /(ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) /(ru2*rcb)
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
c     calculate improper dihedral energy and derivatives
c
               eid   = eid + e

               call atomic_add( deid(1,ialoc),dedxia )
               call atomic_add( deid(2,ialoc),dedyia )
               call atomic_add( deid(3,ialoc),dedzia )
c
               call atomic_add( deid(1,ibloc),dedxib )
               call atomic_add( deid(2,ibloc),dedyib )
               call atomic_add( deid(3,ibloc),dedzib )
c
               call atomic_add( deid(1,icloc),dedxic )
               call atomic_add( deid(2,icloc),dedyic )
               call atomic_add( deid(3,icloc),dedzic )
c
               call atomic_add( deid(1,idloc),dedxid )
               call atomic_add( deid(2,idloc),dedyid )
               call atomic_add( deid(3,idloc),dedzid )
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
