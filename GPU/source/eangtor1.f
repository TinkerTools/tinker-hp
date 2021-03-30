c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine eangtor1  --  angle-torsion energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "eangtor1" calculates the angle-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      module eangtor1_inl
      contains
#include "image.f.inc"
      end module

      subroutine eangtor1
      use atmlst
      use angle
      use angtor
      use atoms
      use bound
      use deriv
      use domdec
      use eangtor1_inl
      use energi
      use group
      use inform ,only: deb_Path
      use math
      use tinheader,only: ti_p
      use torpot
      use tors
      use usage
      use virial
      implicit none
      integer i,k,iangtor
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer iiangtor
      real(t_p) e,e1,e2
      real(t_p) rcb,fgrp
      real(t_p) ddt,dedphi
      real(t_p) rt2,ru2,rtru
      real(t_p) rba2,rcb2,rdc2
      real(t_p) dot,dt
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
      real(t_p) angle1,cosang
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(t_p) terma,termb
      real(t_p) termc,termd
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(r_p) dedxia,dedyia,dedzia
      real(r_p) dedxib,dedyib,dedzib
      real(r_p) dedxic,dedyic,dedzic
      real(r_p) dedxid,dedyid,dedzid
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed
c
      if (deb_Path) print*, "eangtor1"
#ifdef USE_NVSHMEM_CUDA
      ! TODO Remove this check
      ! Implement NVSHMEM Access to anat & itors, tors[123]
      print*, '  FATAL ERROR  '
      print*, 'NVSHMEM feature not implemented inside eangtor1'
      call fatal
#endif
c
c     zero out the angle-torsion energy and first derivatives
c
      eat = 0.0d0
c
c     calculate the angle-torsion energy and first derviatives
c
!$acc parallel loop async
!$acc&         present(eat,deat,angtorglob,iat,anat
!$acc&     ,itors,tors1,tors2,tors3,x,y,z,use,kant,loc)
      do iangtor = 1, nangtorloc
         iiangtor = angtorglob(iangtor)
         i  = iat(1,iiangtor)
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
               rcb = sqrt(rcb2)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
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
               sine2   = 2.0_ti_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3   = cosine*sine2 + sine*cosine2
               phi1  = 1.0_ti_p + (cosine*c1 + sine*s1)
               phi2  = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               phi3  = 1.0_ti_p + (cosine3*c3 + sine3*s3)
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
c
c     set the angle-torsion parameters for the first angle
c
               v1 = kant(1,iiangtor)
               v2 = kant(2,iiangtor)
               v3 = kant(3,iiangtor)
               k  = iat(2,iiangtor)
               dot = xba*xcb + yba*ycb + zba*zcb
               cosang = -dot / sqrt(rba2*rcb2)
               angle1 = radian * acos(cosang)
               dt = angle1 - anat(k)
               e1 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = atorunit * dt * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
cc
cc     scale the interaction based on its group membership
cc
c               if (use_group) then
c                  e1 = e1 * fgrp
c                  dedphi = dedphi * fgrp
c                  ddt = ddt * fgrp
c               end if
c
c     compute derivative components for this interaction
c
               dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
               dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
               dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
               dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c     increment chain rule components for the first angle
c
               terma = -ddt / (rba2*sqrt(rt2))
               termc = ddt / (rcb2*sqrt(rt2))
               dedxia = terma*(zba*yt-yba*zt) + zcb*dedyt - ycb*dedzt
               dedyia = terma*(xba*zt-zba*xt) + xcb*dedzt - zcb*dedxt
               dedzia = terma*(yba*xt-xba*yt) + ycb*dedxt - xcb*dedyt
               dedxib = terma*(yba*zt-zba*yt) + termc*(zcb*yt-ycb*zt)
     &                     + yca*dedzt - zca*dedyt
     &                     + zdc*dedyu - ydc*dedzu
               dedyib = terma*(zba*xt-xba*zt) + termc*(xcb*zt-zcb*xt)
     &                     + zca*dedxt - xca*dedzt
     &                     + xdc*dedzu - zdc*dedxu
               dedzib = terma*(xba*yt-yba*xt) + termc*(ycb*xt-xcb*yt)
     &                     + xca*dedyt - yca*dedxt
     &                     + ydc*dedxu - xdc*dedyu
               dedxic = termc*(ycb*zt-zcb*yt) + zba*dedyt
     &                     - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = termc*(zcb*xt-xcb*zt) + xba*dedzt
     &                     - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = termc*(xcb*yt-ycb*xt) + yba*dedxt
     &                     - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     get the angle-torsion values for the second angle
c
               v1 = kant(4,iiangtor)
               v2 = kant(5,iiangtor)
               v3 = kant(6,iiangtor)
               k = iat(3,iiangtor)
               dot = xcb*xdc + ycb*ydc + zcb*zdc
               cosang = acos(-dot / sqrt(rcb2*rdc2))
               angle1 = radian * cosang
               dt = angle1 - anat(k)
               e2 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = atorunit * dt * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
cc
cc     scale the interaction based on its group membership
cc
c               if (use_group) then
c                  e2 = e2 * fgrp
c                  dedphi = dedphi * fgrp
c                  ddt = ddt * fgrp
c               end if
c
c     compute derivative components for this interaction
c
               dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
               dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
               dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
               dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c     increment chain rule components for the second angle
c
               termb = -ddt / (rcb2*sqrt(ru2))
               termd = ddt / (rdc2*sqrt(ru2))
               dedxia = dedxia + zcb*dedyt - ycb*dedzt
               dedyia = dedyia + xcb*dedzt - zcb*dedxt
               dedzia = dedzia + ycb*dedxt - xcb*dedyt
               dedxib = dedxib + termb*(zcb*yu-ycb*zu) + yca*dedzt
     &                     - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyib = dedyib + termb*(xcb*zu-zcb*xu) + zca*dedxt
     &                     - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzib = dedzib + termb*(ycb*xu-xcb*yu) + xca*dedyt
     &                     - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxic = dedxic + termb*(ycb*zu-zcb*yu)
     &                     + termd*(zdc*yu-ydc*zu) + zba*dedyt
     &                     - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = dedyic + termb*(zcb*xu-xcb*zu)
     &                     + termd*(xdc*zu-zdc*xu) + xba*dedzt
     &                     - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = dedzic + termb*(xcb*yu-ycb*xu)
     &                     + termd*(ydc*xu-xdc*yu) + yba*dedxt
     &                     - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = dedxid + termd*(ydc*zu-zdc*yu)
     &                     + zcb*dedyu - ycb*dedzu
               dedyid = dedyid + termd*(zdc*xu-xdc*zu)
     &                     + xcb*dedzu - zcb*dedxu
               dedzid = dedzid + termd*(xdc*yu-ydc*xu)
     &                     + ycb*dedxu - xcb*dedyu
c
c     increment the angle-torsion energy and gradient
c
               e   = e1 + e2
               eat = eat + e
!$acc atomic
               deat(1,ialoc) = deat(1,ialoc) + dedxia
!$acc atomic
               deat(2,ialoc) = deat(2,ialoc) + dedyia
!$acc atomic
               deat(3,ialoc) = deat(3,ialoc) + dedzia
!$acc atomic
               deat(1,ibloc) = deat(1,ibloc) + dedxib
!$acc atomic
               deat(2,ibloc) = deat(2,ibloc) + dedyib
!$acc atomic
               deat(3,ibloc) = deat(3,ibloc) + dedzib
!$acc atomic
               deat(1,icloc) = deat(1,icloc) + dedxic
!$acc atomic
               deat(2,icloc) = deat(2,icloc) + dedyic
!$acc atomic
               deat(3,icloc) = deat(3,icloc) + dedzic
!$acc atomic
               deat(1,idloc) = deat(1,idloc) + dedxid
!$acc atomic
               deat(2,idloc) = deat(2,idloc) + dedyid
!$acc atomic
               deat(3,idloc) = deat(3,idloc) + dedzid
c
c     increment the internal virial tensor components
c
               if (use_virial) then
            g_vxx =g_vxx + xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
            g_vxy =g_vxy + ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
            g_vxz =g_vxz + zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
            g_vyy =g_vyy + ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
            g_vyz =g_vyz + zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
            g_vzz =g_vzz + zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
               end if
            end if
         end if
      end do
      end
