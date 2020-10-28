c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine etors1  --  torsional energy & derivatives  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "etors1" calculates the torsional potential energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      module etors1gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine etors1gpu
      implicit none
      call etors1agpu
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine etors1a  --  standard torsional energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "etors1a" calculates the torsional potential energy and first
c     derivatives with respect to Cartesian coordinates using a
c     standard sum of Fourier terms
c
c
      subroutine etors1agpu
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use etors1gpu_inl
      use group
      use inform   ,only:deb_Path
      use potent   ,only:use_amd_dih
      use tinheader,only:ti_p,re_p
      use torpot
      use tors
      use usage
      use virial
      use mamd
      implicit none
      integer i,ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer itor
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,dedphi,rcb
      real(t_p) v1,v2,v3,v4,v5,v6
      real(t_p) c1,c2,c3,c4,c5,c6
      real(t_p) s1,s2,s3,s4,s5,s6
      real(t_p) cosine,cosine2
      real(t_p) cosine3,cosine4
      real(t_p) cosine5,cosine6
      real(t_p) sine,sine2,sine3
      real(t_p) sine4,sine5,sine6
      real(t_p) phi1,phi2,phi3
      real(t_p) phi4,phi5,phi6
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xdc,ydc,zdc
      real(t_p) xcb,ycb,zcb
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(t_p) xt,yt,zt,rt2
      real(t_p) xu,yu,zu,ru2
      real(t_p) xtu,ytu,ztu,rtru
      real(t_p) dphi1,dphi2,dphi3
      real(t_p) dphi4,dphi5,dphi6
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(r_p) dedi(3)
      logical proceed
!$acc routine(image_acc) seq

      if (deb_Path) write(*,*) 'etors1gpu'
c
c     calculate the torsional angle energy and first derivatives
c
!$acc parallel loop private(dedi) async
!$acc&  present(det,et,torsglob,loc,x,y,z,vir,use,
#ifndef USE_NVSHMEM_CUDA
!$acc&     itors,tors1,tors2,tors3,tors4,tors5,tors6,
#endif
!$acc&          g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do itor = 1, ntorsloc
         i   = torsglob(itor)
#ifdef USE_NVSHMEM_CUDA
         ipe = (i-1)/ntors_pe
         ind = mod((i-1),ntors_pe) +1
         ia  = d_itors(ipe)%pel(1,ind)
         ib  = d_itors(ipe)%pel(2,ind)
         ic  = d_itors(ipe)%pel(3,ind)
         id  = d_itors(ipe)%pel(4,ind)
#else
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
#endif
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
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
c     set the torsional parameters for this angle
c
#ifdef USE_NVSHMEM_CUDA
               v1 = d_tors1(ipe)%pel(1,ind)
               c1 = d_tors1(ipe)%pel(3,ind)
               s1 = d_tors1(ipe)%pel(4,ind)
               v2 = d_tors2(ipe)%pel(1,ind)
               c2 = d_tors2(ipe)%pel(3,ind)
               s2 = d_tors2(ipe)%pel(4,ind)
               v3 = d_tors3(ipe)%pel(1,ind)
               c3 = d_tors3(ipe)%pel(3,ind)
               s3 = d_tors3(ipe)%pel(4,ind)
               v4 = d_tors4(ipe)%pel(1,ind)
               c4 = d_tors4(ipe)%pel(3,ind)
               s4 = d_tors4(ipe)%pel(4,ind)
               v5 = d_tors5(ipe)%pel(1,ind)
               c5 = d_tors5(ipe)%pel(3,ind)
               s5 = d_tors5(ipe)%pel(4,ind)
               v6 = d_tors6(ipe)%pel(1,ind)
               c6 = d_tors6(ipe)%pel(3,ind)
               s6 = d_tors6(ipe)%pel(4,ind)
#else
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
#endif
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine  - sine*sine
               sine2   = 2.0_ti_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3   = cosine*sine2   + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4   = cosine*sine3   + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5   = cosine*sine4   + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6   = cosine*sine5   + sine*cosine5
               phi1  = 1.0_ti_p + (cosine *c1 + sine *s1)
               phi2  = 1.0_ti_p + (cosine2*c2 + sine2*s2)
               phi3  = 1.0_ti_p + (cosine3*c3 + sine3*s3)
               phi4  = 1.0_ti_p + (cosine4*c4 + sine4*s4)
               phi5  = 1.0_ti_p + (cosine5*c5 + sine5*s5)
               phi6  = 1.0_ti_p + (cosine6*c6 + sine6*s6)
               dphi1 =            (cosine *s1 - sine *c1)
               dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
               dphi4 = 4.0_ti_p * (cosine4*s4 - sine4*c4)
               dphi5 = 5.0_ti_p * (cosine5*s5 - sine5*c5)
               dphi6 = 6.0_ti_p * (cosine6*s6 - sine6*c6)
c
c     calculate torsional energy and master chain rule term
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                       + v4*phi4 + v5*phi5 + v6*phi6)
               dedphi = torsunit * (v1*dphi1 + v2*dphi2 + v3*dphi3
     &                            + v4*dphi4 + v5*dphi5 + v6*dphi6)
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
c     increment the total torsional angle energy and gradient
c
               et = et + e

               dedi(1) = dedxib
               dedi(2) = dedyib
               dedi(3) = dedzib
!$acc atomic update
               det(1,ibloc) = det(1,ibloc) + dedi(1)
!$acc atomic update
               det(2,ibloc) = det(2,ibloc) + dedi(2)
!$acc atomic update
               det(3,ibloc) = det(3,ibloc) + dedi(3)
c
               dedi(1) = dedxia
               dedi(2) = dedyia
               dedi(3) = dedzia
!$acc atomic update
               det(1,ialoc) = det(1,ialoc) + dedi(1)
!$acc atomic update
               det(2,ialoc) = det(2,ialoc) + dedi(2)
!$acc atomic update
               det(3,ialoc) = det(3,ialoc) + dedi(3)
c
               dedi(1) = dedxic
               dedi(2) = dedyic
               dedi(3) = dedzic
!$acc atomic update
               det(1,icloc) = det(1,icloc) + dedi(1)
!$acc atomic update
               det(2,icloc) = det(2,icloc) + dedi(2)
!$acc atomic update
               det(3,icloc) = det(3,icloc) + dedi(3)
c
               dedi(1) = dedxid
               dedi(2) = dedyid
               dedi(3) = dedzid
!$acc atomic update
               det(1,idloc) = det(1,idloc) + dedi(1)
!$acc atomic update
               det(2,idloc) = det(2,idloc) + dedi(2)
!$acc atomic update
               det(3,idloc) = det(3,idloc) + dedi(3)
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
c
c     increment the aMD/GaMD derivs table
c
      if (use_amd_dih) then
!$acc parallel loop collapse(2) async
!$acc&         present(deamdD,det)
         do ia = 1,nloc
            do ib = 1,3
               deamdD(ib,ia) = det(ib,ia)
            end do
         end do
!$acc serial async present(eDaMD,et)
         eDaMD = et
!$acc end serial
      end if
      end
