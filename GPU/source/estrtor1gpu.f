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
      module estrtor1gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine estrtor1gpu
      use atmlst
      use atoms
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use estrtor1gpu_inl
      use group
      use mamd
      use nvshmem
      use potent   ,only:use_amd_dih
      use strtor
      use tinheader
      use torpot
      use tors
      use usage
      use virial
      implicit none
      integer i,k,istrtor,iistrtor
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(r_p) e
      real(r_p) rcb,dr
      real(r_p) ddr,dedphi
      real(r_p) rt2,ru2,rtru
      real(r_p) ddrdx,ddrdy,ddrdz
      real(r_p) xt,yt,zt
      real(r_p) xu,yu,zu
      real(r_p) xtu,ytu,ztu
      real(r_p) v1,v2,v3
      real(r_p) c1,c2,c3
      real(r_p) s1,s2,s3
      real(r_p) sine,cosine
      real(r_p) sine2,cosine2
      real(r_p) sine3,cosine3
      real(r_p) phi1,phi2,phi3
      real(r_p) dphi1,dphi2,dphi3
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(r_p) dedxt,dedyt,dedzt
      real(r_p) dedxu,dedyu,dedzu
      real(r_p) dedxia,dedyia,dedzia
      real(r_p) dedxib,dedyib,dedzib
      real(r_p) dedxic,dedyic,dedzic
      real(r_p) dedxid,dedyid,dedzid
      real(r_p) dedi(3)
      logical proceed
c
c     calculate the stretch-torsion interaction energy term
c
!$acc parallel loop private(dedi)
!$acc&     present(ebt,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) async
#ifdef USE_NVSHMEM_CUDA
!$acc&     default(present) deviceptr(d_bl)
#else
!$acc&     present(debt,bl,ist,strtorglob,itors,loc,x,y,z)
#endif
      do istrtor = 1, nstrtorloc
         iistrtor = strtorglob(istrtor)
         i     = ist(1,iistrtor)
#ifdef USE_NVSHMEM_CUDA
         ipe   = (i-1)/ntors_pe
         ind   = mod((i-1),ntors_pe) +1
         ia    = d_itors(ipe)%pel(1,ind)
         ib    = d_itors(ipe)%pel(2,ind)
         ic    = d_itors(ipe)%pel(3,ind)
         id    = d_itors(ipe)%pel(4,ind)
#else
         ia    = itors(1,i)
         ib    = itors(2,i)
         ic    = itors(3,i)
         id    = itors(4,i)
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
            xt  = yba*zcb - ycb*zba
            yt  = zba*xcb - zcb*xba
            zt  = xba*ycb - xcb*yba
            xu  = ycb*zdc - ydc*zcb
            yu  = zcb*xdc - zdc*xcb
            zu  = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0_re_p) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,iistrtor)
               v2 = kst(2,iistrtor)
               v3 = kst(3,iistrtor)
#ifdef USE_NVSHMEM_CUDA
               c1 = d_tors1(ipe)%pel(3,ind)
               s1 = d_tors1(ipe)%pel(4,ind)
               c2 = d_tors2(ipe)%pel(3,ind)
               s2 = d_tors2(ipe)%pel(4,ind)
               c3 = d_tors3(ipe)%pel(3,ind)
               s3 = d_tors3(ipe)%pel(4,ind)
#else
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
#endif
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2   = 2.0_re_p * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3   = cosine*sine2 + sine*cosine2
               phi1    = 1.0_re_p + (cosine*c1 + sine*s1)
               phi2    = 1.0_re_p + (cosine2*c2 + sine2*s2)
               phi3    = 1.0_re_p + (cosine3*c3 + sine3*s3)
               dphi1   = (cosine*s1 - sine*c1)
               dphi2   = 2.0_re_p * (cosine2*s2 - sine2*c2)
               dphi3   = 3.0_re_p * (cosine3*s3 - sine3*c3)
c
c     calculate the bond-stretch for the central bond
c
               k  = ist(2,iistrtor)
#ifdef USE_NVSHMEM_CUDA
               ipe = (k-1)/nbond_pe
               ind = mod((k-1),nbond_pe) +1
               dr  = rcb - d_bl(ipe)%pel(ind)
#else
               dr = rcb - real(bl(k),r_p)
#endif
               !!!  ---  Warning  --- !!!
               ! In single mode , distance between two strech-torsions atoms 
               ! can reach the ideal.
               ! 'dr' can reach zero and by the occasion make the energy
               ! nill and and force between  c<-->b stretch undefined
               ! To make the following test (NaN) work on pgi you need to compile with "-Kieee" flag
               ! if (ddr-ddr.ne.0.0) print*,ddr,dr,e,rcb,bl(k)
               if (abs(dr)<precm_eps) cycle
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

               !increment the stretch-torsion energy and gradient
               ebt = ebt + e

c             if(dedxia-dedxia.ne.0.0)
c    &           print*,'xa',zcb,dedyt , ycb,dedzt
c             if(dedyia-dedyia.ne.0.0)
c    &           print*,'ya',xcb,dedzt , zcb,dedxt
c             if(dedzia-dedzia.ne.0.0)
c    &           print*,'za',ycb,dedxt , xcb,dedyt
c             if(dedxib-dedxib.ne.0.0)
c    &        print*,'xb',yca,dedzt,zca,dedyt,zdc,dedyu,ydc,dedzu,ddrdx
c             if(dedyib-dedyib.ne.0.0)
c    &        print*,'yb',zca,dedxt,xca,dedzt,xdc,dedzu,zdc,dedxu,ddrdy
c             if(dedzib-dedzib.ne.0.0)
c    &        print*,'zb',xca,dedyt,yca,dedxt,ydc,dedxu,xdc,dedyu,ddrdz
c             if(dedxic-dedxic.ne.0.0)
c    & print*,'xc',zba,dedyt,yba,dedzt,ydb,dedzu,zdb,dedyu,ddrdx,xcb,ddr
c             if(dedyic-dedyic.ne.0.0)
c    & print*,'yc',xba,dedzt,zba,dedxt,zdb,dedxu,xdb,dedzu,ddrdy,ycb,ddr
c             if(dedzic-dedzic.ne.0.0)
c    & print*,'zc',yba,dedxt,xba,dedyt,xdb,dedyu,ydb,dedxu,ddrdz,zcb,ddr
c             if(dedxid-dedxid.ne.0.0)
c    &           print*,'xd',zcb,dedyu,ycb,dedzu
c             if(dedyid-dedyid.ne.0.0)
c    &           print*,'yd',xcb,dedzu , zcb,dedxu
c             if(dedzid-dedzid.ne.0.0)
c    &           print*,'zd',ycb,dedxu , xcb,dedyu

               dedi(1) = dedxia
               dedi(2) = dedyia
               dedi(3) = dedzia
!$acc atomic
               debt(1,ialoc) = debt(1,ialoc) + dedi(1)
!$acc atomic
               debt(2,ialoc) = debt(2,ialoc) + dedi(2)
!$acc atomic
               debt(3,ialoc) = debt(3,ialoc) + dedi(3)
c
               dedi(1) = dedxid
               dedi(2) = dedyid
               dedi(3) = dedzid
!$acc atomic
               debt(1,idloc) = debt(1,idloc) + dedi(1)
!$acc atomic
               debt(2,idloc) = debt(2,idloc) + dedi(2)
!$acc atomic
               debt(3,idloc) = debt(3,idloc) + dedi(3)

c
               dedi(1) = dedxib
               dedi(2) = dedyib
               dedi(3) = dedzib
!$acc atomic
               debt(1,ibloc) = debt(1,ibloc) + dedi(1)
!$acc atomic
               debt(2,ibloc) = debt(2,ibloc) + dedi(2)
!$acc atomic
               debt(3,ibloc) = debt(3,ibloc) + dedi(3)
c
               dedi(1) = dedxic
               dedi(2) = dedyic
               dedi(3) = dedzic
!$acc atomic
               debt(1,icloc) = debt(1,icloc) + dedi(1)
!$acc atomic
               debt(2,icloc) = debt(2,icloc) + dedi(2)
!$acc atomic
               debt(3,icloc) = debt(3,icloc) + dedi(3)
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
!$acc&         present(deamdD,debt)
         do i = 1,nloc
            do k = 1,3
               deamdD(k,i) = deamdD(k,i) + debt(k,i)
            end do
         end do
!$acc serial async present(eDaMD,ebt)
         eDaMD = eDaMD + ebt
!$acc end serial
      end if
      end
