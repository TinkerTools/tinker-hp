#ifndef KER_STRTOR_INC_F
#define KER_STRTOR_INC_F

#include "image.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"

        M_subroutine
     &          ker_strtor(iistrtor,i,ia,ib,ic,id,ist,loc,x,y,z
     &             ,tors1,tors2,tors3,kst,bl,storunit
     &             ,nbloc,n,use_group,use_polymer
     &             ,e,ebt,debt,vxx,vxy,vxz,vyy,vyz,vzz,ver,fea
     &             )
        use tinheader ,only: precm_eps
        use ktrtor    ,only: maxntt
!$acc routine seq
        implicit none
        integer  ,intent(in):: iistrtor,i,nbloc,n,ver,fea
     &           ,ist(4,*),loc(*)
        logical  ,intent(in):: use_group,use_polymer
        real(t_p),intent(in):: storunit,kst(9,*),bl(*)
     &           ,tors1(4,*),tors2(4,*),tors3(4,*)
        real(r_p),intent(in):: x(*),y(*),z(*)
        integer   ia,ib,ic,id
        real(t_p) e
        real(r_p) ebt,debt(3,*)
     &           ,vxx,vxy,vxz,vyy,vyz,vzz

        integer k
        integer grd,ene,act,vir,gamd,grp,plm
#ifdef   USE_NVSHMEM_CUDA
        integer ipe,ind
#endif  
        real(r_p) dr,ddr,dedphi
        real(r_p) rt2,ru2,rtru
        real(r_p) rba,rcb,rdc
        real(t_p) e1,e2,e3
        real(t_p) xt,yt,zt
        real(t_p) xu,yu,zu
        real(t_p) xtu,ytu,ztu
        real(t_p) v1,v2,v3
        real(t_p) c1,c2,c3
        real(t_p) s1,s2,s3
        real(r_p) sine,cosine
        real(r_p) sine2,cosine2
        real(r_p) sine3,cosine3
        real(r_p) phi1,phi2,phi3
        real(r_p) dphi1,dphi2,dphi3
        real(r_p) xia,yia,zia
        real(r_p) xib,yib,zib
        real(r_p) xic,yic,zic
        real(r_p) xid,yid,zid
        real(t_p) xba,yba,zba
        real(t_p) xcb,ycb,zcb
        real(t_p) xdc,ydc,zdc
        real(t_p) xca,yca,zca
        real(t_p) xdb,ydb,zdb
        real(r_p) ddrdx,ddrdy,ddrdz
        real(r_p) dedxt,dedyt,dedzt
        real(r_p) dedxu,dedyu,dedzu
        real(r_p) dedxia,dedyia,dedzia
        real(r_p) dedxib,dedyib,dedzib
        real(r_p) dedxic,dedyic,dedzic
        real(r_p) dedxid,dedyid,dedzid
        real(t_p) fgrp
        integer   iga,igb,igc,igd,gmin,gmax
        real(t_p),parameter::rtiny=0.000001_ti_p
        logical proceed
        parameter(
     &           grd=__use_grd__,
     &           ene=__use_ene__,
     &           act=__use_act__,
     &           vir=__use_vir__,
     &          gamd=__use_gamd__,
     &           grp=__use_groups__,
     &           plm=__use_polymer__
     &           )

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
        if (IAND(fea,plm).NE.0.and.use_polymer) then
           call image_inl (xba,yba,zba)
           call image_inl (xcb,ycb,zcb)
           call image_inl (xdc,ydc,zdc)
        end if
        rba = sqrt(xba*xba + yba*yba + zba*zba)
        rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
        rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
        if (min(rba,rcb,rdc) .ne. 0.0_ti_p) then
           xt  = yba*zcb - ycb*zba
           yt  = zba*xcb - zcb*xba
           zt  = xba*ycb - xcb*yba
           xu  = ycb*zdc - ydc*zcb
           yu  = zcb*xdc - zdc*xcb
           zu  = xcb*ydc - xdc*ycb
           xtu = yt*zu - yu*zt
           ytu = zt*xu - zu*xt
           ztu = xt*yu - xu*yt
           rt2 = max( xt*xt + yt*yt + zt*zt,rtiny )
           ru2 = max( xu*xu + yu*yu + zu*zu,rtiny )
           rtru= sqrt(rt2 * ru2)
           rt2 = 1.0_re_p/(rt2*rcb)
           ru2 = 1.0_re_p/(ru2*rcb)
c
c     chain rule terms for first derivative components
c
           xca = xic - xia
           yca = yic - yia
           zca = zic - zia
           xdb = xid - xib
           ydb = yid - yib
           zdb = zid - zib
           if (IAND(fea,plm).NE.0.and.use_polymer) then
              call image_inl (xca,yca,zca)
              call image_inl (xdb,ydb,zdb)
           end if
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
           cosine  = (xt*xu + yt*yu + zt*zu) / rtru
           sine    = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
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

           !zeros
           e1=0;e2=0;e3=0
           dedxia=0;dedyia=0;dedzia=0
           dedxib=0;dedyib=0;dedzib=0
           dedxic=0;dedyic=0;dedzic=0
           dedxid=0;dedyid=0;dedzid=0
c
c     calculate the bond-stretch for the first bond
c
           k   = ist(2,iistrtor)
           v1  = kst(1,iistrtor)
           v2  = kst(2,iistrtor)
           v3  = kst(3,iistrtor)
#ifdef USE_NVSHMEM_CUDA
           ipe = (k-1)/nbond_pe
           ind = mod((k-1),nbond_pe) +1
           dr  = rba - d_bl(ipe)%pel(ind)
#else
           dr  = rba - real(bl(k),r_p)
#endif
           !!!  ---  Warning  --- !!!
           ! In single mode , distance between two strech-torsions atoms 
           ! can reach the ideal.
           ! 'dr' can reach zero and by the occasion make the energy
           ! nill and and force between  c<-->b stretch undefined
           ! To make the following test (NaN) work on pgi you need to compile with "-Kieee" flag
           ! if (ddr-ddr.ne.0.0) print*,ddr,dr,e,rcb,bl(k)
           if (abs(dr)>precm_eps) then
           e1     = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
           dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
           ddr    = storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rba
c
c     scale the interaction based on its group membership
c
           if (IAND(fea,grp).NE.0.and.use_group) then
              e1 = e1 * fgrp
              dedphi = dedphi * fgrp
              ddr = ddr * fgrp
           end if
           IF (IAND(ver,grd).NE.0) THEN
c
c     compute derivative components for this interaction
c
           ddrdx = xba * ddr
           ddrdy = yba * ddr
           ddrdz = zba * ddr
           dedxt =  dedphi * (yt*zcb - ycb*zt) * (rt2)
           dedyt =  dedphi * (zt*xcb - zcb*xt) * (rt2)
           dedzt =  dedphi * (xt*ycb - xcb*yt) * (rt2)
           dedxu = -dedphi * (yu*zcb - ycb*zu) * (ru2)
           dedyu = -dedphi * (zu*xcb - zcb*xu) * (ru2)
           dedzu = -dedphi * (xu*ycb - xcb*yu) * (ru2)
c
c     determine chain rule components for the first bond
c
           dedxia = zcb*dedyt - ycb*dedzt - ddrdx
           dedyia = xcb*dedzt - zcb*dedxt - ddrdy
           dedzia = ycb*dedxt - xcb*dedyt - ddrdz
           dedxib = yca*dedzt - zca*dedyt + zdc*dedyu
     &            - ydc*dedzu + ddrdx
           dedyib = zca*dedxt - xca*dedzt + xdc*dedzu
     &            - zdc*dedxu + ddrdy
           dedzib = xca*dedyt - yca*dedxt + ydc*dedxu
     &            - xdc*dedyu + ddrdz
           dedxic = zba*dedyt - yba*dedzt + ydb*dedzu
     &            - zdb*dedyu
           dedyic = xba*dedzt - zba*dedxt + zdb*dedxu
     &            - xdb*dedzu
           dedzic = yba*dedxt - xba*dedyt + xdb*dedyu
     &                 - ydb*dedxu
           dedxid = zcb*dedyu - ycb*dedzu
           dedyid = xcb*dedzu - zcb*dedxu
           dedzid = ycb*dedxu - xcb*dedyu
           END IF
           end if
c
c     get the stretch-torsion values for the second bond
c
           k  = ist(3,iistrtor)
           v1 = kst(4,iistrtor)
           v2 = kst(5,iistrtor)
           v3 = kst(6,iistrtor)
#ifdef USE_NVSHMEM_CUDA
           ipe = (k-1)/nbond_pe
           ind = mod((k-1),nbond_pe) +1
           dr  = rcb - d_bl(ipe)%pel(ind)
#else
           dr = rcb - real(bl(k),r_p)
#endif
           if (abs(dr)>precm_eps) then
           e2     = storunit * dr * (v1* phi1 + v2* phi2 + v3* phi3)
           dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
           ddr    = storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rcb
c
c     scale the interaction based on its group membership
c
           if (IAND(fea,grp).NE.0.and.use_group) then
              e2  = e2 * fgrp
              dedphi = dedphi * fgrp
              ddr = ddr * fgrp
           end if
c
c     compute derivative components for this interaction
c
           IF (IAND(ver,grd).NE.0) THEN
           ddrdx = xcb * ddr
           ddrdy = ycb * ddr
           ddrdz = zcb * ddr
           dedxt =  dedphi * (yt*zcb - ycb*zt) * (rt2)
           dedyt =  dedphi * (zt*xcb - zcb*xt) * (rt2)
           dedzt =  dedphi * (xt*ycb - xcb*yt) * (rt2)
           dedxu = -dedphi * (yu*zcb - ycb*zu) * (ru2)
           dedyu = -dedphi * (zu*xcb - zcb*xu) * (ru2)
           dedzu = -dedphi * (xu*ycb - xcb*yu) * (ru2)
c
c     increment chain rule components for the second bond
c
           dedxia = dedxia + zcb*dedyt - ycb*dedzt
           dedyia = dedyia + xcb*dedzt - zcb*dedxt
           dedzia = dedzia + ycb*dedxt - xcb*dedyt
           dedxib = dedxib + yca*dedzt - zca*dedyt + zdc*dedyu
     &                     - ydc*dedzu - ddrdx
           dedyib = dedyib + zca*dedxt - xca*dedzt + xdc*dedzu
     &                     - zdc*dedxu - ddrdy
           dedzib = dedzib + xca*dedyt - yca*dedxt + ydc*dedxu
     &                     - xdc*dedyu - ddrdz
           dedxic = dedxic + zba*dedyt - yba*dedzt + ydb*dedzu
     &                     - zdb*dedyu + ddrdx
           dedyic = dedyic + xba*dedzt - zba*dedxt + zdb*dedxu
     &                     - xdb*dedzu + ddrdy
           dedzic = dedzic + yba*dedxt - xba*dedyt + xdb*dedyu
     &                     - ydb*dedxu + ddrdz
           dedxid = dedxid + zcb*dedyu - ycb*dedzu
           dedyid = dedyid + xcb*dedzu - zcb*dedxu
           dedzid = dedzid + ycb*dedxu - xcb*dedyu
           END IF
           end if

c          if(dedxia-dedxia.ne.0.0)
c    &        print*,'xa',zcb,dedyt , ycb,dedzt
c          if(dedyia-dedyia.ne.0.0)
c    &        print*,'ya',xcb,dedzt , zcb,dedxt
c          if(dedzia-dedzia.ne.0.0)
c    &        print*,'za',ycb,dedxt , xcb,dedyt
c          if(dedxib-dedxib.ne.0.0)
c    &     print*,'xb',yca,dedzt,zca,dedyt,zdc,dedyu,ydc,dedzu,ddrdx
c          if(dedyib-dedyib.ne.0.0)
c    &     print*,'yb',zca,dedxt,xca,dedzt,xdc,dedzu,zdc,dedxu,ddrdy
c          if(dedzib-dedzib.ne.0.0)
c    &     print*,'zb',xca,dedyt,yca,dedxt,ydc,dedxu,xdc,dedyu,ddrdz
c          if(dedxic-dedxic.ne.0.0)
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
c
c     get the stretch-torsion values for the third bond
c
           k  = ist(4,iistrtor)
           v1 = kst(7,iistrtor)
           v2 = kst(8,iistrtor)
           v3 = kst(9,iistrtor)
#ifdef USE_NVSHMEM_CUDA
           ipe= (k-1)/nbond_pe
           ind= mod((k-1),nbond_pe) +1
           dr = rdc - d_bl(ipe)%pel(ind)
#else
           dr = rdc - real(bl(k),r_p)
#endif
           if (abs(dr)>precm_eps) then
           e3     = storunit * dr * (v1* phi1 + v2* phi2 + v3* phi3)
           dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
           ddr    = storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rdc
c
c     scale the interaction based on its group membership
c
           if (IAND(fea,grp).NE.0.and.use_group) then
              e3 = e3 * fgrp
              dedphi = dedphi * fgrp
              ddr = ddr * fgrp
           end if
c
c     compute derivative components for this interaction
c
           IF (IAND(ver,grd).NE.0) THEN
           ddrdx = xdc * ddr
           ddrdy = ydc * ddr
           ddrdz = zdc * ddr
           dedxt =  dedphi * (yt*zcb - ycb*zt) * (rt2)
           dedyt =  dedphi * (zt*xcb - zcb*xt) * (rt2)
           dedzt =  dedphi * (xt*ycb - xcb*yt) * (rt2)
           dedxu = -dedphi * (yu*zcb - ycb*zu) * (ru2)
           dedyu = -dedphi * (zu*xcb - zcb*xu) * (ru2)
           dedzu = -dedphi * (xu*ycb - xcb*yu) * (ru2)
c
c     increment chain rule components for the third bond
c
           dedxia = dedxia + zcb*dedyt - ycb*dedzt
           dedyia = dedyia + xcb*dedzt - zcb*dedxt
           dedzia = dedzia + ycb*dedxt - xcb*dedyt
           dedxib = dedxib + yca*dedzt - zca*dedyt + zdc*dedyu
     &                     - ydc*dedzu
           dedyib = dedyib + zca*dedxt - xca*dedzt + xdc*dedzu
     &                     - zdc*dedxu
           dedzib = dedzib + xca*dedyt - yca*dedxt + ydc*dedxu
     &                     - xdc*dedyu
           dedxic = dedxic + zba*dedyt - yba*dedzt + ydb*dedzu
     &                     - zdb*dedyu - ddrdx
           dedyic = dedyic + xba*dedzt - zba*dedxt + zdb*dedxu
     &                     - xdb*dedzu - ddrdy
           dedzic = dedzic + yba*dedxt - xba*dedyt + xdb*dedyu
     &                     - ydb*dedxu - ddrdz
           dedxid = dedxid + zcb*dedyu - ycb*dedzu + ddrdx
           dedyid = dedyid + xcb*dedzu - zcb*dedxu + ddrdy
           dedzid = dedzid + ycb*dedxu - xcb*dedyu + ddrdz
           END IF
           end if

           e   = e1 + e2 + e3

           !increment the stretch-torsion energy and gradient
           ebt = ebt + e

           IF (IAND(ver,grd).NE.0) THEN
           ia  = loc(ia)
           ib  = loc(ib)
           ic  = loc(ic)
           id  = loc(id)
           call atomic_add( debt(1,ia),dedxia )
           call atomic_add( debt(2,ia),dedyia )
           call atomic_add( debt(3,ia),dedzia )
c
           call atomic_add( debt(1,id),dedxid )
           call atomic_add( debt(2,id),dedyid )
           call atomic_add( debt(3,id),dedzid )
c
           call atomic_add( debt(1,ib),dedxib )
           call atomic_add( debt(2,ib),dedyib )
           call atomic_add( debt(3,ib),dedzib )
c
           call atomic_add( debt(1,ic),dedxic )
           call atomic_add( debt(2,ic),dedyic )
           call atomic_add( debt(3,ic),dedzic )
c
c     increment the internal virial tensor components
c
              IF (IAND(ver,vir).NE.0) THEN
        vxx = vxx + xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
        vxy = vxy + ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
        vxz = vxz + zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
        vyy = vyy + ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
        vyz = vyz + zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
        vzz = vzz + zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
              END IF
           END IF
        end if
        end subroutine
#endif
