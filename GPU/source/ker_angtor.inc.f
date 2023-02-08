#ifndef IJ2_ANGTOR
#define IJ2_ANGTOR
#include "image.f.inc"
#include "convert.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"
      M_subroutine
     &        ker_angtor(iiangtor,i,ia,ib,ic,id,loc,iat,radian
     &           ,atorunit,fgrp,x,y,z,anat,kant,tors1,tors2,tors3
     &           ,use_polymer,use_group,use_virial
     &           ,eat
#ifdef TINKER_CUF
     &           ,deatx,deaty,deatz
#else
     &           ,e,deat
#endif
     &           ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &           ,tver,tfea)
      use tinheader,only: ti_p
      implicit none
      integer  ,intent(in):: iiangtor,i,tver,tfea,loc(*),iat(4,*)
      logical  ,intent(in):: use_polymer,use_group,use_virial
      real(t_p),intent(in):: atorunit,radian,fgrp,x(*),y(*),z(*)
     &         , anat(*),kant(6,*),tors1(4,*),tors2(4,*),tors3(4,*)
      integer  ,intent(inout):: ia,ib,ic,id
#ifdef TINKER_CUF
      ener_rtyp,intent(inout):: eat
      mdyn_rtyp,intent(inout):: deatx(*),deaty(*),deatz(*)
      real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      real(t_p) e
#else
      real(t_p),intent(inout):: e
      real(r_p),intent(inout):: eat,deat(3,*)
     &         , g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
      integer   grd,ene,act,vrl,grp,plm,k
      real(t_p) e1,e2,rcb
      real(t_p) dt,ddt,dot,dedphi,rt2,ru2,rtru,rba2,rcb2,rdc2
      real(t_p) xt,yt,zt,xu,yu,zu,v1,v2,v3,xtu,ytu,ztu
      real(t_p) c1,c2,c3,s1,s2,s3
      real(t_p) sine,cosine,sine2,cosine2,sine3,cosine3
      real(t_p) phi1,phi2,phi3,dphi1,dphi2,dphi3
      real(t_p) angle1,cosang
      real(t_p) terma,termb,termc,termd
      real(t_p) xia,yia,zia,xib,yib,zib,xic,yic,zic
     &         ,xid,yid,zid,xba,yba,zba,xcb,ycb,zcb
     &         ,xdc,ydc,zdc,xca,yca,zca,xdb,ydb,zdb
      real(t_p) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
     &         ,dedxia,dedyia,dedzia,dedxib,dedyib,dedzib
     &         ,dedxic,dedyic,dedzic,dedxid,dedyid,dedzid
      parameter(
     &      grd=__use_grd__,
     &      ene=__use_ene__,
     &      act=__use_act__,
     &      vrl=__use_vir__,
     &      grp=__use_groups__,
     &      plm=__use_polymer__
     &         )
!$acc routine

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
      if (IAND(tfea,plm).NE.0.and.use_polymer) then
         call image_inl (xba,yba,zba)
         call image_inl (xcb,ycb,zcb)
         call image_inl (xdc,ydc,zdc)
      end if
      rba2 = xba*xba + yba*yba + zba*zba
      rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
      rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
      e    = 0
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
         if (use_polymer.and.IAND(tfea,plm).NE.0) then
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
         v1  = kant(1,iiangtor)
         v2  = kant(2,iiangtor)
         v3  = kant(3,iiangtor)
         k   = iat(2,iiangtor)
         dot = xba*xcb + yba*ycb + zba*zcb
         cosang = -dot / sqrt(rba2*rcb2)
         angle1 = radian * acos(cosang)
         dt  = angle1 - anat(k)
         e1  = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
         dedphi = atorunit * dt * (v1*dphi1 + v2*dphi2 + v3*dphi3)
         ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
         if (IAND(tfea,grp).NE.0.and.use_group) then
            e1 = e1 * fgrp
            dedphi = dedphi * fgrp
            ddt = ddt * fgrp
         end if
c
c     compute derivative components for this interaction
c
         IF (IAND(tver,grd).NE.0) THEN
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
     &               + yca*dedzt - zca*dedyt
     &               + zdc*dedyu - ydc*dedzu
         dedyib = terma*(zba*xt-xba*zt) + termc*(xcb*zt-zcb*xt)
     &               + zca*dedxt - xca*dedzt
     &               + xdc*dedzu - zdc*dedxu
         dedzib = terma*(xba*yt-yba*xt) + termc*(ycb*xt-xcb*yt)
     &               + xca*dedyt - yca*dedxt
     &               + ydc*dedxu - xdc*dedyu
         dedxic = termc*(ycb*zt-zcb*yt) + zba*dedyt
     &               - yba*dedzt + ydb*dedzu - zdb*dedyu
         dedyic = termc*(zcb*xt-xcb*zt) + xba*dedzt
     &               - zba*dedxt + zdb*dedxu - xdb*dedzu
         dedzic = termc*(xcb*yt-ycb*xt) + yba*dedxt
     &               - xba*dedyt + xdb*dedyu - ydb*dedxu
         dedxid = zcb*dedyu - ycb*dedzu
         dedyid = xcb*dedzu - zcb*dedxu
         dedzid = ycb*dedxu - xcb*dedyu
         END IF
c
c     get the angle-torsion values for the second angle
c
         v1  = kant(4,iiangtor)
         v2  = kant(5,iiangtor)
         v3  = kant(6,iiangtor)
         k   = iat(3,iiangtor)
         dot = xcb*xdc + ycb*ydc + zcb*zdc
         cosang = acos(-dot / sqrt(rcb2*rdc2))
         angle1 = radian * cosang
         dt  = angle1 - anat(k)
         e2  = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
         dedphi = atorunit * dt * (v1*dphi1 + v2*dphi2 + v3*dphi3)
         ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
         if (IAND(tfea,grp).NE.0.and.use_group) then
            e2 = e2 * fgrp
            dedphi = dedphi * fgrp
            ddt = ddt * fgrp
         end if
c
c     increment the angle-torsion energy
c
         IF (IAND(tver,ene).NE.0) THEN
         e   = e1 + e2
#ifdef TINKER_CUF
         eat = eat + tp2enr(e)
#else
         eat = eat + e
#endif
         END IF
c
c     compute derivative components for this interaction
c
         IF (IAND(tver,grd).NE.0) THEN
         dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
         dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
         dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
         dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
         dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
         dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c     increment chain rule components for the second angle
c
         ia = loc(ia)
         ib = loc(ib)
         ic = loc(ic)
         id = loc(id)
         termb = -ddt / (rcb2*sqrt(ru2))
         termd =  ddt / (rdc2*sqrt(ru2))
         dedxia = dedxia + zcb*dedyt - ycb*dedzt
         dedyia = dedyia + xcb*dedzt - zcb*dedxt
         dedzia = dedzia + ycb*dedxt - xcb*dedyt
         dedxib = dedxib + termb*(zcb*yu-ycb*zu) + yca*dedzt
     &               - zca*dedyt + zdc*dedyu - ydc*dedzu
         dedyib = dedyib + termb*(xcb*zu-zcb*xu) + zca*dedxt
     &               - xca*dedzt + xdc*dedzu - zdc*dedxu
         dedzib = dedzib + termb*(ycb*xu-xcb*yu) + xca*dedyt
     &               - yca*dedxt + ydc*dedxu - xdc*dedyu
         dedxic = dedxic + termb*(ycb*zu-zcb*yu)
     &               + termd*(zdc*yu-ydc*zu) + zba*dedyt
     &               - yba*dedzt + ydb*dedzu - zdb*dedyu
         dedyic = dedyic + termb*(zcb*xu-xcb*zu)
     &               + termd*(xdc*zu-zdc*xu) + xba*dedzt
     &               - zba*dedxt + zdb*dedxu - xdb*dedzu
         dedzic = dedzic + termb*(xcb*yu-ycb*xu)
     &               + termd*(ydc*xu-xdc*yu) + yba*dedxt
     &               - xba*dedyt + xdb*dedyu - ydb*dedxu
         dedxid = dedxid + termd*(ydc*zu-zdc*yu)
     &               + zcb*dedyu - ycb*dedzu
         dedyid = dedyid + termd*(zdc*xu-xdc*zu)
     &               + xcb*dedzu - zcb*dedxu
         dedzid = dedzid + termd*(xdc*yu-ydc*xu)
     &               + ycb*dedxu - xcb*dedyu
c     increment the angle-torsion gradient
#ifdef TINKER_CUF
         call atomic_add_f1( deatx(ia),dedxia )
         call atomic_add_f1( deaty(ia),dedyia )
         call atomic_add_f1( deatz(ia),dedzia )
         call atomic_add_f1( deatx(ib),dedxib )
         call atomic_add_f1( deaty(ib),dedyib )
         call atomic_add_f1( deatz(ib),dedzib )
         call atomic_add_f1( deatx(ic),dedxic )
         call atomic_add_f1( deaty(ic),dedyic )
         call atomic_add_f1( deatz(ic),dedzic )
         call atomic_add_f1( deatx(id),dedxid )
         call atomic_add_f1( deaty(id),dedyid )
         call atomic_add_f1( deatz(id),dedzid )
#else
         call atomic_add( deat(1,ia),dedxia )
         call atomic_add( deat(2,ia),dedyia )
         call atomic_add( deat(3,ia),dedzia )
         call atomic_add( deat(1,ib),dedxib )
         call atomic_add( deat(2,ib),dedyib )
         call atomic_add( deat(3,ib),dedzib )
         call atomic_add( deat(1,ic),dedxic )
         call atomic_add( deat(2,ic),dedyic )
         call atomic_add( deat(3,ic),dedzic )
         call atomic_add( deat(1,id),dedxid )
         call atomic_add( deat(2,id),dedyid )
         call atomic_add( deat(3,id),dedzid )
#endif
         ENDIF
c
c     increment the internal virial tensor components
c
         if (use_virial.and.IAND(tver,vrl).NE.0) then
            g_vxx =g_vxx + xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
            g_vxy =g_vxy + ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
            g_vxz =g_vxz + zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
            g_vyy =g_vyy + ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
            g_vyz =g_vyz + zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
            g_vzz =g_vzz + zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
         end if
      end if
      end subroutine
#endif
