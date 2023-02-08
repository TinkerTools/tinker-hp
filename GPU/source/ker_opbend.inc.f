#include "image.f.inc"
#include "groups.inc.f"
c
c     Comments
c     Zhi Wang, Jun 25, 2019
c
c     The original implementation in Tinker uses ACOS(cosine) to calculate the
c     out-of-plane angle, which is the major source of error in the single
c     precision mode.
c
c     These angles (theta) are usually very small (e.g. 0.001 rad), so the value of
c     variable cosine is very close to 1 (cosine = SQRT(1 - eps**2)). As a result,
c     it is much more accurate to use theta = ASIN(eps) instead of theta =
c     ACOS(cosine) to calculate the angles.
c
#include "atomicOp.inc.f"
        M_subroutine
     &       ker_opbend(ia,ib,ic,id,opbtypInt,loc
     &          ,use_group,use_polymer
     &          ,opbunit,fgrp,force,copb,qopb,popb,sopb,x,y,z
     &          ,eopb
#ifdef TINKER_CUF
     &          ,deopbx,deopby,deopbz
#else
     &          ,e,deopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
     &          ,tver,tfea)
        use angpot    ,only: OPB_W_D_C,OPB_ALLINGER
        use math      ,only: radian
        use tinheader ,only: ti_p
        implicit none
        integer  ,intent(in):: opbtypInt,loc(*),tver,tfea
        integer  ,intent(inout):: ia,ib,ic,id
        logical  ,intent(in):: use_group,use_polymer
        real(t_p),intent(in):: opbunit,copb,qopb,popb,sopb,force,fgrp
        real(r_p),intent(in):: x(*),y(*),z(*)
#ifdef TINKER_CUF
        ener_rtyp,intent(inout)::  eopb
        mdyn_rtyp,intent(inout):: deopb(3,*)
        real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
        real(t_p) e
#else
        real(t_p),intent(  out):: e
        real(r_p),intent(inout):: eopb,deopb(3,*)
     &           ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
        integer   grd,ene,vrl,plm,grp
        real(t_p) angle1,cosine,cosine1,sine,cc,ee,bkk2,term
        real(t_p) dot,deddt,dedcos,dt,dt2,dt3,dt4
        real(r_p) xib,yib,zib ,xid,yid,zid
        real(t_p) xab,yab,zab ,xcb,ycb,zcb ,xdb,ydb,zdb
     &           ,xad,yad,zad ,xcd,ycd,zcd
        real(t_p) rdb2
        real(t_p) rad2,rcd2,rab2,rcb2
        real(t_p) dccdxia,dccdyia,dccdzia ,dccdxic,dccdyic,dccdzic
     &           ,dccdxid,dccdyid,dccdzid ,deedxia,deedyia,deedzia
     &           ,deedxic,deedyic,deedzic ,deedxid,deedyid,deedzid
     &           ,dedxia,dedyia,dedzia ,dedxib,dedyib,dedzib
     &           ,dedxic,dedyic,dedzic ,dedxid,dedyid,dedzid
        parameter(
     &        grd=__use_grd__,ene=__use_ene__,
     &        vrl=__use_vir__,
     &        grp=__use_groups__,plm=__use_polymer__
     &           )
!$acc routine

        xib = x(ib)
        yib = y(ib)
        zib = z(ib)
        xid = x(id)
        yid = y(id)
        zid = z(id)
c
c     compute the out-of-plane bending angle
c
        xab = x(ia) - xib
        yab = y(ia) - yib
        zab = z(ia) - zib
        xcb = x(ic) - xib
        ycb = y(ic) - yib
        zcb = z(ic) - zib
        xdb = xid - xib
        ydb = yid - yib
        zdb = zid - zib
        xad = x(ia) - xid
        yad = y(ia) - yid
        zad = z(ia) - zid
        xcd = x(ic) - xid
        ycd = y(ic) - yid
        zcd = z(ic) - zid
        if (use_polymer.and.IAND(tfea,plm).NE.0) then
           call image_inl(xab,yab,zab)
           call image_inl(xcb,ycb,zcb)
           call image_inl(xdb,ydb,zdb)
           call image_inl(xad,yad,zad)
           call image_inl(xcd,ycd,zcd)
        end if
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
        if (opbtypInt .eq. OPB_W_D_C) then
           rab2 = xab*xab + yab*yab + zab*zab
           rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
           dot  = xab*xcb+yab*ycb+zab*zcb
           cc   = rab2*rcb2 - dot*dot
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
        else if (opbtypInt .eq. OPB_ALLINGER) then
           rad2 = xad*xad + yad*yad + zad*zad
           rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
           dot  = xad*xcd + yad*ycd + zad*zcd
           cc   = rad2*rcd2 - dot*dot
        end if
c
c     find the out-of-plane angle bending energy
c
        ee   = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &       + zdb*(xab*ycb-yab*xcb)
        rdb2 = max(xdb*xdb + ydb*ydb + zdb*zdb, 1e-4)
        if (rdb2.ne.0.0_ti_p .and. cc.ne.0.0_ti_p) then
           !bkk2    = rdb2 - ee*ee/cc
           !cosine1 = sqrt(bkk2/rdb2)
           !cosine  = min(1.0_re_p,max(-1.0_re_p,cosine1))
           !angle1  = radian * acos(cosine)
           sine   = abs(ee) * sqrt(1.0/(cc*rdb2))
           sine   = min(1.0,sine)
           angle1 = radian * asin(sine)
           dt     = angle1
           dt2    = dt* dt
           dt3    = dt2* dt
           dt4    = dt2* dt2
c
c     increment the out-of-plane bending energy
c
           IF (IAND(tver,ene).NE.0) THEN
           e    = opbunit* force* dt2
     &          * (1.0+copb*dt+qopb*dt2+popb*dt3+sopb*dt4)
          !scale the interaction based on its group membership
           if(use_group.and.IAND(tfea,grp).NE.0) e= e*fgrp
#ifdef TINKER_CUF
              eopb = eopb + tp2enr(e)
#else
              eopb = eopb + e
#endif
           END IF
c
           IF (IAND(tver,grd).NE.0) THEN
           deddt = opbunit * force * dt * radian
     &            * (2.0_ti_p + 3.0_ti_p*copb*dt + 4.0_ti_p*qopb*dt2
     &            + 5.0_ti_p*popb*dt3 + 6.0_ti_p*sopb*dt4)
           dedcos = -deddt * sign(1.0_ti_p,ee) / sqrt(cc*rdb2-ee*ee)
          !scale the interaction based on its group membership
           if(use_group.and.IAND(tfea,grp).NE.0) dedcos=dedcos*fgrp
c
c     chain rule terms for first derivative components
c
           if (opbtypInt .eq. OPB_W_D_C) then
              term = ee / cc
              dccdxia = (xab*rcb2-xcb*dot) * term
              dccdyia = (yab*rcb2-ycb*dot) * term
              dccdzia = (zab*rcb2-zcb*dot) * term
              dccdxic = (xcb*rab2-xab*dot) * term
              dccdyic = (ycb*rab2-yab*dot) * term
              dccdzic = (zcb*rab2-zab*dot) * term
              dccdxid = 0.0_ti_p
              dccdyid = 0.0_ti_p
              dccdzid = 0.0_ti_p
           else if (opbtypInt .eq. OPB_ALLINGER) then
              term = ee / cc
              dccdxia = (xad*rcd2-xcd*dot) * term
              dccdyia = (yad*rcd2-ycd*dot) * term
              dccdzia = (zad*rcd2-zcd*dot) * term
              dccdxic = (xcd*rad2-xad*dot) * term
              dccdyic = (ycd*rad2-yad*dot) * term
              dccdzic = (zcd*rad2-zad*dot) * term
              dccdxid = -dccdxia - dccdxic
              dccdyid = -dccdyia - dccdyic
              dccdzid = -dccdzia - dccdzic
           end if
           ia = loc(ia)
           ib = loc(ib)
           ic = loc(ic)
           id = loc(id)
           term    = ee / rdb2
           deedxia = ydb*zcb - zdb*ycb
           deedyia = zdb*xcb - xdb*zcb
           deedzia = xdb*ycb - ydb*xcb
           deedxic = yab*zdb - zab*ydb
           deedyic = zab*xdb - xab*zdb
           deedzic = xab*ydb - yab*xdb
           deedxid = ycb*zab - zcb*yab + xdb*term
           deedyid = zcb*xab - xcb*zab + ydb*term
           deedzid = xcb*yab - ycb*xab + zdb*term
c
c     compute first derivative components for this angle
c
           dedxia = dedcos * (dccdxia+deedxia)
           dedyia = dedcos * (dccdyia+deedyia)
           dedzia = dedcos * (dccdzia+deedzia)
           dedxic = dedcos * (dccdxic+deedxic)
           dedyic = dedcos * (dccdyic+deedyic)
           dedzic = dedcos * (dccdzic+deedzic)
           dedxid = dedcos * (dccdxid+deedxid)
           dedyid = dedcos * (dccdyid+deedyid)
           dedzid = dedcos * (dccdzid+deedzid)
           dedxib =-dedxia - dedxic - dedxid
           dedyib =-dedyia - dedyic - dedyid
           dedzib =-dedzia - dedzic - dedzid
c
c     increment the out-of-plane bending gradient
c
#ifdef TINKER_CUF
           call atomic_add_f1( deopbx(ib),dedxib )
           call atomic_add_f1( deopby(ib),dedyib )
           call atomic_add_f1( deopbz(ib),dedzib )
           call atomic_add_f1( deopbx(ia),dedxia )
           call atomic_add_f1( deopby(ia),dedyia )
           call atomic_add_f1( deopbz(ia),dedzia )
           call atomic_add_f1( deopbx(ic),dedxic )
           call atomic_add_f1( deopby(ic),dedyic )
           call atomic_add_f1( deopbz(ic),dedzic )
           call atomic_add_f1( deopbx(id),dedxid )
           call atomic_add_f1( deopby(id),dedyid )
           call atomic_add_f1( deopbz(id),dedzid )
#else
           call atomic_add( deopb(1,ib),dedxib )
           call atomic_add( deopb(2,ib),dedyib )
           call atomic_add( deopb(3,ib),dedzib )
           call atomic_add( deopb(1,ia),dedxia )
           call atomic_add( deopb(2,ia),dedyia )
           call atomic_add( deopb(3,ia),dedzia )
           call atomic_add( deopb(1,ic),dedxic )
           call atomic_add( deopb(2,ic),dedyic )
           call atomic_add( deopb(3,ic),dedzic )
           call atomic_add( deopb(1,id),dedxid )
           call atomic_add( deopb(2,id),dedyid )
           call atomic_add( deopb(3,id),dedzid )
#endif
           END IF
c
c     increment the internal virial tensor components
c
           IF (IAND(tver,vrl).NE.0) THEN
              g_vxx = g_vxx + xab*dedxia + xcb*dedxic + xdb*dedxid
              g_vxy = g_vxy + yab*dedxia + ycb*dedxic + ydb*dedxid
              g_vxz = g_vxz + zab*dedxia + zcb*dedxic + zdb*dedxid
              g_vyy = g_vyy + yab*dedyia + ycb*dedyic + ydb*dedyid
              g_vyz = g_vyz + zab*dedyia + zcb*dedyic + zdb*dedyid
              g_vzz = g_vzz + zab*dedzia + zcb*dedzic + zdb*dedzid
           END IF
        end if
        end subroutine

