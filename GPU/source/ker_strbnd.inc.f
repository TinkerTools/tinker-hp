#include "image.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"
        M_subroutine
     &       ker_strbnd(i,ia,ib,ic,istrbnd,loc,isb,typeA,aMDwattype
     &          ,use_polymer,use_group,use_amd_wat1
     &          ,stbnunit,fgrp,force1,force2,anat,bl,x,y,z
     &          ,eba,eW1aMD
#ifdef TINKER_CUF
     &          ,debax,debay,debaz,deaMD
#else
     &          ,e,deba,deaMD
#endif
     &          ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &          ,tver,tfea)
        use tinheader,only: ti_p
        use math     ,only: radian
        implicit none
        integer  ,intent(in):: i,istrbnd,loc(*),tver,tfea,isb(3,*)
     &           , typeA(*),aMDwattype(2)
        integer  ,intent(inout):: ia,ib,ic
        logical  ,intent(in   ):: use_group,use_polymer,use_amd_wat1
        real(t_p),intent(in   ):: stbnunit,force1,force2,fgrp
     &           , bl(*),anat(*),x(*),y(*),z(*)
#ifdef TINKER_CUF
        ener_rtyp,intent(inout)::  eopb
        mdyn_rtyp,intent(inout):: deopb(3,*)
        real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
        real(r_p),intent(inout):: eW1aMD,deaMD(3,*)
        real(t_p) e
#else
        real(t_p),intent(  out):: e
        real(r_p),intent(inout):: eba,eW1aMD,deba(3,*),deaMD(3,*)
     &           ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
        integer j,k
        integer grd,ene,vrl,plm,grp,gam
        real(t_p) dr1,dr2,dt,angle1
        real(t_p) dot,cosine
        real(t_p) xia,yia,zia,xib,yib,zib,xic,yic,zic
        real(t_p) xab,yab,zab,xcb,ycb,zcb
        real(t_p) rab,rab2,rcb,rcb2,xp,yp,zp,rp
        real(t_p) term1,term2,termr,term1t,term2t
        real(t_p) ddtdxia,ddtdyia,ddtdzia,ddtdxic,ddtdyic,ddtdzic
        real(t_p) ddrdxia,ddrdyia,ddrdzia,ddrdxic,ddrdyic,ddrdzic
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        parameter(
     &        grd=__use_grd__,ene=__use_ene__,
     &        vrl=__use_vir__,gam=__use_gamd__,
     &        grp=__use_groups__,plm=__use_polymer__
     &           )
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
c
c     compute the value of the bond angle
c
        xab = xia - xib
        yab = yia - yib
        zab = zia - zib
        xcb = xic - xib
        ycb = yic - yib
        zcb = zic - zib
        if (use_polymer.and.IAND(tfea,plm).NE.0) then
           call image_inl (xab,yab,zab)
           call image_inl (xcb,ycb,zcb)
        end if
        rab2 = xab*xab + yab*yab + zab*zab
        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
        if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
           rab = sqrt(rab2)
           rcb = sqrt(rcb2)
           xp  = ycb*zab - zcb*yab
           yp  = zcb*xab - xcb*zab
           zp  = xcb*yab - ycb*xab
           rp  = sqrt(xp*xp + yp*yp + zp*zp)
           rp  = max(rp,0.001_ti_p)
           dot = xab*xcb + yab*ycb + zab*zcb
           cosine = dot / (rab*rcb)
           cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
           angle1 = radian * acos(cosine)
c
c     find chain rule terms for the bond angle deviation
c
           dt      = angle1 - anat(i)
           term1   = -radian / (rab2*rp)
           term2   = radian / (rcb2*rp)
           ddtdxia = term1 * (yab*zp-zab*yp)
           ddtdyia = term1 * (zab*xp-xab*zp)
           ddtdzia = term1 * (xab*yp-yab*xp)
           ddtdxic = term2 * (ycb*zp-zcb*yp)
           ddtdyic = term2 * (zcb*xp-xcb*zp)
           ddtdzic = term2 * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the bond length deviations
c
           j = isb(2,istrbnd)
           k = isb(3,istrbnd)
#ifdef USE_NVSHMEM_CUDA
           ipe = (j-1)/nbond_pe
           ind = mod((j-1),nbond_pe) +1
           dr1 = rab - d_bl(ipe)%pel(ind)
           ipe = (k-1)/nbond_pe
           ind = mod((k-1),nbond_pe) +1
           dr2 = rcb - d_bl(ipe)%pel(ind)
#else
           dr1 = rab - bl(j)
           dr2 = rcb - bl(k)
#endif
           term1   = 1.0_ti_p / rab
           term2   = 1.0_ti_p / rcb
           ddrdxia = term1 * xab
           ddrdyia = term1 * yab
           ddrdzia = term1 * zab
           ddrdxic = term2 * xcb
           ddrdyic = term2 * ycb
           ddrdzic = term2 * zcb
c
c     abbreviations used in defining chain rule terms
c
           term1 = stbnunit * force1
           term2 = stbnunit * force2
           termr = term1*dr1 + term2*dr2
           term1t = term1 * dt
           term2t = term2 * dt
c
c     scale the interaction based on its group membership
c
           if (use_group.and.IAND(tfea,grp).NE.0) then
              termr = termr * fgrp
              term1t = term1t * fgrp
              term2t = term2t * fgrp
           end if
c
c     get the energy and master chain rule terms for derivatives
c
           IF (IAND(tver,ene).NE.0) THEN
          !increment the total stretch-bend energy
              e = termr * dt
#ifdef TINKER_CUF
              eba = eba + tp2enr(e)
#else
              eba = eba + e
#endif
           END IF

           IF (IAND(tver,grd).NE.0) THEN
           ia = loc(ia)
           ib = loc(ib)
           ic = loc(ic)
           dedxia = term1t*ddrdxia + termr*ddtdxia
           dedyia = term1t*ddrdyia + termr*ddtdyia
           dedzia = term1t*ddrdzia + termr*ddtdzia
           dedxic = term2t*ddrdxic + termr*ddtdxic
           dedyic = term2t*ddrdyic + termr*ddtdyic
           dedzic = term2t*ddrdzic + termr*ddtdzic
           dedxib = -dedxia - dedxic
           dedyib = -dedyia - dedyic
           dedzib = -dedzia - dedzic
c
c     Increment the total stretch-bend derivatives
c
#ifdef TINKER_CUF
           call atomic_add( debax(ib),dedxib )
           call atomic_add( debay(ib),dedyib )
           call atomic_add( debaz(ib),dedzib )
           call atomic_add( debax(ia),dedxia )
           call atomic_add( debay(ia),dedyia )
           call atomic_add( debaz(ia),dedzia )
           call atomic_add( debax(ic),dedxic )
           call atomic_add( debay(ic),dedyic )
           call atomic_add( debaz(ic),dedzic )
#else
           call atomic_add( deba(1,ib),dedxib )
           call atomic_add( deba(2,ib),dedyib )
           call atomic_add( deba(3,ib),dedzib )
           call atomic_add( deba(1,ia),dedxia )
           call atomic_add( deba(2,ia),dedyia )
           call atomic_add( deba(3,ia),dedzia )
           call atomic_add( deba(1,ic),dedxic )
           call atomic_add( deba(2,ic),dedyic )
           call atomic_add( deba(3,ic),dedzic )
#endif
           END IF
c
c     aMD storage if waters are considered
c
           if (use_amd_wat1.and.IAND(tfea,gam).NE.0) then
           if (typeA(ia) == aMDwattype(1) .or. typeA(ib)
     $     == aMDwattype(1) .or. typeA(ic) == aMDwattype(1)) then
              eW1aMD = eW1aMD + e
              call atomic_add( deaMD(1,ia),dedxia )
              call atomic_add( deaMD(2,ia),dedyia )
              call atomic_add( deaMD(3,ia),dedzia )
              call atomic_add( deaMD(1,ib),dedxib )
              call atomic_add( deaMD(2,ib),dedyib )
              call atomic_add( deaMD(3,ib),dedzib )
              call atomic_add( deaMD(1,ic),dedxic )
              call atomic_add( deaMD(2,ic),dedyic )
              call atomic_add( deaMD(3,ic),dedzic )
           end if
           end if
c
c     increment the internal virial tensor components
c
           IF (IAND(tver,vrl).NE.0) THEN
              g_vxx = g_vxx + xab*dedxia + xcb*dedxic
              g_vxy = g_vxy + yab*dedxia + ycb*dedxic
              g_vxz = g_vxz + zab*dedxia + zcb*dedxic
              g_vyy = g_vyy + yab*dedyia + ycb*dedyic
              g_vyz = g_vyz + zab*dedyia + zcb*dedyic
              g_vzz = g_vzz + zab*dedzia + zcb*dedzic
           END IF
        end if
        end subroutine

