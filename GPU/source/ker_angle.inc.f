#include "image.f.inc"
#include "groups.inc.f"
#include "convert.f.inc"
#include "atomicOp.inc.f"
c
c     compute the bond angle bending energy and gradient
c
        M_subroutine
     &          ker_angle(i,ia,ib,ic,loc,ideal,force
     &               ,angunit,cang,pang,sang,qang,fgrp
     &               ,use_group,use_polymer,angtypI
     &               ,x,y,z,afld
     &               ,ea
#if TINKER_CUF       
     &               ,deax,deay,deaz
#else                
     &               ,e,dea
#endif               
     &               ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,tver,tfea
     &               ,fmat_ps,dfmat_ps)
        use angpot   ,only: ANG_HARMONIC,ANG_LINEAR,ANG_FOURIER
     &               , ANG_IN_PLANE, ANG_PS, c5z_ps, idx_ps
        use math     ,only: radian
        use tinheader,only: ti_p
        implicit none
        integer  ,intent(in):: i,loc(*)
     &           , angtypI,tver,tfea
        logical  ,intent(in):: use_group,use_polymer
        real(t_p),intent(in):: ideal,force,angunit
     &           ,cang,pang,sang,qang,fgrp
        real(t_p),intent(in):: x(*),y(*),z(*),afld(*)
        integer  ,intent(inout):: ia,ib,ic
#ifdef TINKER_CUF
        ener_rtyp,intent(inout)::  ea
        mdyn_rtyp,intent(inout):: deax(*),deay(*),deaz(*)
        real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
        real(t_p) e
#else
        real(t_p),intent(  out):: e
        real(r_p),intent(inout):: dea(3,*),ea
     &           , g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
        real(t_p),intent(inout):: fmat_ps(15,3,*),dfmat_ps(15,3,*)
        integer   grd,ene,act,vir,gamd,grp,plm,j
        real(t_p) xia,yia,zia,xib,yib,zib,xic,yic,zic,xab,yab,zab
     &           ,xcb,ycb,zcb,rab2,rcb2,rab,rcb,xp,yp,zp,rp,dot,cosine
     &           ,angle1,r_e,x1,x2,x3,drab,drcb,gaussterm
     &           ,dgaussdrab,dgaussdrcb,term,dtermdrab,dtermdrcb
     &           ,dt,dt2,dt3,dt4,factor,sine,deddt,fold,dedrab,dedrcb
     &           ,terma,termc,dedxia,dedyia,dedzia
     &           ,dedxib,dedyib,dedzib,dedxic,dedyic,dedzic
        parameter(
     &           grd=__use_grd__,
     &           ene=__use_ene__,
     &           act=__use_act__,
     &           vir=__use_vir__,
     &          gamd=__use_gamd__,
     &           grp=__use_groups__,
     &           plm=__use_polymer__
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
        e    = 0
        if (rab2.ne.0.0 .and. rcb2.ne.0.0) then
           xp     = ycb*zab - zcb*yab
           yp     = zcb*xab - xcb*zab
           zp     = xcb*yab - ycb*xab
           rp     = sqrt(xp*xp + yp*yp + zp*zp)
           rp     = max(rp,0.000001)
           dot    = xab*xcb + yab*ycb + zab*zcb
           rab = sqrt(rab2)
           rcb = sqrt(rcb2)
           cosine = dot / (rab*rcb)
           cosine = min(1.0,max(-1.0,cosine))
           angle1 = radian * acos(cosine)
c
c     get the energy and master chain rule term for derivatives
c
           if (angtypI .eq. ANG_HARMONIC) then
              dt  = angle1 - ideal
              dt2 = dt * dt
              dt3 = dt2 * dt
              dt4 = dt2 * dt2
              e   = angunit * force * dt2
     &             * (1.0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
              deddt = angunit * force * dt * radian
     &          * (2.0 + 3.0_ti_p*cang*dt + 4.0_ti_p*qang*dt2
     &                    + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
              dedrab = 0.0_ti_p
              dedrcb = 0.0_ti_p
           else if (angtypI .eq. ANG_LINEAR) then
              factor = 2.0_ti_p * angunit * radian**2
              sine   = sqrt(1.0_ti_p-cosine*cosine)
              e      = factor * force * (1.0_ti_p+cosine)
              deddt  = -factor * force * sine
              dedrab = 0.0_ti_p
              dedrcb = 0.0_ti_p
           else if (angtypI .eq. ANG_FOURIER) then
              fold   = afld(i)
              factor = 2.0_ti_p * angunit * (radian/fold)**2
              cosine = cos((fold*angle1-ideal)/radian)
              sine   = sin((fold*angle1-ideal)/radian)
              e      = factor * force * (1.0_ti_p+cosine)
              deddt  = -factor * force * fold * sine
              dedrab = 0.0_ti_p
              dedrcb = 0.0_ti_p
           elseif (angtypI .eq. ANG_PS) then
              r_e = afld(i)
              x3 = cosine - cos(ideal/radian)
              sine = sin(angle1/radian)
              drab = rab - r_e
              drcb = rcb - r_e 
              gaussterm = exp(-force*(drab**2+drcb**2))
              dgaussdrab = -2.0_ti_p*force*drab*gaussterm
              dgaussdrcb = -2.0_ti_p*force*drcb*gaussterm
              x1 = drab/r_e
              x2 = drcb/r_e
              fmat_ps(1,1,i)=1.0_ti_p
              fmat_ps(1,2,i)=1.0_ti_p
              fmat_ps(1,3,i)=1.0_ti_p
              dfmat_ps(1,1,i)=0.0_ti_p
              dfmat_ps(1,2,i)=0.0_ti_p
              dfmat_ps(1,3,i)=0.0_ti_p
!$acc loop seq
              do j=2,15
                fmat_ps(j,1,i)=fmat_ps(j-1,1,i)*x1
                fmat_ps(j,2,i)=fmat_ps(j-1,2,i)*x2
                fmat_ps(j,3,i)=fmat_ps(j-1,3,i)*x3
                !deriv with respect to rab
                dfmat_ps(j,1,i)=dfmat_ps(j-1,1,i)*x1
     &               + fmat_ps(j-1,1,i)/r_e
                !deriv with respect to rcb
                dfmat_ps(j,2,i)=dfmat_ps(j-1,2,i)*x2 
     &              + fmat_ps(j-1,2,i)/r_e
                !deriv with respect to angle
                dfmat_ps(j,3,i)=dfmat_ps(j-1,3,i)*x3 
     &              - fmat_ps(j-1,3,i)*sine
              enddo
              e=0.0_ti_p
              deddt=0._ti_p
              dedrab=0._ti_p
              dedrcb=0._ti_p
!$acc loop seq
              do j=2,245
                term=fmat_ps(idx_ps(j,1),1,i)
     &                    *fmat_ps(idx_ps(j,2),2,i) 
     &               +fmat_ps(idx_ps(j,2),1,i)
     &                    *fmat_ps(idx_ps(j,1),2,i)
                dtermdrab=
     &              dfmat_ps(idx_ps(j,1),1,i)
     &                     *fmat_ps(idx_ps(j,2),2,i) 
     &             +dfmat_ps(idx_ps(j,2),1,i)
     &                     *fmat_ps(idx_ps(j,1),2,i)
                dtermdrcb=
     &              fmat_ps(idx_ps(j,1),1,i)
     &                     *dfmat_ps(idx_ps(j,2),2,i) 
     &             +fmat_ps(idx_ps(j,2),1,i)
     &                     *dfmat_ps(idx_ps(j,1),2,i)

                e=e+c5z_ps(j)*term*fmat_ps(idx_ps(j,3),3,i)
                deddt=deddt+c5z_ps(j)*term
     &              *dfmat_ps(idx_ps(j,3),3,i)
                dedrab=dedrab+c5z_ps(j)
     &              *dtermdrab*fmat_ps(idx_ps(j,3),3,i)
                dedrcb=dedrcb+c5z_ps(j)
     &              *dtermdrcb*fmat_ps(idx_ps(j,3),3,i)
              enddo
              deddt = deddt*gaussterm
              dedrab = dedrab*gaussterm + e*dgaussdrab
              dedrcb = dedrcb*gaussterm + e*dgaussdrcb
              e =  e*gaussterm ! + c5z_ps(1)
           end if
c
c     scale the interaction based on its group membership
c
           if (use_group.and.IAND(tfea,grp).NE.0) then
              e     =  e    * fgrp
              deddt = deddt * fgrp
              if (dedrab .ne. 0.0_ti_p) dedrab = dedrab * fgrp
              if (dedrcb .ne. 0.0_ti_p) dedrcb = dedrcb * fgrp
           end if
c
c     increment the total bond angle energy and derivatives
c
#ifdef TINKER_CUF
           IF (IAND(tver,ene).ne.0) ea = ea + tp2enr(e)
#else
           IF (IAND(tver,ene).ne.0) ea = ea + e
#endif

           IF (IAND(tver,grd).NE.0) THEN
c
c     compute derivative components for this interaction
c
              ia = loc(ia)
              ib = loc(ib)
              ic = loc(ic)
              terma  = -deddt / (rab2*rp)
              termc  =  deddt / (rcb2*rp)
              dedxia = terma * (yab*zp-zab*yp)
              dedyia = terma * (zab*xp-xab*zp)
              dedzia = terma * (xab*yp-yab*xp)
              dedxic = termc * (ycb*zp-zcb*yp)
              dedyic = termc * (zcb*xp-xcb*zp)
              dedzic = termc * (xcb*yp-ycb*xp)
              dedxib = -dedxia - dedxic
              dedyib = -dedyia - dedyic
              dedzib = -dedzia - dedzic
              if (dedrab .ne. 0.0_ti_p) then
                term = dedrab/rab
                dedxia = dedxia + term * xab
                dedxib = dedxib - term * xab
                dedyia = dedyia + term * yab
                dedyib = dedyib - term * yab
                dedzia = dedzia + term * zab
                dedzib = dedzib - term * zab
              end if

              if (dedrcb .ne. 0.0_ti_p) then
                term = dedrcb/rcb
                dedxic = dedxic + term * xcb
                dedxib = dedxib - term * xcb
                dedyic = dedyic + term * ycb
                dedyib = dedyib - term * ycb
                dedzic = dedzic + term * zcb
                dedzib = dedzib - term * zcb
              end if
#ifdef TINKER_CUF
              call atomic_add_f1 ( deax(ib),dedxib )
              call atomic_add_f1 ( deay(ib),dedyib )
              call atomic_add_f1 ( deaz(ib),dedzib )
              call atomic_add_f1 ( deax(ia),dedxia )
              call atomic_add_f1 ( deay(ia),dedyia )
              call atomic_add_f1 ( deaz(ia),dedzia )
              call atomic_add_f1 ( deax(ic),dedxic )
              call atomic_add_f1 ( deay(ic),dedyic )
              call atomic_add_f1 ( deaz(ic),dedzic )
#else
              call atomic_add ( dea(1,ib),dedxib )
              call atomic_add ( dea(2,ib),dedyib )
              call atomic_add ( dea(3,ib),dedzib )
              call atomic_add ( dea(1,ia),dedxia )
              call atomic_add ( dea(2,ia),dedyia )
              call atomic_add ( dea(3,ia),dedzia )
              call atomic_add ( dea(1,ic),dedxic )
              call atomic_add ( dea(2,ic),dedyic )
              call atomic_add ( dea(3,ic),dedzic )
#endif
           END IF
c
c     increment the internal virial tensor components
c
           IF (IAND(tver,vir).NE.0) THEN
              g_vxx = g_vxx + xab*dedxia + xcb*dedxic
              g_vxy = g_vxy + yab*dedxia + ycb*dedxic
              g_vxz = g_vxz + zab*dedxia + zcb*dedxic
              g_vyy = g_vyy + yab*dedyia + ycb*dedyic
              g_vyz = g_vyz + zab*dedyia + zcb*dedyic
              g_vzz = g_vzz + zab*dedzia + zcb*dedzic
           END IF
        end if
        end subroutine

        M_subroutine
     &          ker_angle_plan(i,ia,ib,ic,id,loc
     &               ,ideal,force,angunit,cang,pang,sang,qang,fgrp
     &               ,use_group,use_polymer,angtypI
     &               ,x,y,z,afld
     &               ,use_amd_wat1,typeAtm,aMDwattype,eW1aMD
     &               ,ea
#if TINKER_CUF       
     &               ,deax,deay,deaz,deW1aMD
#else                
     &               ,e,dea,deW1aMD
#endif               
     &               ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,tver,tfea)
        use angpot   ,only: ANG_HARMONIC,ANG_LINEAR,ANG_FOURIER
     &               , ANG_IN_PLANE
        use math     ,only: radian
        use tinheader,only: ti_p
        implicit none
        integer  ,intent(in):: i,loc(*)
     &           ,tver,tfea,angtypI,typeAtm(*),aMDwattype(*)
        logical  ,intent(in):: use_group,use_amd_wat1,use_polymer
        real(t_p),intent(in):: ideal,force,angunit,fgrp
     &           ,cang,pang,sang,qang
        real(t_p),intent(in):: x(*),y(*),z(*),afld(*)
        real(r_p),intent(inout):: eW1aMD,deW1aMD(3,*)
        integer  ,intent(inout):: ia,ib,ic,id
#ifdef TINKER_CUF
        ener_rtyp,intent(inout)::  ea
        mdyn_rtyp,intent(inout):: deax(*),deay(*),deaz(*)
        real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
        real(t_p) e
#else
        real(t_p),intent(  out):: e
        real(r_p),intent(inout):: dea(3,*),ea
     &           , g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
        integer   :: grd,ene,vir,act,gamd,grp,plm
        real(t_p) :: xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid
     &           ,xab,yab,zab,xcb,ycb,zcb,xad,yad,zad,xbd,ybd,zbd
     &           ,xcd,ycd,zcd,rab2,rcb2,xp,yp,zp,rp,dot,cosine,angle1
     &           ,xip,yip,zip,xap,yap,zap,xcp,ycp,zcp,rap2,rcp2
     &           ,dt,dt2,dt3,dt4,factor,sine,deddt,fold
     &           ,term,terma,termc,dedxia,dedyia,dedzia
     &           ,dpdxia,dpdyia,dpdzia,dedxib,dedyib,dedzib
     &           ,dedxic,dedyic,dedzic,dedxid,dedyid,dedzid
     &           ,dedxip,dedyip,dedzip,dpdxic,dpdyic,dpdzic
        real(t_p) rt2,ptrt2,delta,delta2
        real(t_p) xt,yt,zt,xm,ym,zm,rm
        parameter(
     &           grd=__use_grd__,
     &           ene=__use_ene__,
     &           act=__use_act__,
     &           vir=__use_vir__,
     &          gamd=__use_gamd__,
     &           grp=__use_groups__,
     &           plm=__use_polymer__
     &           )
!$acc routine

       !compute the projected in-plane angle energy and gradient
        xid = x(id)
        yid = y(id)
        zid = z(id)
        xia = x(ia)
        yia = y(ia)
        zia = z(ia)
        xib = x(ib)
        yib = y(ib)
        zib = z(ib)
        xic = x(ic)
        yic = y(ic)
        zic = z(ic)
        xad = xia - xid
        yad = yia - yid
        zad = zia - zid
        xbd = xib - xid
        ybd = yib - yid
        zbd = zib - zid
        xcd = xic - xid
        ycd = yic - yid
        zcd = zic - zid
        if (use_polymer.and.IAND(tfea,plm).NE.0) then
           call image_inl(xad,yad,zad)
           call image_inl(xbd,ybd,zbd)
           call image_inl(xcd,ycd,zcd)
        end if
        xt  = yad*zcd - zad*ycd
        yt  = zad*xcd - xad*zcd
        zt  = xad*ycd - yad*xcd
        rt2 = xt*xt + yt*yt + zt*zt
        delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
        xip = xib + xt*delta
        yip = yib + yt*delta
        zip = zib + zt*delta
        xap = xia - xip
        yap = yia - yip
        zap = zia - zip
        xcp = xic - xip
        ycp = yic - yip
        zcp = zic - zip
        if (IAND(tfea,plm).NE.0.and.use_polymer) then
           call image_inl (xap,yap,zap)
           call image_inl (xcp,ycp,zcp)
        end if
        rap2 = xap*xap + yap*yap + zap*zap
        rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
        e    = 0
        if (rap2.ne.0.0_ti_p .and. rcp2.ne.0.0_ti_p) then
           xm = ycp*zap - zcp*yap
           ym = zcp*xap - xcp*zap
           zm = xcp*yap - ycp*xap
           rm = sqrt(xm*xm + ym*ym + zm*zm)
           rm = max(rm,0.000001_ti_p)
           dot = xap*xcp + yap*ycp + zap*zcp
           cosine = dot / sqrt(rap2*rcp2)
           cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
           angle1 = radian * acos(cosine)
c
c     get the energy and master chain rule term for derivatives
c
           dt  = angle1 - ideal
           dt2 = dt * dt
           dt3 = dt2 * dt
           dt4 = dt2 * dt2
           e   = angunit * force * dt2
     &           * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
           deddt = angunit * force * dt * radian
     &         * (2.0_ti_p + 3.0_ti_p*cang*dt + 4.0_ti_p*qang*dt2
     &                 + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
c
c     scale the interaction based on its group membership
c
           if (IAND(tfea,grp).NE.0.and.use_group) then
              e     =  e    * fgrp
              deddt = deddt * fgrp
           end if

           ! increment the total bond angle energy and derivatives
#ifdef TINKER_CUF
           IF(IAND(tver,ene).NE.0) ea = ea + tp2enr(e)
#else
           IF(IAND(tver,ene).NE.0) ea = ea + e
#endif

           IF (IAND(tver,grd).NE.0) THEN
c
c     chain rule terms for first derivative components
c
           terma  = -deddt / (rap2*rm)
           termc  =  deddt / (rcp2*rm)
           dedxia = terma * (yap*zm-zap*ym)
           dedyia = terma * (zap*xm-xap*zm)
           dedzia = terma * (xap*ym-yap*xm)
           dedxic = termc * (ycp*zm-zcp*ym)
           dedyic = termc * (zcp*xm-xcp*zm)
           dedzic = termc * (xcp*ym-ycp*xm)
           dedxip = -dedxia - dedxic
           dedyip = -dedyia - dedyic
           dedzip = -dedzia - dedzic
           ia     = loc(ia)
           ib     = loc(ib)
           ic     = loc(ic)
           id     = loc(id)
c
c     chain rule components for the projection of the central atom
c
           delta2 = 2.0_ti_p * delta
           ptrt2  = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
           term   = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
           dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
           term   = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
           dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
           term   = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
           dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
           term   = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
           dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
           term   = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
           dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
           term   = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
           dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c     compute derivative components for this interaction
c
           dedxia =  dedxia + dpdxia
           dedyia =  dedyia + dpdyia
           dedzia =  dedzia + dpdzia
           dedxib =  dedxip
           dedyib =  dedyip
           dedzib =  dedzip
           dedxic =  dedxic + dpdxic
           dedyic =  dedyic + dpdyic
           dedzic =  dedzic + dpdzic
           dedxid = -dedxia - dedxib - dedxic
           dedyid = -dedyia - dedyib - dedyic
           dedzid = -dedzia - dedzib - dedzic
 
#ifdef TINKER_CUF
           call atomic_add_f1 ( deax(ib),dedxib )
           call atomic_add_f1 ( deay(ib),dedyib )
           call atomic_add_f1 ( deaz(ib),dedzib )
           call atomic_add_f1 ( deax(ia),dedxia )
           call atomic_add_f1 ( deay(ia),dedyia )
           call atomic_add_f1 ( deaz(ia),dedzia )
           call atomic_add_f1 ( deax(ic),dedxic )
           call atomic_add_f1 ( deay(ic),dedyic )
           call atomic_add_f1 ( deaz(ic),dedzic )
           call atomic_add_f1 ( deax(id),dedxid )
           call atomic_add_f1 ( deay(id),dedyid )
           call atomic_add_f1 ( deaz(id),dedzid )
#else
           call atomic_add ( dea(1,ib),dedxib )
           call atomic_add ( dea(2,ib),dedyib )
           call atomic_add ( dea(3,ib),dedzib )
           call atomic_add ( dea(1,ia),dedxia )
           call atomic_add ( dea(2,ia),dedyia )
           call atomic_add ( dea(3,ia),dedzia )
           call atomic_add ( dea(1,ic),dedxic )
           call atomic_add ( dea(2,ic),dedyic )
           call atomic_add ( dea(3,ic),dedzic )
           call atomic_add ( dea(1,id),dedxid )
           call atomic_add ( dea(2,id),dedyid )
           call atomic_add ( dea(3,id),dedzid )
#endif
           END IF
c
c     aMD storage if waters are considered
c
           if (use_amd_wat1.and.IAND(tfea,gamd).NE.0) then  ! Check amd
           if (typeAtm(ia) == aMDwattype(1) .or. typeAtm(ib)
     &     == aMDwattype(1) .or. typeAtm(ic) == aMDwattype(1))
     &     then   !Check type
              eW1aMD = eW1aMD + e
              call atomic_add_m ( deW1aMD(1,ib),dedxib )
              call atomic_add_m ( deW1aMD(2,ib),dedyib )
              call atomic_add_m ( deW1aMD(3,ib),dedzib )
              call atomic_add_m ( deW1aMD(1,ia),dedxia )
              call atomic_add_m ( deW1aMD(2,ia),dedyia )
              call atomic_add_m ( deW1aMD(3,ia),dedzia )
              call atomic_add_m ( deW1aMD(1,ic),dedxic )
              call atomic_add_m ( deW1aMD(2,ic),dedyic )
              call atomic_add_m ( deW1aMD(3,ic),dedzic )
           end if  !Check type
           end if  !Check amd
 
           ! Increment the internal virial tensor components
           IF (IAND(tver,vir).NE.0) THEN
              g_vxx = g_vxx + xad*dedxia + xbd*dedxib + xcd*dedxic
              g_vxy = g_vxy + yad*dedxia + ybd*dedxib + ycd*dedxic
              g_vxz = g_vxz + zad*dedxia + zbd*dedxib + zcd*dedxic
              g_vyy = g_vyy + yad*dedyia + ybd*dedyib + ycd*dedyic
              g_vyz = g_vyz + zad*dedyia + zbd*dedyib + zcd*dedyic
              g_vzz = g_vzz + zad*dedzia + zbd*dedzib + zcd*dedzic
           END IF
        end if
        end subroutine
