c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine egeom1  --  restraint energy & derivatives  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "egeom1" calculates the energy and first derivatives
c     with respect to Cartesian coordinates due to restraints
c     on positions, distances, angles and torsions as well as
c     Gaussian basin and spherical droplet restraints
c
c
#include "tinker_precision.h"
      module egeom1gpu_inl
        contains
#include "image.f.inc"
#include "midpointimage.f.inc"
      end module

      subroutine egeom1gpu
      use atmlst
      use atmtyp
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use egeom1gpu_inl
      use group
      use inter
      use kgeoms
      use math
      use molcul
      use tinheader ,only:ti_p,re_p
      use usage
      use USampling ,only:USdt=>timestep,US_save,Rd_save
     &              ,cpt_wh,US_enable,USwrite,step,step_save
      use virial
      implicit none
      integer i,j,k,sl
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer iglob,kglob
      integer inpfix,indfix,inafix,intfix,ingfix,inchir,kloc
      logical docompute,save_US
      real(t_p) xk,yk,zk
      real(t_p) e,xr,yr,zr
      real(t_p) de,dt,dt2,deddt
      real(t_p) r,r2,r6,r12
      real(r_p) dedx,dedy,dedz
      real(t_p) angle,target
      real(t_p) dot,force
      real(t_p) cosine,sine
      real(t_p) terma,termc
      real(t_p) rab2,rcb2
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xab,yab,zab
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      real(t_p) xca,yca,zca
      real(t_p) xdb,ydb,zdb
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xcd,ycd,zcd
      real(t_p) xp,yp,zp,rp
      real(t_p) xt,yt,zt
      real(t_p) xu,yu,zu
      real(t_p) xtu,ytu,ztu
      real(t_p) rt2,ru2,rtru
      real(t_p) rcb,dedphi
      real(t_p) dedxt,dedyt,dedzt
      real(t_p) dedxu,dedyu,dedzu
      real(r_p) dedxia,dedyia,dedzia
      real(r_p) dedxib,dedyib,dedzib
      real(r_p) dedxic,dedyic,dedzic
      real(r_p) dedxid,dedyid,dedzid
      real(t_p) df1,df2
      real(t_p) af1,af2
      real(t_p) tf1,tf2,t1,t2
      real(t_p) gf1,gf2
      real(t_p) weigh,ratio
      real(t_p) weigha,weighb
      real(t_p) xcm,ycm,zcm
      real(t_p) cf1,cf2,vol
      real(t_p) c1,c2,c3
      real(t_p) xi,yi,zi,ri
      real(t_p) a,b,buffer,term
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed,intermol
c
c
c     zero out the restraint energy term and first derivatives
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'egeom1gpu'
      eg = 0.0_ti_p
      save_US = .false.
      if (US_enable) save_US=(US_enable.and.mod(step,USwrite).eq.0
     &                                 .and.step.ne.step_save)

c
!$acc data present(mass,xpfix,ypfix,zpfix,kpfix,    
!$acc&  use,npfixglob,ndfixglob,nafixglob,ntfixglob,ngfixglob,
!$acc&  nchirglob,igrp,kgrp,loc,x,y,z)


!$acc parallel vector_length(32) default(present) async
!$acc&     present(eg,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,einter)
!$acc& reduction(+:eg,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

!$acc loop gang vector
      do inpfix = 1, npfixloc
         i = npfixglob(inpfix)
         ia = ipfix(i)
         ialoc = loc(ia)
         proceed = .true.
         if (proceed)  proceed = (use(ia))
         if (proceed) then
            xr = 0.0_ti_p
            yr = 0.0_ti_p
            zr = 0.0_ti_p
            if (kpfix(1,i) .ne. 0)  xr = x(ia) - xpfix(i)
            if (kpfix(2,i) .ne. 0)  yr = y(ia) - ypfix(i)
            if (kpfix(3,i) .ne. 0)  zr = z(ia) - zpfix(i)
            if (use_bounds)  call image_inl (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            force = pfix(1,i)
            dt = max(0.0_ti_p,r-pfix(2,i))
            dt2 = dt * dt
            e = force * dt2
            if (r .eq. 0.0_ti_p)  r = 1.0_ti_p
            de = 2.0_ti_p * force * dt / r
c
c     compute chain rule terms needed for derivatives
c
            dedx = de * xr
            dedy = de * yr
            dedz = de * zr
c
c     increment the total energy and first derivatives
c
            eg = eg + e
!$acc atomic update
            deg(1,ialoc) = deg(1,ialoc) + dedx
!$acc atomic update
            deg(2,ialoc) = deg(2,ialoc) + dedy
!$acc atomic update
            deg(3,ialoc) = deg(3,ialoc) + dedz
c
c     increment the internal virial tensor components
c
            g_vxx = g_vxx + xr * dedx
            g_vxy = g_vxy + yr * dedx
            g_vxz = g_vxz + zr * dedx
            g_vyy = g_vyy + yr * dedy
            g_vyz = g_vyz + zr * dedy
            g_vzz = g_vzz + zr * dedz
         end if
      end do


c
c     get energy and derivatives for distance restraint terms
c
!$acc loop gang vector
      do indfix = 1, ndfixloc
         i = ndfixglob(indfix)
         ia = idfix(1,i)
         ialoc = loc(ia)
         ib = idfix(2,i)
         ibloc = loc(ib)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib))
         if (proceed) then
            xr = x(ia) - x(ib)
            yr = y(ia) - y(ib)
            zr = z(ia) - z(ib)
            intermol = (molcule(ia) .ne. molcule(ib))
            if (use_bounds .and. intermol) call image_inl(xr,yr,zr)
c            if (use_bounds)  call image (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            force = dfix(1,i)
            df1 = dfix(2,i)
            df2 = dfix(3,i)
            target = r
            if (r .lt. df1)  target = df1
            if (r .gt. df2)  target = df2
            dt = r - target
            dt2 = dt * dt
            e = force * dt2
            if (r .eq. 0.0_ti_p)  r = 1.0_ti_p
            de = 2.0_ti_p * force * dt / r
            if (save_US) then
               !print*,"save us in array"
               sl = (cpt_wh-1)*ndfix*4 + (i-1)*4
               Rd_save(sl+1) = i
               Rd_save(sl+2) = step*USdt
               Rd_save(sl+3) = r
               Rd_save(sl+4) = e
            end if
c
c     compute chain rule terms needed for derivatives
c
            dedx = de * xr
            dedy = de * yr
            dedz = de * zr
c
c     increment the total energy and first derivatives
c
            eg = eg + e
!$acc atomic update
            deg(1,ialoc) = deg(1,ialoc) + dedx
!$acc atomic update
            deg(2,ialoc) = deg(2,ialoc) + dedy
!$acc atomic update
            deg(3,ialoc) = deg(3,ialoc) + dedz
c
!$acc atomic update
            deg(1,ibloc) = deg(1,ibloc) - dedx
!$acc atomic update
            deg(2,ibloc) = deg(2,ibloc) - dedy
!$acc atomic update
            deg(3,ibloc) = deg(3,ibloc) - dedz
c
c     increment the internal virial tensor components
c
            g_vxx = g_vxx + xr * dedx
            g_vxy = g_vxy + yr * dedx
            g_vxz = g_vxz + zr * dedx
            g_vyy = g_vyy + yr * dedy
            g_vyz = g_vyz + zr * dedy
            g_vzz = g_vzz + zr * dedz
c
c     increment the total intermolecular energy
c
c           if (intermol) then
c              einter = einter + e
c           end if
         end if
      end do


c
c     get energy and derivatives for angle restraint terms
c

!$acc loop gang vector
      do inafix = 1, nafixloc
         i = nafixglob(inafix)
         ia = iafix(1,i)
         ialoc = loc(ia)
         ib = iafix(2,i)
         ibloc = loc(ib)
         ic = iafix(3,i)
         icloc = loc(ic)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               rp = max(rp,0.000001_ti_p)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle = radian * acos(cosine)
               force = afix(1,i)
               af1 = afix(2,i)
               af2 = afix(3,i)
               target = angle
               if (angle .lt. af1)  target = af1
               if (angle .gt. af2)  target = af2
               dt = angle - target
               dt2 = dt * dt
               e = force * dt2
               deddt = 2.0_ti_p * force * dt * radian
c
c     compute derivative components for this interaction
c
               terma = -deddt / (rab2*rp)
               termc = deddt / (rcb2*rp)
               dedxia = terma * (yab*zp-zab*yp)
               dedyia = terma * (zab*xp-xab*zp)
               dedzia = terma * (xab*yp-yab*xp)
               dedxic = termc * (ycb*zp-zcb*yp)
               dedyic = termc * (zcb*xp-xcb*zp)
               dedzic = termc * (xcb*yp-ycb*xp)
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the overall energy term and derivatives
c
               eg = eg + e
!$acc atomic update  
               deg(1,ialoc) = deg(1,ialoc) + dedxia
!$acc atomic update  
               deg(2,ialoc) = deg(2,ialoc) + dedyia
!$acc atomic update  
               deg(3,ialoc) = deg(3,ialoc) + dedzia
c
!$acc atomic update  
               deg(1,ibloc) = deg(1,ibloc) + dedxib
!$acc atomic update  
               deg(2,ibloc) = deg(2,ibloc) + dedyib
!$acc atomic update  
               deg(3,ibloc) = deg(3,ibloc) + dedzib
c
!$acc atomic update  
               deg(1,icloc) = deg(1,icloc) + dedxic
!$acc atomic update  
               deg(2,icloc) = deg(2,icloc) + dedyic
!$acc atomic update  
               deg(3,icloc) = deg(3,icloc) + dedzic
c
c     increment the internal virial tensor components
c
               g_vxx = g_vxx + xab*dedxia + xcb*dedxic
               g_vxy = g_vxy + yab*dedxia + ycb*dedxic
               g_vxz = g_vxz + zab*dedxia + zcb*dedxic
               g_vyy = g_vyy + yab*dedyia + ycb*dedyic
               g_vyz = g_vyz + zab*dedyia + zcb*dedyic
               g_vzz = g_vzz + zab*dedzia + zcb*dedzic
c
c     increment the total intermolecular energy
c
c              if (intermol) then
c                 einter = einter + e
c              end if
            end if
         end if
      end do
c
c     get energy and derivatives for torsion restraint terms
c
!$acc loop gang vector
      do intfix = 1, ntfixloc
         i = ntfixglob(intfix)
         ia = itfix(1,i)
         ialoc = loc(ia)
         ib = itfix(2,i)
         ibloc = loc(ib)
         ic = itfix(3,i)
         icloc = loc(ic)
         id = itfix(4,i)
         idloc = loc(id)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
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
               force = tfix(1,i)
               tf1 = tfix(2,i)
               tf2 = tfix(3,i)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0_ti_p) then
                     t1 = t1 - 360.0_ti_p
                  else if (t1 .lt. -180.0_ti_p) then
                     t1 = t1 + 360.0_ti_p
                  end if
                  if (t2 .gt. 180.0_ti_p) then
                     t2 = t2 - 360.0_ti_p
                  else if (t2 .lt. -180.0_ti_p) then
                     t2 = t2 + 360.0_ti_p
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               dt = angle - target
               if (dt .gt. 180.0_ti_p) then
                  dt = dt - 360.0_ti_p
               else if (dt .lt. -180.0_ti_p) then
                  dt = dt + 360.0_ti_p
               end if
               dt2 = dt * dt
               e = force * dt2
               dedphi = 2.0_ti_p * radian * force * dt
c
c     chain rule terms for first derivative components
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute derivative components for this interaction
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
c     increment the overall energy term and derivatives
c
               eg = eg + e
!$acc atomic update
               deg(1,ialoc) = deg(1,ialoc) + dedxia
!$acc atomic update
               deg(2,ialoc) = deg(2,ialoc) + dedyia
!$acc atomic update
               deg(3,ialoc) = deg(3,ialoc) + dedzia
c
!$acc atomic update
               deg(1,ibloc) = deg(1,ibloc) + dedxib
!$acc atomic update
               deg(2,ibloc) = deg(2,ibloc) + dedyib
!$acc atomic update
               deg(3,ibloc) = deg(3,ibloc) + dedzib
c
!$acc atomic update
               deg(1,icloc) = deg(1,icloc) + dedxic
!$acc atomic update
               deg(2,icloc) = deg(2,icloc) + dedyic
!$acc atomic update
               deg(3,icloc) = deg(3,icloc) + dedzic
c
!$acc atomic update
               deg(1,idloc) = deg(1,idloc) + dedxid
!$acc atomic update
               deg(2,idloc) = deg(2,idloc) + dedyid
!$acc atomic update
               deg(3,idloc) = deg(3,idloc) + dedzid
c
c     increment the internal virial tensor components
c
           g_vxx = g_vxx + xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
           g_vxy = g_vxy + ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
           g_vxz = g_vxz + zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
           g_vyy = g_vyy + ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
           g_vyz = g_vyz + zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
           g_vzz = g_vzz + zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
c
c     increment the total intermolecular energy
c
c              if (molcule(ia).ne.molcule(ib) .or.
c    &             molcule(ia).ne.molcule(ic) .or.
c    &             molcule(ia).ne.molcule(id)) then
c                 einter = einter + e
c              end if
            end if
         end if
      end do
c
c     get energy and derivatives for group distance restraint terms
c
!$acc loop gang private(r,e)
      do ingfix = 1, ngfixloc
         i = ngfixglob(ingfix)
         ia = igfix(1,i)
         ialoc = loc(ia)
         ib = igfix(2,i)
         ibloc = loc(ib)
         xcm = 0.0_ti_p
         ycm = 0.0_ti_p
         zcm = 0.0_ti_p
!$acc loop vector
         do j = igrp(1,ia), igrp(2,ia)
           k = kgrp(j)
           weigh = mass(k)
           xcm = xcm + x(k)*weigh
           ycm = ycm + y(k)*weigh
           zcm = zcm + z(k)*weigh
         end do
         weigha = max(1.0_ti_p,grpmass(ia))
         xr = xcm / weigha
         yr = ycm / weigha
         zr = zcm / weigha
         xcm = 0.0_ti_p
         ycm = 0.0_ti_p
         zcm = 0.0_ti_p
!$acc loop vector
         do j = igrp(1,ib), igrp(2,ib)
           k = kgrp(j)
           weigh = mass(k)
           xcm = xcm + x(k)*weigh
           ycm = ycm + y(k)*weigh
           zcm = zcm + z(k)*weigh
         end do
         weighb = max(1.0_ti_p,grpmass(ib))
         xr = xr - xcm/weighb
         yr = yr - ycm/weighb
         zr = zr - zcm/weighb
         intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &               molcule(kgrp(igrp(1,ib))))
         if (use_bounds .and. intermol) call image_inl(xr,yr,zr)
c         if (use_bounds)  call image (xr,yr,zr)
         r = sqrt(xr*xr + yr*yr + zr*zr)

         force = gfix(1,i)
         gf1 = gfix(2,i)
         gf2 = gfix(3,i)
         target = r
         if (r .lt. gf1)  target = gf1
         if (r .gt. gf2)  target = gf2
         dt = r - target
         dt2 = dt * dt
         e = force * dt2
         if (r .eq. 0.0_ti_p)  r = 1.0_ti_p
         de = 2.0_ti_p * force * dt / r
         ! Save data for Fred

!$acc loop vector
         do j = 1,1
            if (save_US) then
               !print*,"save us in array"
               sl = (cpt_wh-1)*ngfix*4 + (i-1)*4
               US_save(sl+1) = i
               US_save(sl+2) = step*USdt
               US_save(sl+3) = r
               US_save(sl+4) = e
            end if
         end do
c
c     compute chain rule terms needed for derivatives
c
         dedx = de * xr
         dedy = de * yr
         dedz = de * zr
!$acc loop vector
         do j = igrp(1,ia), igrp(2,ia)
            k = kgrp(j)
            kloc = loc(k)
            ratio = mass(k) / weigha
!$acc atomic update  
            deg(1,kloc) = deg(1,kloc) + dedx*ratio
!$acc atomic update  
            deg(2,kloc) = deg(2,kloc) + dedy*ratio
!$acc atomic update  
            deg(3,kloc) = deg(3,kloc) + dedz*ratio
         end do
!$acc loop vector
         do j = igrp(1,ib), igrp(2,ib)
            k = kgrp(j)
            kloc = loc(k)
            ratio = mass(k) / weighb
!$acc atomic update 
            deg(1,kloc) = deg(1,kloc) - dedx*ratio
!$acc atomic update 
            deg(2,kloc) = deg(2,kloc) - dedy*ratio
!$acc atomic update 
            deg(3,kloc) = deg(3,kloc) - dedz*ratio
         end do

!$acc loop vector
         do j = 1,1
c
c     increment the total energy and first derivatives
c
         eg = eg + e
c
c     increment the internal virial tensor components
c
         g_vxx = g_vxx + xr * dedx
         g_vxy = g_vxy + yr * dedx
         g_vxz = g_vxz + zr * dedx
         g_vyy = g_vyy + yr * dedy
         g_vyz = g_vyz + zr * dedy
         g_vzz = g_vzz + zr * dedz
         end do
c
c     increment the total intermolecular energy
c
c        if (intermol) then
c           einter = einter + e
c        end if
      end do
c
c     get energy and derivatives for chirality restraint terms
c
!$acc loop gang vector
      do inchir = 1, nchirloc
         i = nchirglob(inchir)
         ia = ichir(1,i)
         ialoc = loc(ia)
         ib = ichir(2,i)
         ibloc = loc(ib)
         ic = ichir(3,i)
         icloc = loc(ic)
         id = ichir(4,i)
         idloc = loc(id)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
         if (proceed) then
            xad = x(ia) - x(id)
            yad = y(ia) - y(id)
            zad = z(ia) - z(id)
            xbd = x(ib) - x(id)
            ybd = y(ib) - y(id)
            zbd = z(ib) - z(id)
            xcd = x(ic) - x(id)
            ycd = y(ic) - y(id)
            zcd = z(ic) - z(id)
            c1 = ybd*zcd - zbd*ycd
            c2 = ycd*zad - zcd*yad
            c3 = yad*zbd - zad*ybd
            vol = xad*c1 + xbd*c2 + xcd*c3
            force = chir(1,i)
            cf1 = chir(2,i)
            cf2 = chir(3,i)
            target = vol
            if (vol .lt. min(cf1,cf2))  target = min(cf1,cf2)
            if (vol .gt. max(cf1,cf2))  target = max(cf1,cf2)
            dt = vol - target
            dt2 = dt * dt
            e = force * dt2
            deddt = 2.0_ti_p * force * dt
c
c     compute derivative components for this interaction
c
            dedxia = deddt * (ybd*zcd - zbd*ycd)
            dedyia = deddt * (zbd*xcd - xbd*zcd)
            dedzia = deddt * (xbd*ycd - ybd*xcd)
            dedxib = deddt * (zad*ycd - yad*zcd)
            dedyib = deddt * (xad*zcd - zad*xcd)
            dedzib = deddt * (yad*xcd - xad*ycd)
            dedxic = deddt * (yad*zbd - zad*ybd)
            dedyic = deddt * (zad*xbd - xad*zbd)
            dedzic = deddt * (xad*ybd - yad*xbd)
            dedxid = -dedxia - dedxib - dedxic
            dedyid = -dedyia - dedyib - dedyic
            dedzid = -dedzia - dedzib - dedzic
c
c     increment the overall energy term and derivatives
c
            eg = eg + e
!$acc atomic update  
            deg(1,ialoc) = deg(1,ialoc) + dedxia
!$acc atomic update  
            deg(2,ialoc) = deg(2,ialoc) + dedyia
!$acc atomic update  
            deg(3,ialoc) = deg(3,ialoc) + dedzia
c
!$acc atomic update  
            deg(1,ibloc) = deg(1,ibloc) + dedxib
!$acc atomic update  
            deg(2,ibloc) = deg(2,ibloc) + dedyib
!$acc atomic update  
            deg(3,ibloc) = deg(3,ibloc) + dedzib
c
!$acc atomic update  
            deg(1,icloc) = deg(1,icloc) + dedxic
!$acc atomic update  
            deg(2,icloc) = deg(2,icloc) + dedyic
!$acc atomic update  
            deg(3,icloc) = deg(3,icloc) + dedzic
c
!$acc atomic update  
            deg(1,idloc) = deg(1,idloc) + dedxid
!$acc atomic update  
            deg(2,idloc) = deg(2,idloc) + dedyid
!$acc atomic update  
            deg(3,idloc) = deg(3,idloc) + dedzid
c
c     increment the internal virial tensor components
c
            g_vxx = g_vxx + xad*dedxia + xbd*dedxib + xcd*dedxic
            g_vxy = g_vxy + yad*dedxia + ybd*dedxib + ycd*dedxic
            g_vxz = g_vxz + zad*dedxia + zbd*dedxib + zcd*dedxic
            g_vyy = g_vyy + yad*dedyia + ybd*dedyib + ycd*dedyic
            g_vyz = g_vyz + zad*dedyia + zbd*dedyib + zcd*dedyic
            g_vzz = g_vzz + zad*dedzia + zbd*dedzib + zcd*dedzic
c
c     increment the total intermolecular energy
c
c           if (molcule(ia).ne.molcule(ib) .or.
c    &          molcule(ia).ne.molcule(ic) .or.
c    &          molcule(ia).ne.molcule(id)) then
c              einter = einter + e
c           end if
         end if
      end do
!$acc end parallel

      if(US_enable) then
         step_save=step
         if (save_US) cpt_wh = cpt_wh + 1
      end if

c
c     get energy and derivatives for a Gaussian basin restraint
c
c      print *,"use_basin",use_basin
      if (use_basin) then
!$acc parallel loop default(present) async
!$acc&         present(eg,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     reduction(+:eg,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         do i = 1, nloc
            iglob = glob(i)
            xi = x(iglob)
            yi = y(iglob)
            zi = z(iglob)
!$acc loop vector
            do k = 1, nbloc
               kglob = glob(k)
               xk = x(kglob)
               yk = y(kglob)
               zk = z(kglob)
               proceed = .true.
               if (proceed)  proceed = (use(iglob) .or. use(kglob))
               if (proceed) then
                  if (kglob.le.iglob) cycle

                  xr = xi - xk
                  yr = yi - yk
                  zr = zi - zk
                  call midpoint_inl(xk,yk,zk,xr,yr,zr,docompute)
                  if (.not.(docompute)) cycle
                  r2 = xr*xr + yr*yr + zr*zr
                  term = -width * r2
                  e = 0.0_ti_p
                  if (term .gt. -50.0_ti_p)  e = depth * exp(term)
                  de = -2.0_ti_p * width * e
                  e = e - depth
c
c     compute chain rule terms needed for derivatives
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy term and derivatives
c
                  eg = eg + e
!$acc atomic update  
                  deg(1,i) = deg(1,i) + dedx
!$acc atomic update  
                  deg(2,i) = deg(2,i) + dedy
!$acc atomic update  
                  deg(3,i) = deg(3,i) + dedz
c
!$acc atomic update  
                  deg(1,k) = deg(1,k) - dedx
!$acc atomic update  
                  deg(2,k) = deg(2,k) - dedy
!$acc atomic update  
                  deg(3,k) = deg(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  g_vxx = g_vxx + xr * dedx
                  g_vxy = g_vxy + yr * dedx
                  g_vxz = g_vxz + zr * dedx
                  g_vyy = g_vyy + yr * dedy
                  g_vyz = g_vyz + zr * dedy
                  g_vzz = g_vzz + zr * dedz
               end if
            end do
         end do
      end if

c
c     get energy and derivatives for a spherical droplet restraint
c
c      print *,"use-wall",use_wall
      if (use_wall) then
         buffer = 2.5_ti_p
         a = 2048.0_ti_p
         b = 64.0_ti_p
!$acc parallel loop default(present) async
!$acc&         present(eg,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     reduction(+:eg,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         do i = 1, nloc
            iglob = glob(i)
            proceed = .true.
            if (proceed)  proceed = (use(iglob))
            if (proceed) then
               xi = x(iglob)
               yi = y(iglob)
               zi = z(iglob)
               ri = sqrt(xi**2 + yi**2 + zi**2)
               r = rwall + buffer - ri
               r2 = r * r
               r6 = r2 * r2 * r2
               r12 = r6 * r6
               e = a/r12 - b/r6
               if (ri .eq. 0.0_ti_p)  ri = 1.0_ti_p
               de = (12.0_ti_p*a/r12 - 6.0_ti_p*b/r6) / (r*ri)
c
c     compute chain rule terms needed for derivatives
c
               dedx = de * xi
               dedy = de * yi
               dedz = de * zi
c
c     increment the overall energy term and derivatives
c
               eg = eg + e
!$acc atomic update
               deg(1,i) = deg(1,i) + dedx
!$acc atomic update
               deg(2,i) = deg(2,i) + dedy
!$acc atomic update
               deg(3,i) = deg(3,i) + dedz
c
c     increment the internal virial tensor components
c
               xr = r * xi/ri
               yr = r * yi/ri
               zr = r * zi/ri
               g_vxx = g_vxx + xr * dedx
               g_vxy = g_vxy + yr * dedx
               g_vxz = g_vxz + zr * dedx
               g_vyy = g_vyy + yr * dedy
               g_vyz = g_vyz + zr * dedy
               g_vzz = g_vzz + zr * dedz
            end if
         end do
      end if

!$acc end data
c
      return
      end
