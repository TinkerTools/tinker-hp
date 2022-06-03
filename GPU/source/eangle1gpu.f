c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangle1  --  angle bend energy and derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangle1" calculates the angle bending potential energy and
c     the first derivatives with respect to Cartesian coordinates;
c     projected in-plane angles at trigonal centers, special linear
c     or Fourier angle bending terms are optionally used
c
c
c
c    Comment
c    Zhi Wang, July 1, 2019
c    
c    The original Tinker implementation has following expressions
c    xip = xib + xt * delta
c    xap = xia - xip
c    
c    And they were reorganized to
c    xap = xia - xib - xt * delta
c    for higher accuracy in the single precision mode.
c    
c    Consider an example where
c    xia = 33.553368, xib = 34.768604
c    xt = 0.33142909, delta = 0.0044494048,
c    the later expression gives a better numerical result.
c
#include "tinker_precision.h"
      module eangle1gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine eangle1gpu
      use angle
      use angpot
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eangle1gpu_inl
      use group
      use inform,only: deb_Path
      use math
      use potent,only:use_amd_wat1
      use usage
      use virial
      use timestat,only:timer_enter,timer_exit,timer_eangle
      use tinheader
      use mamd
      implicit none
      integer i,ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      integer iangle
      real(t_p) e,ideal,force
      real(t_p) fold,factor,dot
      real(t_p) cosine,sine
      real(t_p) angle1
      real(t_p) dt,dt2,dt3,dt4
      real(t_p) deddt,term
      real(t_p) terma,termc
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) xp,yp,zp,rp
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xcd,ycd,zcd
      real(t_p) xip,yip,zip
      real(t_p) xap,yap,zap
      real(t_p) xcp,ycp,zcp
      real(t_p) rab2,rcb2
      real(t_p) rap2,rcp2
      real(t_p) xt,yt,zt
      real(t_p) rt2,ptrt2
      real(t_p) xm,ym,zm,rm
      real(t_p) delta,delta2
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(t_p) dedxip,dedyip,dedzip
      real(t_p) dpdxia,dpdyia,dpdzia
      real(t_p) dpdxic,dpdyic,dpdzic
      real(r_p) dedi(3),dpdi(3)
      logical proceed

      if(deb_Path) write(*,*) 'eangle1gpu'
      call timer_enter( timer_eangle )
c
c     calculate the bond angle bending energy term
c
!     Force reduction directive because implicit reduction is not
!     detected by PGI-20 compiler

!$acc parallel loop private(dedi) async
#ifdef USE_NVSHMEM_CUDA
!$acc&     present(x,y,z,loc,use,angleglob,anat,
#else
!$acc&     present(x,y,z,loc,use,iang,angleglob,anat,
#endif
!$acc&     ak,afld,angtypI,dea,vir,ea,eW1aMD)
!$acc&     present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&  reduction(+:ea,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do iangle = 1, nangleloc
         i  = angleglob(iangle)
#ifdef USE_NVSHMEM_CUDA
         ipe   =     (i-1)/nangle_pe
         ind   = mod((i-1),nangle_pe) +1
         ia    = d_iang(ipe)%pel(1,ind)
         ib    = d_iang(ipe)%pel(2,ind)
         ic    = d_iang(ipe)%pel(3,ind)
         id    = d_iang(ipe)%pel(4,ind)
#else
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
#endif
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         ideal = anat(i)
         force = ak(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (angtypI(i) .eq. ANG_IN_PLANE) then
            idloc = loc(id)
            if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                               use(ic) .or. use(id))
         else
            if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
         end if
c
c     get the coordinates of the atoms in the angle
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
c
c     compute the bond angle bending energy and gradient
c
            if (angtypI(i) .ne. ANG_IN_PLANE) then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               if (use_polymer) then
                  call image_inl (xab,yab,zab)
                  call image_inl (xcb,ycb,zcb)
               end if
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
                  angle1 = radian * acos(cosine)
c
c     get the energy and master chain rule term for derivatives
c
                  if (angtypI(i) .eq. ANG_HARMONIC) then
                     dt = angle1 - ideal
                     dt2 = dt * dt
                     dt3 = dt2 * dt
                     dt4 = dt2 * dt2
                     e = angunit * force * dt2
     &                   * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                     deddt = angunit * force * dt * radian
     &                * (2.0_ti_p + 3.0_ti_p*cang*dt + 4.0_ti_p*qang*dt2
     &                          + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
                  else if (angtypI(i) .eq. ANG_LINEAR) then
                     factor = 2.0_ti_p * angunit * radian**2
                     sine = sqrt(1.0_ti_p-cosine*cosine)
                     e = factor * force * (1.0_ti_p+cosine)
                     deddt = -factor * force * sine
                  else if (angtypI(i) .eq. ANG_FOURIER) then
                     fold = afld(i)
                     factor = 2.0_ti_p * angunit * (radian/fold)**2
                     cosine = cos((fold*angle1-ideal)/radian)
                     sine = sin((fold*angle1-ideal)/radian)
                     e = factor * force * (1.0_ti_p+cosine)
                     deddt = -factor * force * fold * sine
                  end if
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
c     increment the total bond angle energy and derivatives
c
                  ea = ea + e

                  dedi(1) = dedxib
                  dedi(2) = dedyib
                  dedi(3) = dedzib
!$acc atomic update 
                  dea(1,ibloc) = dea(1,ibloc) + dedi(1)
!$acc atomic update 
                  dea(2,ibloc) = dea(2,ibloc) + dedi(2)
!$acc atomic update 
                  dea(3,ibloc) = dea(3,ibloc) + dedi(3)
c
                  dedi(1) = dedxia
                  dedi(2) = dedyia
                  dedi(3) = dedzia
!$acc atomic update 
                  dea(1,ialoc) = dea(1,ialoc) + dedi(1)
!$acc atomic update 
                  dea(2,ialoc) = dea(2,ialoc) + dedi(2)
!$acc atomic update 
                  dea(3,ialoc) = dea(3,ialoc) + dedi(3)
c
                  dedi(1) = dedxic
                  dedi(2) = dedyic
                  dedi(3) = dedzic
!$acc atomic update 
                  dea(1,icloc) = dea(1,icloc) + dedi(1)
!$acc atomic update 
                  dea(2,icloc) = dea(2,icloc) + dedi(2)
!$acc atomic update 
                  dea(3,icloc) = dea(3,icloc) + dedi(3)
c
c     increment the internal virial tensor components
c
                  g_vxx = g_vxx + xab*dedxia + xcb*dedxic
                  g_vxy = g_vxy + yab*dedxia + ycb*dedxic
                  g_vxz = g_vxz + zab*dedxia + zcb*dedxic
                  g_vyy = g_vyy + yab*dedyia + ycb*dedyic
                  g_vyz = g_vyz + zab*dedyia + zcb*dedyic
                  g_vzz = g_vzz + zab*dedzia + zcb*dedzic
               end if
c
c     compute the projected in-plane angle energy and gradient
c
            else
               xid = x(id)
               yid = y(id)
               zid = z(id)
               xad = xia - xid
               yad = yia - yid
               zad = zia - zid
               xbd = xib - xid
               ybd = yib - yid
               zbd = zib - zid
               xcd = xic - xid
               ycd = yic - yid
               zcd = zic - zid
               if (use_polymer) then
                  call image_inl (xad,yad,zad)
                  call image_inl (xbd,ybd,zbd)
                  call image_inl (xcd,ycd,zcd)
               end if
               xt = yad*zcd - zad*ycd
               yt = zad*xcd - xad*zcd
               zt = xad*ycd - yad*xcd
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
               if (use_polymer) then
                  call image_inl (xap,yap,zap)
                  call image_inl (xcp,ycp,zcp)
               end if
               rap2 = xap*xap + yap*yap + zap*zap
               rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
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
                  dt = angle1 - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                  deddt = angunit * force * dt * radian
     &                * (2.0_ti_p + 3.0_ti_p*cang*dt + 4.0_ti_p*qang*dt2
     &                        + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
c
c     chain rule terms for first derivative components
c
                  terma = -deddt / (rap2*rm)
                  termc = deddt / (rcp2*rm)
                  dedxia = terma * (yap*zm-zap*ym)
                  dedyia = terma * (zap*xm-xap*zm)
                  dedzia = terma * (xap*ym-yap*xm)
                  dedxic = termc * (ycp*zm-zcp*ym)
                  dedyic = termc * (zcp*xm-xcp*zm)
                  dedzic = termc * (xcp*ym-ycp*xm)
                  dedxip = -dedxia - dedxic
                  dedyip = -dedyia - dedyic
                  dedzip = -dedzia - dedzic
c
c     chain rule components for the projection of the central atom
c
                  delta2 = 2.0_ti_p * delta
                  ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
                  term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
                  dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
                  term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
                  dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
                  term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
                  dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
                  term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
                  dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
                  term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
                  dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
                  term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
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
c
c     increment the total bond angle energy and derivatives
c
                  ea = ea + e

                  dedi(1) = dedxib
                  dedi(2) = dedyib
                  dedi(3) = dedzib
!$acc atomic update
                  dea(1,ibloc) = dea(1,ibloc) + dedi(1)
!$acc atomic update
                  dea(2,ibloc) = dea(2,ibloc) + dedi(2)
!$acc atomic update
                  dea(3,ibloc) = dea(3,ibloc) + dedi(3)
c
                  dedi(1) = dedxia
                  dedi(2) = dedyia
                  dedi(3) = dedzia
!$acc atomic update
                  dea(1,ialoc) = dea(1,ialoc) + dedi(1)
!$acc atomic update
                  dea(2,ialoc) = dea(2,ialoc) + dedi(2)
!$acc atomic update
                  dea(3,ialoc) = dea(3,ialoc) + dedi(3)
c
                  dedi(1) = dedxic
                  dedi(2) = dedyic
                  dedi(3) = dedzic
!$acc atomic update
                  dea(1,icloc) = dea(1,icloc) + dedi(1)
!$acc atomic update
                  dea(2,icloc) = dea(2,icloc) + dedi(2)
!$acc atomic update
                  dea(3,icloc) = dea(3,icloc) + dedi(3)
c
                  dedi(1) = dedxid
                  dedi(2) = dedyid
                  dedi(3) = dedzid
!$acc atomic update
                  dea(1,idloc) = dea(1,idloc) + dedi(1)
!$acc atomic update
                  dea(2,idloc) = dea(2,idloc) + dedi(2)
!$acc atomic update
                  dea(3,idloc) = dea(3,idloc) + dedi(3)
c
c     aMD storage if waters are considered
c
                  if (use_amd_wat1) then  ! Check amd
                  if (type(ia) == aMDwattype(1) .or. type(ib)
     $            == aMDwattype(1) .or. type(ic) == aMDwattype(1))
     $            then   !Check type
                     eW1aMD = eW1aMD + e
                     dedi(1) = dedxia
                     dedi(2) = dedyia
                     dedi(3) = dedzia
!$acc atomic
                     deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedi(1)
!$acc atomic
                     deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedi(2)
!$acc atomic
                     deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedi(3)
                     dedi(1) = dedxib
                     dedi(2) = dedyib
                     dedi(3) = dedzib
!$acc atomic
                     deW1aMD(1,ibloc) = deW1aMD(1,ibloc) + dedi(1)
!$acc atomic
                     deW1aMD(2,ibloc) = deW1aMD(2,ibloc) + dedi(2)
!$acc atomic
                     deW1aMD(3,ibloc) = deW1aMD(3,ibloc) + dedi(3)
                     dedi(1) = dedxic
                     dedi(2) = dedyic
                     dedi(3) = dedzic
!$acc atomic
                     deW1aMD(1,icloc) = deW1aMD(1,icloc) + dedi(1)
!$acc atomic
                     deW1aMD(2,icloc) = deW1aMD(2,icloc) + dedi(2)
!$acc atomic
                     deW1aMD(3,icloc) = deW1aMD(3,icloc) + dedi(3)
                  end if  !Check type
                  end if  !Check amd
c
c     increment the internal virial tensor components
c
                g_vxx = g_vxx + xad*dedxia + xbd*dedxib + xcd*dedxic
                g_vxy = g_vxy + yad*dedxia + ybd*dedxib + ycd*dedxic
                g_vxz = g_vxz + zad*dedxia + zbd*dedxib + zcd*dedxic
                g_vyy = g_vyy + yad*dedyia + ybd*dedyib + ycd*dedyic
                g_vyz = g_vyz + zad*dedyia + zbd*dedyib + zcd*dedyic
                g_vzz = g_vzz + zad*dedzia + zbd*dedzib + zcd*dedzic
               end if
            end if
         end if
      end do
      call timer_exit( timer_eangle )
      end
