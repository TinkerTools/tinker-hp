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
#include "tinker_precision.h"
      subroutine eangle1
      use angle
      use angpot
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use group
      use math
      use mamd
      use potent,only:use_amd_wat1
      use usage
      use tinheader
      use virial
      implicit none
      integer i,ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
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
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed
c
!$acc update host(dea,vir)
c
c     zero out energy and first derivative components
c
      ea = 0.0_ti_p
c
c     calculate the bond angle bending energy term
c
      do iangle = 1, nangleloc
         i = angleglob(iangle)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
         ideal = anat(i)
         force = ak(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (angtyp(i) .eq. 'IN-PLANE') then
            if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                                 use(ic) .or. use(id))
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
            if (angtyp(i) .ne. 'IN-PLANE') then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               if (use_polymer) then
                  call image (xab,yab,zab)
                  call image (xcb,ycb,zcb)
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
                  if (angtyp(i) .eq. 'HARMONIC') then
                     dt = angle1 - ideal
                     dt2 = dt * dt
                     dt3 = dt2 * dt
                     dt4 = dt2 * dt2
                     e = angunit * force * dt2
     &                   * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                     deddt = angunit * force * dt * radian
     &                     * (2.0_ti_p + 3.0_ti_p*cang*dt + 
     &                                                4.0_ti_p*qang*dt2
     &                     + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
                  else if (angtyp(i) .eq. 'LINEAR') then
                     factor = 2.0_ti_p * angunit * radian**2
                     sine = sqrt(1.0_ti_p-cosine*cosine)
                     e = factor * force * (1.0_ti_p+cosine)
                     deddt = -factor * force * sine
                  else if (angtyp(i) .eq. 'FOURIER') then
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
                  dea(1,ibloc) = dea(1,ibloc) + dedxib
                  dea(2,ibloc) = dea(2,ibloc) + dedyib
                  dea(3,ibloc) = dea(3,ibloc) + dedzib
c
                  dea(1,ialoc) = dea(1,ialoc) + dedxia
                  dea(2,ialoc) = dea(2,ialoc) + dedyia
                  dea(3,ialoc) = dea(3,ialoc) + dedzia
c
                  dea(1,icloc) = dea(1,icloc) + dedxic
                  dea(2,icloc) = dea(2,icloc) + dedyic
                  dea(3,icloc) = dea(3,icloc) + dedzic
c
c     increment the internal virial tensor components
c
                  vxx = xab*dedxia + xcb*dedxic
                  vyx = yab*dedxia + ycb*dedxic
                  vzx = zab*dedxia + zcb*dedxic
                  vyy = yab*dedyia + ycb*dedyic
                  vzy = zab*dedyia + zcb*dedyic
                  vzz = zab*dedzia + zcb*dedzic
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
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
                  call image (xad,yad,zad)
                  call image (xbd,ybd,zbd)
                  call image (xcd,ycd,zcd)
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
                  call image (xap,yap,zap)
                  call image (xcp,ycp,zcp)
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
                  dedxia = dedxia + dpdxia
                  dedyia = dedyia + dpdyia
                  dedzia = dedzia + dpdzia
                  dedxib = dedxip
                  dedyib = dedyip
                  dedzib = dedzip
                  dedxic = dedxic + dpdxic
                  dedyic = dedyic + dpdyic
                  dedzic = dedzic + dpdzic
                  dedxid = -dedxia - dedxib - dedxic
                  dedyid = -dedyia - dedyib - dedyic
                  dedzid = -dedzia - dedzib - dedzic
c
c     increment the total bond angle energy and derivatives
c
                  ea = ea + e
                  dea(1,ibloc) = dea(1,ibloc) + dedxib
                  dea(2,ibloc) = dea(2,ibloc) + dedyib
                  dea(3,ibloc) = dea(3,ibloc) + dedzib
c
                  dea(1,ialoc) = dea(1,ialoc) + dedxia
                  dea(2,ialoc) = dea(2,ialoc) + dedyia
                  dea(3,ialoc) = dea(3,ialoc) + dedzia
c
                  dea(1,icloc) = dea(1,icloc) + dedxic
                  dea(2,icloc) = dea(2,icloc) + dedyic
                  dea(3,icloc) = dea(3,icloc) + dedzic
c
                  dea(1,idloc) = dea(1,idloc) + dedxid
                  dea(2,idloc) = dea(2,idloc) + dedyid
                  dea(3,idloc) = dea(3,idloc) + dedzid
c
c     aMD storage if waters are considered
c
                  if (use_amd_wat1) then
                  if (type(ia) == aMDwattype(1) .or. type(ib)
     $            == aMDwattype(1) .or. type(ic) == aMDwattype(1))
     $            then
                     eW1aMD = eW1aMD + e
                     deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedxia
                     deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedyia
                     deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedzia
                     deW1aMD(1,ibloc) = deW1aMD(1,ibloc) + dedxib
                     deW1aMD(2,ibloc) = deW1aMD(2,ibloc) + dedyib
                     deW1aMD(3,ibloc) = deW1aMD(3,ibloc) + dedzib
                     deW1aMD(1,icloc) = deW1aMD(1,icloc) + dedxic
                     deW1aMD(2,icloc) = deW1aMD(2,icloc) + dedyic
                     deW1aMD(3,icloc) = deW1aMD(3,icloc) + dedzic
                  end if
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xad*dedxia + xbd*dedxib + xcd*dedxic
                  vyx = yad*dedxia + ybd*dedxib + ycd*dedxic
                  vzx = zad*dedxia + zbd*dedxib + zcd*dedxic
                  vyy = yad*dedyia + ybd*dedyib + ycd*dedyic
                  vzy = zad*dedyia + zbd*dedyib + zcd*dedyic
                  vzz = zad*dedzia + zbd*dedzib + zcd*dedzic
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
               end if
            end if
         end if
      end do
!$acc update device(dea,vir)
      return
      end
