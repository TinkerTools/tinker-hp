!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine eopdist1  --  out-of-plane dist energy & derivs  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "eopdist1" computes the out-of-plane distance potential
!     energy and first derivatives at trigonal centers via
!     the central atom height
!
!
#include "tinker_macro.h"
subroutine eopdist1
   use angpot
   use atmlst
   use atoms
   use bound
   use deriv
   use domdec
   use energi
   use group
   use opdist
   use tinheader ,only:ti_p,re_p
   use usage
   use virial
   implicit none
   integer i,ia,ib,ic,id,iopdist
   integer ialoc,ibloc,icloc,idloc
   real(t_p) e,force
   real(t_p) dot,deddt
   real(t_p) dt,dt2,dt3,dt4
   real(t_p) xia,yia,zia
   real(t_p) xib,yib,zib
   real(t_p) xic,yic,zic
   real(t_p) xid,yid,zid
   real(t_p) xad,yad,zad
   real(t_p) xbd,ybd,zbd
   real(t_p) xcd,ycd,zcd
   real(t_p) xt,yt,zt
   real(t_p) rt2,drt2
   real(t_p) xtd,ytd,ztd
   real(t_p) dedxia,dedyia,dedzia
   real(t_p) dedxib,dedyib,dedzib
   real(t_p) dedxic,dedyic,dedzic
   real(t_p) dedxid,dedyid,dedzid
   real(t_p) vxx,vyy,vzz
   real(t_p) vyx,vzx,vzy
   real(t_p) fgrp
   logical proceed
!
!$acc update host(deopb,deopd,vir)
!
!     zero out out-of-plane energy and first derivatives
!
   eopd = 0.0_ti_p
!
!     calculate the out-of-plane distance energy and derivatives
!
   do iopdist = 1, nopdistloc
      i = opdistglob(iopdist)
      ia = iopd(1,i)
      ib = iopd(2,i)
      ic = iopd(3,i)
      id = iopd(4,i)
      ialoc = loc(ia)
      ibloc = loc(ib)
      icloc = loc(ic)
      idloc = loc(id)
      force = opdk(i)
!
!     decide whether to compute the current interaction
!
      if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
      proceed = (use(ia) .or. use(ib) .or.&
         &use(ic) .or. use(id))
!
!     get the coordinates of the defining atoms
!
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
!
!     compute the out-of-plane distance for central atom
!
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
         xt = ybd*zcd - zbd*ycd
         yt = zbd*xcd - xbd*zcd
         zt = xbd*ycd - ybd*xcd
         rt2 = xt*xt + yt*yt + zt*zt
         dot = xt*xad + yt*yad + zt*zad
         drt2 = dot / rt2
         dt2 = dot * drt2
         dt = sqrt(dt2)
         dt3 = dt2 * dt
         dt4 = dt2 * dt2
!
!     find the out-of-plane energy and master chain rule terms
!
         e = opdunit * force * dt2&
            &* (1.0_ti_p+copd*dt+qopd*dt2+popd*dt3+sopd*dt4)
         deddt = opdunit * force * drt2&
            &* (2.0_ti_p + 3.0_ti_p*copd*dt + 4.0_ti_p*qopd*dt2&
            &+ 5.0_ti_p*popd*dt3 + 6.0_ti_p*sopd*dt4)

         if(use_group) then
            call groups(fgrp,ia,ib,ic,id,0,0)
            e = e*fgrp
            deddt = deddt*fgrp
         endif
!
!     scale the interaction based on its group membership
!
         if (use_group) then
            e     =  e    * fgrp
            deddt = deddt * fgrp
         end if
!
!     chain rule terms for first derivative components
!
         xtd = xad - xt*drt2
         ytd = yad - yt*drt2
         ztd = zad - zt*drt2
!
!     compute derivative components for this interaction
!
         dedxia = deddt * xt
         dedyia = deddt * yt
         dedzia = deddt * zt
         dedxib = deddt * (ycd*ztd-zcd*ytd)
         dedyib = deddt * (zcd*xtd-xcd*ztd)
         dedzib = deddt * (xcd*ytd-ycd*xtd)
         dedxic = deddt * (zbd*ytd-ybd*ztd)
         dedyic = deddt * (xbd*ztd-zbd*xtd)
         dedzic = deddt * (ybd*xtd-xbd*ytd)
!
!     get some derivative components by difference
!
         dedxid = -dedxia - dedxib - dedxic
         dedyid = -dedyia - dedyib - dedyic
         dedzid = -dedzia - dedzib - dedzic
!
!     increment the out-of-plane distance energy and gradient
!
         eopd = eopd + e
         deopd(1,ia) = deopd(1,ia) + dedxia
         deopd(2,ia) = deopd(2,ia) + dedyia
         deopd(3,ia) = deopd(3,ia) + dedzia
!
         deopb(1,ialoc) = deopb(1,ialoc) + dedxia
         deopb(2,ialoc) = deopb(2,ialoc) + dedyia
         deopb(3,ialoc) = deopb(3,ialoc) + dedzia
!
         deopb(1,ibloc) = deopb(1,ibloc) + dedxib
         deopb(2,ibloc) = deopb(2,ibloc) + dedyib
         deopb(3,ibloc) = deopb(3,ibloc) + dedzib
!
         deopb(1,icloc) = deopb(1,icloc) + dedxic
         deopb(2,icloc) = deopb(2,icloc) + dedyic
         deopb(3,icloc) = deopb(3,icloc) + dedzic
!
         deopb(1,idloc) = deopb(1,idloc) + dedxid
         deopb(2,idloc) = deopb(2,idloc) + dedyid
         deopb(3,idloc) = deopb(3,idloc) + dedzid
!
!     increment the internal virial tensor components
!
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
   end do
!$acc update device(deopb,deopd,eopd,vir)
   return
end
