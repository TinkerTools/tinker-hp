c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eopdist1  --  out-of-plane dist energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eopdist1" computes the out-of-plane distance potential
c     energy and first derivatives at trigonal centers via
c     the central atom height
c
c
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
      use usage
      use virial
      implicit none
      integer i,ia,ib,ic,id,iopdist
      integer ialoc,ibloc,icloc,idloc
      real*8 e,force
      real*8 dot,deddt
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xt,yt,zt
      real*8 rt2,drt2
      real*8 xtd,ytd,ztd
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fgrp
      logical proceed
c
c
c     zero out out-of-plane energy and first derivatives
c
      eopd = 0.0d0
c
c     calculate the out-of-plane distance energy and derivatives
c
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
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
         proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the defining atoms
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
c
c     compute the out-of-plane distance for central atom
c
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
c
c     find the out-of-plane energy and master chain rule terms
c
            e = opdunit * force * dt2
     &             * (1.0d0+copd*dt+qopd*dt2+popd*dt3+sopd*dt4)
            deddt = opdunit * force * drt2
     &                 * (2.0d0 + 3.0d0*copd*dt + 4.0d0*qopd*dt2
     &                     + 5.0d0*popd*dt3 + 6.0d0*sopd*dt4)
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     chain rule terms for first derivative components
c
            xtd = xad - xt*drt2
            ytd = yad - yt*drt2
            ztd = zad - zt*drt2
c
c     compute derivative components for this interaction
c
            dedxia = deddt * xt
            dedyia = deddt * yt
            dedzia = deddt * zt
            dedxib = deddt * (ycd*ztd-zcd*ytd)
            dedyib = deddt * (zcd*xtd-xcd*ztd)
            dedzib = deddt * (xcd*ytd-ycd*xtd)
            dedxic = deddt * (zbd*ytd-ybd*ztd)
            dedyic = deddt * (xbd*ztd-zbd*xtd)
            dedzic = deddt * (ybd*xtd-xbd*ytd)
c
c     get some derivative components by difference
c
            dedxid = -dedxia - dedxib - dedxic
            dedyid = -dedyia - dedyib - dedyic
            dedzid = -dedzia - dedzib - dedzic
c
c     increment the out-of-plane distance energy and gradient
c
            eopd = eopd + e
            deopd(1,ia) = deopd(1,ia) + dedxia
            deopd(2,ia) = deopd(2,ia) + dedyia
            deopd(3,ia) = deopd(3,ia) + dedzia
c
            deopb(1,ialoc) = deopb(1,ialoc) + dedxia
            deopb(2,ialoc) = deopb(2,ialoc) + dedyia
            deopb(3,ialoc) = deopb(3,ialoc) + dedzia
c
            deopb(1,ibloc) = deopb(1,ibloc) + dedxib
            deopb(2,ibloc) = deopb(2,ibloc) + dedyib
            deopb(3,ibloc) = deopb(3,ibloc) + dedzib
c
            deopb(1,icloc) = deopb(1,icloc) + dedxic
            deopb(2,icloc) = deopb(2,icloc) + dedyic
            deopb(3,icloc) = deopb(3,icloc) + dedzic
c
            deopb(1,idloc) = deopb(1,idloc) + dedxid
            deopb(2,idloc) = deopb(2,idloc) + dedyid
            deopb(3,idloc) = deopb(3,idloc) + dedzid
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
      end do
      return
      end
