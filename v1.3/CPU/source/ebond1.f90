!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine ebond1  --  bond stretch energy & derivatives  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "ebond1" calculates the bond stretching energy and
!     first derivatives with respect to Cartesian coordinates
!
!
!> @brief 
!> calculates the bond stretching energy and
!> first derivatives with respect to Cartesian coordinates
!> @param no params
subroutine ebond1
   use atmlst
   use atoms
   use bndpot
   use bond
   use bound
   use deriv
   use domdec
   use energi
   use group
   use iounit
   use inform
   use usage
   use virial
   implicit none
   integer i,ia,ib,ialoc,ibloc
   integer ibond
   real*8 e,ideal,force,ba2
   real*8 expterm,bde
   real*8 dt,dt2,deddt
   real*8 de,dedx,dedy,dedz
   real*8 xab,yab,zab,rab
   real*8 vxx,vyy,vzz
   real*8 vyx,vzx,vzy
   real*8 fgrp
   logical proceed
!
   if (deb_Path) write(iout,*), 'ebond1 '
!
!
!     zero out the bond energy and first derivatives
!
   eb = 0.0d0
!
!     calculate the bond stretch energy and first derivatives
!
   do ibond = 1, nbondloc
      i = bndglob(ibond)
      ia = ibnd(1,i)
      ib = ibnd(2,i)
      ialoc = loc(ia)
      ibloc = loc(ib)
      ideal = bl(i)
      force = bk(i)
!
!     decide whether to compute the current interaction
!
      proceed = (use(ia) .or. use(ib))
      if (use_group) call groups(fgrp,ia,ib,0,0,0,0)
!
!     compute the value of the bond length deviation
!
      if (proceed) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         if (use_polymer)  call image (xab,yab,zab)
         rab = sqrt(xab*xab + yab*yab + zab*zab)
         dt = rab - ideal
!
!     harmonic potential uses Taylor expansion of Morse potential
!     through the fourth power of the bond length deviation
!
         if (bndtyp(i) .eq. 'HARMONIC') then
            dt2 = dt * dt
            e = bndunit * force * dt2 * (1.0d0+cbnd*dt+qbnd*dt2)
            deddt = 2.0d0 * bndunit * force * dt&
            &* (1.0d0+1.5d0*cbnd*dt+2.0d0*qbnd*dt2)
!
!     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
!     with the approximations alpha = sqrt(ForceConst/BDE) = -2
!     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
!
         else if (bndtyp(i) .eq. 'MORSE') then
            expterm = exp(-ba(i)*dt)
            bde = bndunit * force / ba(i)**2
            e = bde * (1.0d0-expterm)**2
            deddt = 2.0d0 * bde * ba(i) * (1.0d0-expterm) * expterm
!
!     Morse potential expanded to 4th order
!
         else if (bndtyp(i) .eq. 'MORSE4') then
            dt2 = dt * dt
            ba2 = 7.d0/12.d0*ba(i) * ba(i)
            e = bndunit * force * dt2 * (1.0d0-ba(i)*dt&
            &+ ba2*dt2)
            deddt = bndunit * force * dt * (2.0d0-3.d0*ba(i)*dt&
            &+ 4.d0*ba2*dt2)
         end if
!
!     scale the interaction based on its group membership
!
         if (use_group) then
            e = e * fgrp
            deddt = deddt * fgrp
         end if
!
!     compute chain rule terms needed for derivatives
!
         if (rab .eq. 0.0d0) then
            de = 0.0d0
         else
            de = deddt / rab
         end if
         dedx = de * xab
         dedy = de * yab
         dedz = de * zab
!
!     increment the total bond energy and first derivatives
!
         eb = eb + e
         deb(1,ialoc) = deb(1,ialoc) + dedx
         deb(2,ialoc) = deb(2,ialoc) + dedy
         deb(3,ialoc) = deb(3,ialoc) + dedz
!
         deb(1,ibloc) = deb(1,ibloc) - dedx
         deb(2,ibloc) = deb(2,ibloc) - dedy
         deb(3,ibloc) = deb(3,ibloc) - dedz
!
!     increment the internal virial tensor components
!
         vxx = xab * dedx
         vyx = yab * dedx
         vzx = zab * dedx
         vyy = yab * dedy
         vzy = zab * dedy
         vzz = zab * dedz
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
