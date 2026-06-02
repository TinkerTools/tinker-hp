!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine ebond  --  bond stretch potential energy  ##
!     ##                                                       ##
!     ###########################################################
!
!
!     "ebond" calculates the bond stretching energy
!
!
!> @brief 
!> calculates the bond stretching energy
!> @param no params
subroutine ebond
   use atmlst
   use atoms
   use bndpot
   use bond
   use bound
   use energi
   use group
   use inform
   use iounit
   use usage
   implicit none
   integer i,ia,ib,ibond
   real*8 e,ideal,force
   real*8 expterm,bde
   real*8 dt,dt2
   real*8 xab,yab,zab,rab
   real*8 fgrp
   logical proceed
!
   if (deb_Path) write(iout,*), 'ebond '
!
!
!     zero out the bond stretching energy
!
   eb = 0.0d0
!
!     calculate the bond stretching energy term
!
   do ibond = 1, nbondloc
      i = bndglob(ibond)
      ia = ibnd(1,i)
      ib = ibnd(2,i)
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
!
!     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
!     with the approximations alpha = sqrt(ForceConst/BDE) = -2
!     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
!
         else if (bndtyp(i) .eq. 'MORSE') then
            expterm = exp(-ba(i)*dt)
            bde = bndunit * force / ba(i)**2
            e = bde * (1.0d0-expterm)**2
!
!     Morse potential expanded to 4th order
!
         else if (bndtyp(i) .eq. 'MORSE4') then
            dt2 = dt * dt
            e = bndunit * force * dt2 * (1.0d0-ba(i)*dt&
            &+ 7.d0/12.d0*ba(i)*ba(i)*dt2)
         end if
!
!     scale the interaction based on its group membership
!
         if (use_group)  e = e * fgrp
!
!     increment the total bond stretching energy
!
         eb = eb + e
      end if
   end do
   return
end
