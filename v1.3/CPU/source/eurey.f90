!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine eurey  --  Urey-Bradley potential energy  ##
!     ##                                                       ##
!     ###########################################################
!
!
!     "eurey" calculates the Urey-Bradley 1-3 interaction energy
!
!
!> @brief 
!> calculates the Urey-Bradley 1-3 interaction energy
!> @param no params
subroutine eurey
   use atmlst
   use atoms
   use bound
   use energi
   use group
   use inform
   use iounit
   use urey
   use urypot
   use usage
   implicit none
   integer i,ia,ic,iurey
   real*8 e,ideal,force
   real*8 dt,dt2
   real*8 xac,yac,zac,rac
   real*8 fgrp
   logical proceed
!
   if (deb_Path) write(iout,*), 'eurey '
!
!
!     zero out the Urey-Bradley interaction energy
!
   eub = 0.0d0
!
!     calculate the Urey-Bradley 1-3 energy term
!
   do iurey = 1, nureyloc
      i = ureyglob(iurey)
      ia = iury(1,i)
      ic = iury(3,i)
      ideal = ul(i)
      force = uk(i)
!
!     decide whether to compute the current interaction
!
      if (use_group)  call groups (fgrp,ia,ic,0,0,0,0)
      proceed = (use(ia) .or. use(ic))
!
!     compute the value of the 1-3 distance deviation
!
      if (proceed) then
         xac = x(ia) - x(ic)
         yac = y(ia) - y(ic)
         zac = z(ia) - z(ic)
         if (use_polymer)  call image (xac,yac,zac)
         rac = sqrt(xac*xac + yac*yac + zac*zac)
!
!     calculate the Urey-Bradley energy for this interaction
!
         if (ureytyp(i) == 'ANGREP') then
            e = ureyunit * force * exp(-rac/ideal)
         elseif (ureytyp(i) == 'UREYQUAR') then
            dt  = ideal / rac
            dt2 = dt * dt
            e   = ureyunit *force * (dt2 - 1.0d0)**2
         else
            dt = rac - ideal
            dt2 = dt * dt
            e = ureyunit * force * dt2 * (1.0d0+cury*dt+qury*dt2)
         endif
!
!     scale the interaction based on its group membership
!
         if (use_group)  e = e * fgrp
!
!     increment the total Urey-Bradley energy
!
         eub = eub + e
      end if
   end do
   return
end
