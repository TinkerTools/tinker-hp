c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eurey  --  Urey-Bradley potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eurey" calculates the Urey-Bradley 1-3 interaction energy
c
c
      subroutine eurey
      use atmlst
      use atoms
      use bound
      use energi
      use group
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
c
c
c     zero out the Urey-Bradley interaction energy
c
      eub = 0.0d0
c
c     calculate the Urey-Bradley 1-3 energy term
c
      do iurey = 1, nureyloc
         i = ureyglob(iurey)
         ia = iury(1,i)
         ic = iury(3,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ic,0,0,0,0)
         proceed = (use(ia) .or. use(ic))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xac = x(ia) - x(ic)
            yac = y(ia) - y(ic)
            zac = z(ia) - z(ic)
            if (use_polymer)  call image (xac,yac,zac)
            rac = sqrt(xac*xac + yac*yac + zac*zac)
c
c     calculate the Urey-Bradley energy for this interaction
c
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
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total Urey-Bradley energy
c
            eub = eub + e
         end if
      end do
      return
      end
