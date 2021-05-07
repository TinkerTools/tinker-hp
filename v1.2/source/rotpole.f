c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole
      use atmlst
      use mpole
      implicit none
      integer i,iipole,iglob
      real*8 a(3,3)
c
c
c     rotate the atomic multipoles at each site in turn
c
      do i = 1, npolebloc
         iipole = poleglob(i)
         iglob = ipole(iipole)
         call rotmat (iipole,iglob,a)
         call rotsite (iipole,a)
      end do
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call rotmat (iipole,iglob,a)
        call rotsite (iipole,a)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotmat  --  find global frame rotation matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotmat" finds the rotation matrix that converts from the local
c     coordinate system to the global frame at a multipole site
c
c
      subroutine rotmat (iipole,iglob,a)
      use atoms
      use bound
      use mpole
      implicit none
      integer iipole,iglob
      integer ix,iy,iz
      real*8 r,dot
      real*8 random
      real*8 xi,yi,zi
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      real*8 a(3,3)
c
c      if (polsiz(iglob).eq.0) then
c        a = 0d0
c        return
c      end if
c
c     get coordinates and frame definition for the multipole site
c
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      ix = xaxis(iipole)
      iy = yaxis(iipole)
      iz = zaxis(iipole)
c
c     use the identity matrix as the default rotation matrix
c
      a(1,1) = 1.0d0
      a(2,1) = 0.0d0
      a(3,1) = 0.0d0
      a(1,3) = 0.0d0
      a(2,3) = 0.0d0
      a(3,3) = 1.0d0
c
c     Z-Only method rotation matrix elements for z-axis only
c
      if (polaxe(iipole) .eq. 'Z-Only') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = random ()
         dy = random ()
         dz = random ()
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     Z-then-X method rotation matrix elements for z- and x-axes
c
      else if (polaxe(iipole) .eq. 'Z-then-X') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         if (use_polymer) call image(dx,dy,dz)
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     Bisector method rotation matrix elements for z- and x-axes
c
      else if (polaxe(iipole) .eq. 'Bisector') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     Z-Bisect method rotation matrix elements for z- and x-axes
c
      else if (polaxe(iipole) .eq. 'Z-Bisect') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx = dx / r
         dy = dy / r
         dz = dz / r
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     3-Fold method rotation matrix elements for z- and x-axes
c
      else if (polaxe(iipole) .eq. '3-Fold') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         if (use_polymer) call image(dx,dy,dz)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx3 = dx / r
         dy3 = dy / r
         dz3 = dz / r
         dx = dx1 + dx2 + dx3
         dy = dy1 + dy2 + dy3
         dz = dz1 + dz2 + dz3
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotsite  --  rotate multipoles at single site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotsite" computes the atomic multipoles at a specified site
c     in the global coordinate frame by applying a rotation matrix
c
c
      subroutine rotsite (iipole,a)
      use mpole
      implicit none
      integer i,j,k,m
      integer iipole
      real*8 a(3,3)
      real*8 m2(3,3)
      real*8 r2(3,3)
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1,iipole) = pole(1,iipole)
c
c     rotate the dipoles to the global coordinate frame
c
      do i = 2, 4
         rpole(i,iipole) = 0.0d0
         do j = 2, 4
            rpole(i,iipole) = rpole(i,iipole)+pole(j,iipole)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles to the global coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = pole(k,iipole)
            r2(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               r2(i,j) = r2(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     r2(i,j) = r2(i,j) + a(i,k)*a(j,m)*m2(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 5
      do i = 1, 3
         do j = 1, 3
            rpole(k,iipole) = r2(i,j)
            k = k + 1
         end do
      end do
      return
      end
