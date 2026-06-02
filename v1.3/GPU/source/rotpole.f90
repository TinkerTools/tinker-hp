!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "rotpole" constructs the set of atomic multipoles in the global
!     frame by applying the correct rotation matrix for each site
!
!
#include "tinker_precision.h"
subroutine rotpole
   use atmlst
   use mpole
   implicit none
   integer i,iipole,iglob
   real(t_p) a(3,3)
!
!
!     rotate the atomic multipoles at each site in turn
!
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
!$acc update device(rpole(:,:))
   return
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine rotmat  --  find global frame rotation matrix  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "rotmat" finds the rotation matrix that converts from the local
!     coordinate system to the global frame at a multipole site
!
!
subroutine rotmat (iipole,iglob,a)
   use atoms
   use mpole
   use random_mod
   use tinheader
   implicit none
   integer iipole,iglob
   integer ix,iy,iz
   real(t_p) r,dot
   real(t_p) xi,yi,zi
   real(t_p) dx,dy,dz
   real(t_p) dx1,dy1,dz1
   real(t_p) dx2,dy2,dz2
   real(t_p) dx3,dy3,dz3
   real(t_p) a(3,3)
!
!
!     get coordinates and frame definition for the multipole site
!
   xi = x(iglob)
   yi = y(iglob)
   zi = z(iglob)
   ix = xaxis(iipole)
   iy = yaxis(iipole)
   iz = zaxis(iipole)
!
!     use the identity matrix as the default rotation matrix
!
   a(1,1) = 1.0_ti_p
   a(2,1) = 0.0_ti_p
   a(3,1) = 0.0_ti_p
   a(1,3) = 0.0_ti_p
   a(2,3) = 0.0_ti_p
   a(3,3) = 1.0_ti_p
!
!     Z-Only method rotation matrix elements for z-axis only
!
   if (polaxe(iipole) .eq. 'Z-Only') then
      dx = x(iz) - xi
      dy = y(iz) - yi
      dz = z(iz) - zi
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
!
!     Z-then-X method rotation matrix elements for z- and x-axes
!
   else if (polaxe(iipole) .eq. 'Z-then-X') then
      dx = x(iz) - xi
      dy = y(iz) - yi
      dz = z(iz) - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      a(1,3) = dx / r
      a(2,3) = dy / r
      a(3,3) = dz / r
      dx = x(ix) - xi
      dy = y(ix) - yi
      dz = z(ix) - zi
      dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
      dx = dx - dot*a(1,3)
      dy = dy - dot*a(2,3)
      dz = dz - dot*a(3,3)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      a(1,1) = dx / r
      a(2,1) = dy / r
      a(3,1) = dz / r
!
!     Bisector method rotation matrix elements for z- and x-axes
!
   else if (polaxe(iipole) .eq. 'Bisector') then
      dx = x(iz) - xi
      dy = y(iz) - yi
      dz = z(iz) - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dx1 = dx / r
      dy1 = dy / r
      dz1 = dz / r
      dx = x(ix) - xi
      dy = y(ix) - yi
      dz = z(ix) - zi
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
!
!     Z-Bisect method rotation matrix elements for z- and x-axes
!
   else if (polaxe(iipole) .eq. 'Z-Bisect') then
      dx = x(iz) - xi
      dy = y(iz) - yi
      dz = z(iz) - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      a(1,3) = dx / r
      a(2,3) = dy / r
      a(3,3) = dz / r
      dx = x(ix) - xi
      dy = y(ix) - yi
      dz = z(ix) - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dx1 = dx / r
      dy1 = dy / r
      dz1 = dz / r
      dx = x(iy) - xi
      dy = y(iy) - yi
      dz = z(iy) - zi
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
!
!     3-Fold method rotation matrix elements for z- and x-axes
!
   else if (polaxe(iipole) .eq. '3-Fold') then
      dx = x(iz) - xi
      dy = y(iz) - yi
      dz = z(iz) - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dx1 = dx / r
      dy1 = dy / r
      dz1 = dz / r
      dx = x(ix) - xi
      dy = y(ix) - yi
      dz = z(ix) - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dx2 = dx / r
      dy2 = dy / r
      dz2 = dz / r
      dx = x(iy) - xi
      dy = y(iy) - yi
      dz = z(iy) - zi
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
!
!     finally, find rotation matrix elements for the y-axis
!
   a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
   a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
   a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
   return
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine rotsite  --  rotate multipoles at single site  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "rotsite" computes the atomic multipoles at a specified site
!     in the global coordinate frame by applying a rotation matrix
!
!
subroutine rotsite (iipole,a)
   use mpole
   use tinheader
   implicit none
   integer i,j,k,m
   integer iipole
   real(t_p) a(3,3)
   real(t_p) m2(3,3)
   real(t_p) r2(3,3)
!
!     monopoles have the same value in any coordinate frame
!
   rpole(1,iipole) = pole(1,iipole)
!
!     rotate the dipoles to the global coordinate frame
!
   do i = 2, 4
      rpole(i,iipole) = 0.0_ti_p
      do j = 2, 4
         rpole(i,iipole) = rpole(i,iipole)+pole(j,iipole)*a(i-1,j-1)
      end do
   end do
!
!     rotate the quadrupoles to the global coordinate frame
!
   k = 5
   do i = 1, 3
      do j = 1, 3
         m2(i,j) = pole(k,iipole)
         r2(i,j) = 0.0_ti_p
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
