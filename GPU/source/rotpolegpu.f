c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
#include "tinker_precision.h"
      subroutine rotpolegpu
      use atmlst
      use mpole
      use potent ,only: use_pmecore
      use utilgpu,only: dir_queue,rec_queue,def_queue
      implicit none
      integer i,iipole,iglob
      real(t_p) a(3,3)
      interface
       subroutine rot_mat_site(nk,poleglobvec)
         integer,intent(in)::nk
         integer,intent(in)::poleglobvec(:)
       end subroutine
      end interface
c
c     rotate the atomic multipoles at each site in turn
c
      ! FIXME
      ! There is a bug here on device
      ! With 2 devices on dhfr2 withtout restart and respa1 integrator
      ! at 8 fs of timestep in double precision
      ! device 0 loses polerecglob adress when entering rot_mat_site
      ! when deleting second call the it is poleglob's adress which is
      ! lost at runtime

      def_queue=rec_queue
      call rot_mat_site(npolebloc,poleglob)
      if (npolebloc.ne.npole.or.use_pmecore)
     &   call rot_mat_site(npolerecloc,polerecglob)
      def_queue=dir_queue

      end
c
c     "rot_mat_site" finds the rotation matrix that converts from 
c     the local coordinate system to the global frame at a multipole site
c
c     "rotsite" computes the atomic multipoles at a specified site
c     in the global coordinate frame by applying a rotation matrix
c
      subroutine rot_mat_site(nk,poleglobvec)
      use atoms
      use mpole
      use random_mod
      use utilgpu
      implicit none
      integer,intent(in):: nk
      integer,intent(in):: poleglobvec(:)
      integer ii,iipole,iglob
      integer ix,iy,iz
      integer i,j,k,m
      integer::kk=0,ksave
      integer axetyp
      real(t_p) r,dot,invr
      real(t_p) xi,yi,zi
      real(t_p) dx,dy,dz
      real(t_p) dx1,dy1,dz1
      real(t_p) dx2,dy2,dz2
      real(t_p) dx3,dy3,dz3
      real(t_p) rpole_i_iglob
      real(t_p) a(3,3)
      real(t_p) m2(3,3)
      real(t_p) r2(3,3)
c
c     get coordinates and frame definition for the multipole site
c

#ifdef _OPENACC
!$acc data present(samplevec)
      if (host_rand_platform) then
         do ii = 1, 3*nZ_Onlyloc
            samplevec(ii) = random ()
         end do
!$acc update device(samplevec(1:3*nZ_Onlyloc)) async(def_queue)
      else
         if (nZ_Onlyloc.gt.0) then
            call randomgpu(samplevec(1),3*nZ_Onlyloc)
         end if
      end if
!$acc end data
#endif
!$acc parallel loop private(a,m2,r2) copyin(kk)
!$acc&         present(poleglobvec,ipole,ipolaxe,x,y,z,rpole,
!$acc&   pole,xaxis,yaxis,zaxis)
!$acc&         async(def_queue)
      do ii = 1, nk
         iipole = poleglobvec(ii)
         iglob  = ipole(iipole)
         axetyp = ipolaxe(iglob) 
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)
         ix     = xaxis(iipole)
         iy     = yaxis(iipole)
         iz     = zaxis(iipole)
c
c     use the identity matrix as the default rotation matrix
c
         a(1,1) = 1.0_ti_p
         a(2,1) = 0.0_ti_p
         a(3,1) = 0.0_ti_p
         a(1,3) = 0.0_ti_p
         a(2,3) = 0.0_ti_p
         a(3,3) = 1.0_ti_p
c        
c     Z-Only method rotation matrix elements for z-axis only
c        
         if (btest(axetyp,3)) then  !'Z-Only'
            dx  = x(iz) - xi
            dy  = y(iz) - yi
            dz  = z(iz) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,3) = dx * invr
            a(2,3) = dy * invr
            a(3,3) = dz * invr
#ifdef _OPENACC
!$acc atomic capture
            ksave = kk
            kk    = kk + 1
!$acc end atomic
            dx  = samplevec(3*ksave+1)
            dy  = samplevec(3*ksave+2)
            dz  = samplevec(3*ksave+3)
#else
            dx  = random ()
            dy  = random ()
            dz  = random ()
#endif
            dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
            dx  = dx - dot*a(1,3)
            dy  = dy - dot*a(2,3)
            dz  = dz - dot*a(3,3)
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,1) = dx * invr
            a(2,1) = dy * invr
            a(3,1) = dz * invr
c
c     Z-then-X method rotation matrix elements for z- and x-axes
c
         else if (btest(axetyp,4)) then  !'Z-then-X'
            dx  = x(iz) - xi
            dy  = y(iz) - yi
            dz  = z(iz) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,3) = dx * invr
            a(2,3) = dy * invr
            a(3,3) = dz * invr
            dx  = x(ix) - xi
            dy  = y(ix) - yi
            dz  = z(ix) - zi
            dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
            dx  = dx - dot*a(1,3)
            dy  = dy - dot*a(2,3)
            dz  = dz - dot*a(3,3)
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,1) = dx * invr
            a(2,1) = dy * invr
            a(3,1) = dz * invr
c
c     Bisector method rotation matrix elements for z- and x-axes
c
         else if (btest(axetyp,1)) then  !'Bisector'
            dx  = x(iz) - xi
            dy  = y(iz) - yi
            dz  = z(iz) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx1 = dx * invr
            dy1 = dy * invr
            dz1 = dz * invr
            dx  = x(ix) - xi
            dy  = y(ix) - yi
            dz  = z(ix) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx2 = dx * invr
            dy2 = dy * invr
            dz2 = dz * invr
            dx  = dx1 + dx2
            dy  = dy1 + dy2
            dz  = dz1 + dz2
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,3) = dx * invr
            a(2,3) = dy * invr
            a(3,3) = dz * invr
            dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
            dx  = dx2 - dot*a(1,3)
            dy  = dy2 - dot*a(2,3)
            dz  = dz2 - dot*a(3,3)
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,1) = dx * invr
            a(2,1) = dy * invr
            a(3,1) = dz * invr
c
c     Z-Bisect method rotation matrix elements for z- and x-axes
c
         else if (btest(axetyp,2)) then  !'Z-Bisect'
            dx  = x(iz) - xi
            dy  = y(iz) - yi
            dz  = z(iz) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,3) = dx * invr
            a(2,3) = dy * invr
            a(3,3) = dz * invr
            dx  = x(ix) - xi
            dy  = y(ix) - yi
            dz  = z(ix) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx1 = dx * invr
            dy1 = dy * invr
            dz1 = dz * invr
            dx  = x(iy) - xi
            dy  = y(iy) - yi
            dz  = z(iy) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx2 = dx * invr
            dy2 = dy * invr
            dz2 = dz * invr
            dx  = dx1 + dx2
            dy  = dy1 + dy2
            dz  = dz1 + dz2
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx  = dx * invr
            dy  = dy * invr
            dz  = dz * invr
            dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
            dx  = dx - dot*a(1,3)
            dy  = dy - dot*a(2,3)
            dz  = dz - dot*a(3,3)
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,1) = dx * invr
            a(2,1) = dy * invr
            a(3,1) = dz * invr
c        
c     3-Fold method rotation matrix elements for z- and x-axes
c        
         else if (btest(axetyp,0)) then  !'3-Fold'
            dx  = x(iz) - xi
            dy  = y(iz) - yi
            dz  = z(iz) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx1 = dx * invr
            dy1 = dy * invr
            dz1 = dz * invr
            dx  = x(ix) - xi
            dy  = y(ix) - yi
            dz  = z(ix) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx2 = dx * invr
            dy2 = dy * invr
            dz2 = dz * invr
            dx  = x(iy) - xi
            dy  = y(iy) - yi
            dz  = z(iy) - zi
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            dx3 = dx * invr
            dy3 = dy * invr
            dz3 = dz * invr
            dx  = dx1 + dx2 + dx3
            dy  = dy1 + dy2 + dy3
            dz  = dz1 + dz2 + dz3
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,3) = dx * invr
            a(2,3) = dy * invr
            a(3,3) = dz * invr
            dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
            dx  = dx2 - dot*a(1,3)
            dy  = dy2 - dot*a(2,3)
            dz  = dz2 - dot*a(3,3)
            r   = sqrt(dx*dx + dy*dy + dz*dz)
            invr= 1.0_ti_p/r
            a(1,1) = dx * invr
            a(2,1) = dy * invr
            a(3,1) = dz * invr
         end if
c        
c     finally, find rotation matrix elements for the y-axis
c        
         a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
         a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
         a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
c
c
c     "rotsite"
c
c     monopoles have the same value in any coordinate frame
c
c
         rpole(1,iglob) = pole(1,iglob)
c
c     rotate the dipoles to the global coordinate frame
c
!$acc    loop seq
         do i = 2, 4
            rpole_i_iglob = 0.0_ti_p
!$acc       loop seq
            do j = 2, 4
               rpole_i_iglob = rpole_i_iglob+pole(j,iglob)*a(i-1,j-1)
            end do
            rpole(i,iglob) = rpole_i_iglob
         end do
c
c     rotate the quadrupoles to the global coordinate frame
c
         k = 5
!$acc    loop seq collapse(2)
         do i = 1, 3
            do j = 1, 3
               m2(i,j) = pole(k,iglob)
               r2(i,j) = 0.0_ti_p
               k = k + 1
            end do
         end do
!$acc    loop seq collapse(2)
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
!$acc    loop seq collapse(2)
         do i = 1, 3
            do j = 1, 3
               rpole(k,iglob) = r2(i,j)
               k = k + 1
            end do
         end do
      end do

      end subroutine
