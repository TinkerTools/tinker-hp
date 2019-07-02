c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdrest  --  stop system translation & rotation  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdrest" finds and removes any translational or rotational
c     kinetic energy of the overall system center of mass
c
c
      subroutine mdrest (istep)
      use atmtyp
      use atoms
      use bound
      use domdec
      use group
      use inform
      use iounit
      use mdstuf
      use moldyn
      use units
      use mpi
      implicit none
      integer i,j,k,istep,iglob,ierr
      real*8 etrans,erot
      real*8 weigh,totmass,eps
      real*8 xx,yy,zz,xy,xz,yz
      real*8 xtot,ytot,ztot
      real*8 xdel,ydel,zdel
      real*8 mang(3),vang(3)
      real*8 vtot(3),tensor(3,3)
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)
c
c
c     check steps between center of mass motion removal
c
      if (.not.dorest)  return
      if (mod(istep,irest) .ne. 0)  return
c
c     zero out the total mass and overall linear velocity
c
      totmass = 0.0d0
      do j = 1, 3
         vtot(j) = 0.0d0
      end do
c
c     compute linear velocity of the system center of mass
c
      do i = 1, nloc
         iglob = glob(i)
         weigh = mass(iglob)
         totmass = totmass + weigh
         do j = 1, 3
            vtot(j) = vtot(j) + v(j,iglob)*weigh
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,vtot,3,MPI_REAL8,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,totmass,1,MPI_REAL8,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)
c
c     compute translational kinetic energy of overall system
c
      etrans = 0.0d0
      do j = 1, 3
         vtot(j) = vtot(j) / totmass
         etrans = etrans + vtot(j)**2
      end do
      etrans = 0.5d0 * etrans * totmass / convert
c
c     find the center of mass coordinates of the overall system
c
      if (.not. use_bounds) then
         xtot = 0.0d0
         ytot = 0.0d0
         ztot = 0.0d0
            do i = 1, nloc
               iglob = glob(i)
               weigh = mass(iglob)
               xtot = xtot + x(iglob)*weigh
               ytot = ytot + y(iglob)*weigh
               ztot = ztot + z(iglob)*weigh
            end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,xtot,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ytot,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ztot,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         xtot = xtot / totmass
         ytot = ytot / totmass
         ztot = ztot / totmass
c
c     compute the angular momentum of the overall system
c
         do j = 1, 3
            mang(j) = 0.0d0
         end do
         do i = 1, nloc
               iglob = glob(i)
               weigh = mass(iglob)
               mang(1) = mang(1) + (y(iglob)*v(3,iglob)-z(iglob)*
     $          v(2,iglob))*weigh
               mang(2) = mang(2) + (z(iglob)*v(1,iglob)-x(iglob)*
     $          v(3,iglob))*weigh
               mang(3) = mang(3) + (x(iglob)*v(2,iglob)-y(iglob)*
     $          v(1,iglob))*weigh
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,mang,3,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
c
         mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
         mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
         mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass
c
c     calculate the moment of inertia tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         do i = 1, nloc
            iglob = glob(i)
            weigh = mass(iglob)
            xdel = x(iglob) - xtot
            ydel = y(iglob) - ytot
            zdel = z(iglob) - ztot
            xx = xx + xdel*xdel*weigh
            xy = xy + xdel*ydel*weigh
            xz = xz + xdel*zdel*weigh
            yy = yy + ydel*ydel*weigh
            yz = yz + ydel*zdel*weigh
            zz = zz + zdel*zdel*weigh
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,xx,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xy,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xz,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,yy,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,yz,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zz,1,MPI_REAL8,MPI_SUM,
     $      MPI_COMM_WORLD,ierr)
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
c
c     fix to avoid singularity for one- or two-body systems
c
         if (n .le. 2) then
            eps = 0.000001d0
            tensor(1,1) = tensor(1,1) + eps
            tensor(2,2) = tensor(2,2) + eps
            tensor(3,3) = tensor(3,3) + eps
         end if
c
c     diagonalize the moment of inertia tensor
c
         call invert (3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
         erot = 0.0d0
         do i = 1, 3
            vang(i) = 0.0d0
            do j = 1, 3
               vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
         end do
         erot = 0.5d0 * erot / convert
      end if
c
c     eliminate any translation of the overall system
c
      do i = 1, nloc
         iglob = glob(i)
         do j = 1, 3
            v(j,iglob) = v(j,iglob) - vtot(j)
         end do
      end do
c
c     print the translational velocity of the overall system
c
      if (debug) then
         write (iout,10)  (vtot(i),i=1,3),etrans
   10    format (' System Linear Velocity :  ',3d12.2,
     &           /,' Translational Kinetic Energy :',10x,f12.4,
     &              ' Kcal/mole')
      end if
      return
      end
