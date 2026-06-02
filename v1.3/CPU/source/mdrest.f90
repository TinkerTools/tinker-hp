!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine mdrest  --  stop system translation & rotation  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "mdrest" finds and removes any translational or rotational
!     kinetic energy of the overall system center of mass
!
!
!> @brief 
!> finds and removes any translational or rotational
!> kinetic energy of the overall system center of mass
!> @param[in]  istep: index of timestep
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
   integer i,j,k,m,istep,iglob,ierr
   real*8 weigh,eps
   real*8 xx,yy,zz,xy,xz,yz
   real*8 xdel,ydel,zdel
   real*8 mang(3)
   real*8 tensor(3,3)
   real*8, allocatable :: totmass(:)
   real*8, allocatable :: etrans(:)
   real*8, allocatable :: erot(:)
   real*8, allocatable :: xtot(:)
   real*8, allocatable :: ytot(:)
   real*8, allocatable :: ztot(:)
   real*8, allocatable :: vtot(:,:)
   real*8, allocatable :: vang(:,:)
!
   if (deb_Path) write(iout,*), 'mdrest '
!
!
!  check steps between center of mass motion removal
!
   if (.not.dorest)  return
   if (mod(istep,irest) .ne. 0)  return
!
!     perform dynamic allocation of some local arrays
!
   allocate (totmass(0:ngrp))
   allocate (etrans(0:ngrp))
   allocate (vtot(3,0:ngrp))
   if (.not.use_bounds .or. ngrp.ne.0) then
      allocate (erot(0:ngrp))
      allocate (xtot(0:ngrp))
      allocate (ytot(0:ngrp))
      allocate (ztot(0:ngrp))
      allocate (vang(3,0:ngrp))
   end if
!
!  zero out the total mass and overall linear velocity
!
   do i = 0, ngrp
     totmass(i) = 0.0d0
     do j = 1, 3
        vtot(j,i) = 0.0d0
     end do
   end do
!
!     compute linear velocity of the system center of mass
!
      do i = 0, ngrp
         do k = igrp(1,i), igrp(2,i)
            m = kgrp(k)
            if (repart(m).ne.rank) cycle
            weigh = mass(m)
            totmass(i) = totmass(i) + weigh
            do j = 1, 3
               vtot(j,i) = vtot(j,i) + v(j,m)*weigh
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,vtot,3*(ngrp+1),MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,totmass,ngrp+1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
!
!     compute translational kinetic energy of overall system
!
      do i  = 0, ngrp
        etrans(i) = 0.0d0
        do j = 1, 3
           vtot(j,i) = vtot(j,i) / totmass(i)
           etrans(i) = etrans(i) + vtot(j,i)**2
        end do
        etrans(i) = 0.5d0 * etrans(i) * totmass(i) / convert
      end do
!
!     find the center of mass coordinates of each atom group
!
      if (.not.use_bounds .or. ngrp.ne.0) then
         do i = 0, ngrp
            xtot(i) = 0.0d0
            ytot(i) = 0.0d0
            ztot(i) = 0.0d0
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               if (repart(m).ne.rank) cycle
               weigh = mass(m)
               xtot(i) = xtot(i) + x(m)*weigh
               ytot(i) = ytot(i) + y(m)*weigh
               ztot(i) = ztot(i) + z(m)*weigh
            end do
           call MPI_ALLREDUCE(MPI_IN_PLACE,xtot(i),1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,ytot(i),1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,ztot(i),1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
            xtot(i) = xtot(i) / totmass(i)
            ytot(i) = ytot(i) / totmass(i)
            ztot(i) = ztot(i) / totmass(i)
!
!     compute the angular momentum of each atom group
!
            do j = 1, 3
               mang(j) = 0.0d0
            end do
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               if (repart(m).ne.rank) cycle
               weigh = mass(m)
               mang(1) = mang(1) + (y(m)*v(3,m)-z(m)*v(2,m))*weigh
               mang(2) = mang(2) + (z(m)*v(1,m)-x(m)*v(3,m))*weigh
               mang(3) = mang(3) + (x(m)*v(2,m)-y(m)*v(1,m))*weigh
            end do
           call MPI_ALLREDUCE(MPI_IN_PLACE,mang,3,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
            mang(1) = mang(1) - (ytot(i)*vtot(3,i)-ztot(i)*vtot(2,i))*totmass(i)
            mang(2) = mang(2) - (ztot(i)*vtot(1,i)-xtot(i)*vtot(3,i))*totmass(i)
            mang(3) = mang(3) - (xtot(i)*vtot(2,i)-ytot(i)*vtot(1,i))*totmass(i)
!
!     calculate the moment of inertia tensor
!
            xx = 0.0d0
            xy = 0.0d0
            xz = 0.0d0
            yy = 0.0d0
            yz = 0.0d0
            zz = 0.0d0
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               if (repart(m).ne.rank) cycle
               weigh = mass(m)
               xdel = x(m) - xtot(i)
               ydel = y(m) - ytot(i)
               zdel = z(m) - ztot(i)
               xx = xx + xdel*xdel*weigh
               xy = xy + xdel*ydel*weigh
               xz = xz + xdel*zdel*weigh
               yy = yy + ydel*ydel*weigh
               yz = yz + ydel*zdel*weigh
               zz = zz + zdel*zdel*weigh
            end do
           call MPI_ALLREDUCE(MPI_IN_PLACE,xx,1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,xy,1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,xz,1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,yy,1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,yz,1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,zz,1,MPI_REAL8,MPI_SUM,COMM_TINKER,ierr)
            tensor(1,1) = yy + zz
            tensor(2,1) = -xy
            tensor(3,1) = -xz
            tensor(1,2) = -xy
            tensor(2,2) = xx + zz
            tensor(3,2) = -yz
            tensor(1,3) = -xz
            tensor(2,3) = -yz
            tensor(3,3) = xx + yy
!
!     fix to avoid singularity for one- or two-body groups
!
            if (igrp(2,i)-igrp(1,i) .le. 2) then
               eps = 0.000001d0
               tensor(1,1) = tensor(1,1) + eps
               tensor(2,2) = tensor(2,2) + eps
               tensor(3,3) = tensor(3,3) + eps
            end if
!
!     diagonalize the moment of inertia tensor
!
            call invert (3,tensor)
!
!     compute angular velocity and rotational kinetic energy
!
            erot(i) = 0.0d0
            do k = 1, 3
               vang(k,i) = 0.0d0
               do j = 1, 3
                  vang(k,i) = vang(k,i) + tensor(k,j)*mang(j)
               end do
               erot(i) = erot(i) + vang(k,i)*mang(k)
            end do
            erot(i) = 0.5d0 * erot(i) / convert
         end do
      end if
!
!     eliminate any translation of each atom group
!
      do i = 0, ngrp
         do k = igrp(1,i), igrp(2,i)
            m = kgrp(k)
            if (repart(m).ne.rank) cycle
            do j = 1, 3
               v(j,m) = v(j,m) - vtot(j,i)
            end do
         end do
      end do
!
!     print the translational velocity of each atom group
!
      if (debug.and.rank.eq.0) then
         write (iout,10)
   10    format ()
         if (ngrp .eq. 0) then
            write (iout,20)  (vtot(i,0),i=1,3),etrans(0)
   20       format (' System Linear Velocity :  ',3d12.2,&
                   /,' Translational Kinetic Energy :',10x,f12.4,&
                      ' Kcal/mole')
         else
            do i = 0, ngrp
               write (iout,30)  i,(vtot(j,i),j=1,3),etrans(i)
   30          format (' Group',i4,' Linear Velocity :  ',3d12.2,&
                      /,' Translational Kinetic Energy :',10x,f12.4,&
                         ' Kcal/mole')
            end do
         end if
      end if
!
!     eliminate any rotation about each group center of mass
!
      if (.not.use_bounds .or. ngrp.ne.0) then
         do i = 0, ngrp
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               if (repart(m).ne.rank) cycle
               xdel = x(m) - xtot(i)
               ydel = y(m) - ytot(i)
               zdel = z(m) - ztot(i)
               v(1,m) = v(1,m) - vang(2,i)*zdel + vang(3,i)*ydel
               v(2,m) = v(2,m) - vang(3,i)*xdel + vang(1,i)*zdel
               v(3,m) = v(3,m) - vang(1,i)*ydel + vang(2,i)*xdel
            end do
         end do
!
!     print the angular velocity of each atom group
!
         if (debug.and.rank.eq.0) then
            if (ngrp .eq. 0) then
               write (iout,40)  (vang(j,0),j=1,3),erot(0)
   40          format (' System Angular Velocity : ',3d12.2,&
                      /,' Rotational Kinetic Energy :',13x,f12.4,&
                         ' Kcal/mole')
            else
               do i = 0, ngrp
                  write (iout,50)  i,(vang(j,i),j=1,3),erot(i)
   50             format (' Group',i4,' Angular Velocity : ',3d12.2,&
                         /,' Rotational Kinetic Energy :',13x,f12.4,&
                            ' Kcal/mole')
               end do
            end if
         end if
      end if
!
!     perform deallocation of some local arrays
!
      deallocate (totmass)
      deallocate (etrans)
      deallocate (vtot)
      if (.not.use_bounds .or. ngrp.ne.0) then
         deallocate (erot)
         deallocate (xtot)
         deallocate (ytot)
         deallocate (ztot)
         deallocate (vang)
      end if
      return
end
