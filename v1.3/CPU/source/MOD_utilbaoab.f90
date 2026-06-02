module utilbaoab
   implicit none

contains

!> @brief 
!> apply the A-block of the BAOAB integrator: update of positions given velocities
!> @param[in] dt: time step
   subroutine apply_a_block(dt)
      use atoms
      use atmlst
      use inform
      use iounit
      use moldyn
      use usage
      use domdec
      use freeze
      implicit none
      real*8, intent(in) :: dt
      integer :: i,iglob
!
      if (deb_Path) write(iout,*), 'apply_a_block '
!

      if (use_rattle) then
         call save_atoms_pos
      endif

      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            x(iglob) = x(iglob) + v(1,iglob)*dt
            y(iglob) = y(iglob) + v(2,iglob)*dt
            z(iglob) = z(iglob) + v(3,iglob)*dt
         end if
      end do
!
      if (use_rattle) then
         call rattle(dt)
         call rattle2(dt)
      endif
   end subroutine apply_a_block

!> @brief 
!> apply the O-block of the BAOAB integrator: analytical integration of Orstein-Uhlenbeck process to update velocities
!> @param[in] dt: time step
   subroutine apply_o_block(dt)
      use atoms
      use atmtyp
      use bath
      use inform
      use iounit
      use usage
      use domdec
      use freeze
      use moldyn
      use qtb, only:qtbrandom,qtb_thermostat
      use units
      use langevin, only: gamma
      implicit none
      real*8, intent(in) :: dt
      real*8 :: a1,a2,Rn
      integer :: i,j,iglob
      real*8 :: normal
!
      if (deb_Path) write(iout,*), 'apply_o_block '
!

      if (use_rattle) then
         call rattle2(dt)
      endif

      if(qtb_thermostat) then
         call qtbrandom()
      else
         a1 = exp(-gamma*dt)
         a2 = sqrt((1-a1**2)*boltzmann*kelvin)

         do i=1,nloc
            iglob = glob(i)
            if (use(iglob)) then
               do j=1,3
                  Rn=normal()
                  v(j,iglob) = a1*v(j,iglob) +&
                  &a2*Rn/sqrt(mass(iglob))
               enddo
            endif
         enddo
      endif
!
      if (use_rattle) then
         call rattle2(dt)
      endif
   end subroutine apply_o_block

!> @brief 
!> apply the B-block of the BAOAB integrator for Langevin Piston: update of velocities given acceleration (forces)
!> @param[in] dt: time step
!> @param[in] pres: current pressure
!> @param[in] stress: current stress tensor
   subroutine apply_b_piston(dt,pres,stress)
      use bath
      use inform
      use iounit
      use units
      implicit none
      real*8, intent(in) :: dt,pres,stress(3,3)
      integer :: i
!
      if (deb_Path) write(iout,*), 'apply_b_piston '
!

      if (anisotrop) then
         do i=1,3
            aextbox(i) = convert*(extvol*(stress(i,i)-atmsph)/prescon&
            &+gasconst*kelvin)/masspiston
            vextbox(i) = vextbox(i) + dt*aextbox(i)
         enddo

      else
         aextvol = 3.0*convert*(extvol*(pres-atmsph)/prescon&
         &+gasconst*kelvin)/masspiston
         vextvol = vextvol + dt*aextvol
      endif

   end subroutine apply_b_piston

!> @brief 
!> apply the O-block of the BAOAB integrator for the Langevin-Piston: analytical integration of Orstein-Uhlenbeck process to update velocities
!> @param[in] dt: time step
   subroutine apply_o_piston(dt)
      use bath
      use inform
      use iounit
      use units
      use mpi
      use domdec
      implicit none
      real*8, intent(in) :: dt
      real*8 :: a1piston,a2piston
      integer :: ierr,i,naxis
      interface
         function normal()
            real*8 normal
         end function normal
      end interface
!
      if (deb_Path) write(iout,*), 'apply_o_piston '
!

      if (rank.eq.0) then
         a1piston = exp(-gammapiston*dt)
         a2piston = sqrt((1.d0-a1piston**2)*&
         &boltzmann*kelvin/masspiston)
         if(anisotrop) then
            do i=1,3
               vextbox(i)=a1piston*vextbox(i)+a2piston*normal()
            enddo
         else
            vextvol=a1piston*vextvol + a2piston*normal()
         endif
      end if
      if(nproc>1) then
         if (anisotrop) then
            call MPI_BCAST(vextbox,3,MPI_REAL8,0,COMM_TINKER,ierr)
         else
            call MPI_BCAST(vextvol,1,MPI_REAL8,0,COMM_TINKER,ierr)
         endif
      endif
      if (anisotrop) then
         temppiston = 0
         naxis = 0
         do i=1,3
            if (.not. freeze_axis(i)) then
               naxis = naxis + 1
               temppiston = temppiston + masspiston*vextbox(i)**2
            end if
         end do
         temppiston = temppiston/(naxis*boltzmann)
      else
         temppiston = masspiston*vextvol**2/boltzmann
      endif

   end subroutine apply_o_piston

!> @brief 
!> apply the A-block of the BAOAB integrator for the Langevin-Piston: update of positions given velocities
!> @param[in] dt: time step
   subroutine apply_a_piston(dt,istep,A_full)
      use atoms
      use atmlst
      use inform
      use iounit
      use moldyn
      use usage
      use domdec
      use freeze
      use bath
      implicit none
      real*8, intent(in) :: dt
      integer, intent(in) :: istep
      logical, intent(in) :: A_full
      integer :: i,iglob
      real*8 :: scale(3)
      real*8 :: tau(3)
!
      if (deb_Path) write(iout,*), 'apply_a_piston '
!

      if (use_rattle) then
         call save_atoms_pos
      endif

      extvolold = extvol
      tau=0.d0
      if(anisotrop) then
         scale(:) = exp(dt*vextbox(:))
         extbox(:)=extbox(:)*scale(:)
         if (A_full) tau=sinh(dt*vextbox)/vextbox
      else
         scale(:) = exp(dt*vextvol)
         if (A_full) tau=sinh(dt*vextvol)/vextvol
      endif
      extvol = extvol*scale(1)*scale(2)*scale(3)

      do i =1,3
         if (freeze_axis(i)) then
            scale(i)=1.d0
            tau(i)=dt
         end if
      enddo
      call rescale_box(istep,scale)

      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            x(iglob) = x(iglob)*scale(1) + v(1,iglob)*tau(1)
            y(iglob) = y(iglob)*scale(2) + v(2,iglob)*tau(2)
            z(iglob) = z(iglob)*scale(3) + v(3,iglob)*tau(3)
            v(:,iglob)= v(:,iglob)/scale(:)
            !a(:,iglob)= a(:,iglob)/scale
         end if
      end do
!
      if (use_rattle) then
         call rattle(dt)
         call rattle2(dt)
      endif
   end subroutine apply_a_piston



end module utilbaoab
