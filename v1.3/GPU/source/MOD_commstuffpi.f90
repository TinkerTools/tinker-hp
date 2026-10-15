!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"

module commstuffpi
   use parallel, only: commposrec, migrate_polymer_atoms, &
      &exchange_polymer_positions, exchange_polymer_recip_positions, &
      &exchange_polymer_gradient_state, exchange_polymer_normal_modes
   implicit none
   private :: exchange_polymer_gradient_state, exchange_polymer_normal_modes

contains

   subroutine full_gradient_prepare_pi(polymer,polymer_ctr&
      &,istep,respa)
      use beads
      use domdec
      use atomsMirror
      use cutoff, only: use_list
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
      integer, intent(in) :: istep
      logical, intent(in) :: respa
      integer :: iloc,i

      !! UPDATE DIRECT SPACE BEFORE GRADIENT CALCULATION !!
      call update_direct_space_pi(polymer)

      !! COMMUNICATE FOR GRADIENT CALCULATIONS !!
      !! PROPAGATION REPARTITION -> GRADIENT REPARTITION !!
      call exchange_polymer_gradient_state(polymer,.TRUE.,send_vel=nproc>1)
      if(contract .AND. polymer_ctr%nbeads>1) then
         !! PREPARE AND COMMUNICATE CONTRACTED POLYMER !!
         call contract_polymer(polymer,polymer_ctr)
         call exchange_polymer_gradient_state(polymer_ctr,.FALSE.)
         !! REASSIGN DOFs ACCORDING TO CENTROID POSITION !!
         !! -> COMMUNICATE ALSO CONTRACTED POLYMER !!
         call migrate_polymer_atoms(polymer,polymer_ctr)
         !! GET NECESSARY CONTRACTED DOFs FROM NEIGHBOURING PROCS !!
         call exchange_polymer_positions(polymer_ctr,.FALSE.,.FALSE.)
      else
         !! REASSIGN DOFs ACCORDING TO CENTROID POSITION !!
         call migrate_polymer_atoms(polymer)
      endif
      !! GET NECESSARY DOFs FROM NEIGHBOURING PROCS !!
      call exchange_polymer_positions(polymer,.TRUE.,.FALSE.)

      !! LOAD CENTROID !!
      call load_bead(polymer,0)

      !! REASSIGN AND COMMUNICATE POSITIONS FOR RECIPROCAL SPACE !!
      call reassignpme(.false.)
      call commposrec
      if(nproc>1) then
         !SYNCHRONIZE CENTROID
!$acc parallel loop async default(present)
         DO iloc=1,nbloc
            i=glob(iloc)
            polymer%eigpos(1,i,1)=x(i)
            polymer%eigpos(2,i,1)=y(i)
            polymer%eigpos(3,i,1)=z(i)
         ENDDO
         ! REORGANIZE POSITIONS AND VELOCITIES AFTER REASSIGN
         if(respa) then
            call exchange_polymer_normal_modes(polymer&
               &,polymer%pos,polymer%vel,polymer%forces)
         else
            call exchange_polymer_normal_modes(polymer&
               &,polymer%pos,polymer%vel)
         endif
         call update_normal_modes_pi(polymer)
      endif

      if(.not. centroid_recip) then
         if(contract .AND. nbeads_ctr>1) then
            call exchange_polymer_recip_positions(polymer_ctr,.FALSE.)
         else
            call exchange_polymer_recip_positions(polymer,.FALSE.)
         endif
      endif

      !! CENTROID IS LOADED !!
      call reCast_position
      call reinitnl(istep)
      if(respa) then
         call mechanicsteprespa(istep,.FALSE.)
         call allocsteprespa(.false.)
      else
         call mechanicstep(istep)
         call allocstep
      endif
      if (use_list) call nblist(istep)

   end subroutine full_gradient_prepare_pi

   subroutine fast_gradient_prepare_pi(polymer,istep)
      use beads
      use domdec
      use atomsMirror
      use cutoff, only: use_list
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer, intent(in) :: istep
      integer :: iloc,i

      !! UPDATE DIRECT SPACE BEFORE GRADIENT CALCULATION !!
      call update_direct_space_pi(polymer)

      !! COMMUNICATE FOR GRADIENT CALCULATIONS !!
      !! PROPAGATION REPARTITION -> GRADIENT REPARTITION !!
      !! (don't send velocities as there is no reassign) !!
      call exchange_polymer_gradient_state(polymer,.TRUE.,.false.)
      !! GET NECESSARY DOFs FROM NEIGHBOURING PROCS (RESPA FAST) !!
      call exchange_polymer_positions(polymer,.TRUE.,.TRUE.)

      !! LOAD CENTROID !!
      call load_bead(polymer,0)
      call mechanicsteprespa(istep,.true.)
      call allocsteprespa(.true.)

   end subroutine fast_gradient_prepare_pi









   subroutine commforcesblocpi(derivs,nbeadsproc,fast)
      use beads
      use deriv
      use domdec
      use inform,only: deb_Path
      use potent
      use sizes
      use timestat
      use mdevice  ,only: reciprocal_work
      use parallel  ,only: accumulate,halo_large,halo_neighbour
      implicit none
      integer, intent(in) :: nbeadsproc
      real(r_p), target, intent(inout) :: derivs(3,nbloc,*)
      logical, intent(in) :: fast
      real(r_p), pointer :: flat(:)
      integer k, halo

!
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      if (deb_Path) write(*,*) '   >> commforcesblocpi'

      halo = merge(halo_neighbour,halo_large,fast)
      do k = 1, nbeadsproc
         flat(1:3*nbloc) => derivs(:,:,k)
         call accumulate(flat,3,halo,channel=reciprocal_work)
      end do

      if (deb_Path) write(*,*) '   << commforcesblocpi'
   end subroutine commforcesblocpi






   subroutine save_deriv_rec(temprec)
      use potent ,pa=>PotentialAll
      use domdec
      use deriv
      implicit none
      real(r_p), intent(inout) :: temprec(3,nlocrec2)
      integer :: i,j

      if (pa) then
!$acc parallel loop collapse(2) async &
!$acc   present(demrec,deprec,decrec,dedsprec)
         do i = 1, nlocrec2
            do j = 1, 3
               temprec(j,i)= demrec(j,i)&
                  &+ deprec(j,i)&
                  &+ decrec(j,i)&
                  &+ dedsprec(j,i)
               demrec(j,i)=0.d0
               deprec(j,i)=0.d0
               decrec(j,i)=0.d0
               dedsprec(j,i)=0.d0
            end do
         end do
      else
         if (use_polar) then
!$acc parallel loop collapse(2) async &
!$acc      present(demrec,deprec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)= demrec(j,i) + deprec(j,i)
                  demrec(j,i)=0.d0
                  deprec(j,i)=0.d0
               end do
            end do
         else if (use_mpole) then
!$acc parallel loop collapse(2) async &
!$acc      present(demrec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)= demrec(j,i)
                  demrec(j,i)=0.d0
               end do
            end do
         end if
         if (use_charge) then
!$acc parallel loop collapse(2) async &
!$acc      present(decrec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)=  decrec(j,i)
                  decrec(j,i)=0.d0
               end do
            end do
         end if
         if (use_disprec) then
!$acc parallel loop collapse(2) async &
!$acc      present(decrec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)=temprec(j,i) + dedsprec(j,i)
                  dedsprec(j,i)=0.d0
               end do
            end do
         end if
      end if


   end subroutine save_deriv_rec






end module commstuffpi
