!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (parallel) parallel_communication
   use mdevice, only: channel_t
   use atomsMirror, only: x, y, z
   use atoms, only: n
   use domdec
   use group
   use mpole, only: polerecloc, poleloc
   use potent, only: use_pmecore
   use mpi
   use status, only: clear_status, set_status, status_invalid_argument, status_t
   use lifecycle_observer, only: LifecycleTicket, record_lifecycle_event
   implicit none

   integer, save :: active_polar_ticket_count = 0
   logical, save :: inject_polar_finish_failure = .false.

contains

module subroutine parallel_ctf(xr, yr, zr)
!$acc routine seq
   use boxes
   use cell
   implicit none
   real(t_p), intent(inout) :: xr, yr, zr

   zr = zr / gamma_term
   yr = (yr - zr*beta_term) / gamma_sin
   xr = xr - yr*gamma_cos - zr*beta_cos
end subroutine parallel_ctf

module subroutine parallel_image(xr, yr, zr)
!$acc routine seq
   use cell
   implicit none
   real(t_p), intent(inout) :: xr, yr, zr

   if (abs(xr).gt.xcell2) xr = xr-sign(xcell,xr)
   if (abs(yr).gt.ycell2) yr = yr-sign(ycell,yr)
   if (abs(zr).gt.zcell2) zr = zr-sign(zcell,zr)
end subroutine parallel_image

module subroutine activate_domain(platform)
   type(parallel_t), target, intent(inout) :: platform

   active_domain => platform%domain
   active_platform => platform
end subroutine activate_domain

module subroutine deactivate_domain(platform)
   type(parallel_t), target, intent(inout) :: platform

   if (associated(active_domain, platform%domain)) nullify(active_domain)
   if (associated(active_platform, platform)) nullify(active_platform)
end subroutine deactivate_domain

module subroutine create_domain_subcommunicators(direct_count, process_rank)
   integer, intent(in) :: direct_count
   integer, intent(out) :: process_rank
   integer :: total_group, direct_group, reciprocal_group, ierr, i
   integer, allocatable :: direct_ranks(:)

   if (.not.associated(active_platform)) then
      call abort_on_mpi_error(MPI_ERR_COMM, 'parallel platform is not active')
      return
   end if
   if (direct_count < 0 .or. direct_count > active_platform%domain%size) then
      call abort_on_mpi_error(MPI_ERR_ARG, 'invalid direct-space communicator size')
      return
   end if
   if (active_platform%domain%direct_comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(active_platform%domain%direct_comm, ierr)
      call abort_on_mpi_error(ierr, 'MPI_Comm_free(direct communicator)')
   end if
   if (active_platform%domain%recip_comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(active_platform%domain%recip_comm, ierr)
      call abort_on_mpi_error(ierr, 'MPI_Comm_free(reciprocal communicator)')
   end if

   allocate(direct_ranks(active_platform%domain%size))
   call MPI_Comm_group(active_platform%domain%communicator, total_group, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Comm_group(domain communicator)')
   do i = 1, direct_count
      direct_ranks(i) = i - 1
   end do
   call MPI_Group_incl(total_group, direct_count, direct_ranks, direct_group, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Group_incl(direct ranks)')
   call MPI_Group_excl(total_group, direct_count, direct_ranks, reciprocal_group, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Group_excl(reciprocal ranks)')
   call MPI_Comm_create(active_platform%domain%communicator, direct_group, &
      &active_platform%domain%direct_comm, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Comm_create(direct communicator)')
   call MPI_Comm_create(active_platform%domain%communicator, reciprocal_group, &
      &active_platform%domain%recip_comm, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Comm_create(reciprocal communicator)')
   call MPI_Group_free(direct_group, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Group_free(direct group)')
   call MPI_Group_free(reciprocal_group, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Group_free(reciprocal group)')
   call MPI_Group_free(total_group, ierr)
   call abort_on_mpi_error(ierr, 'MPI_Group_free(total group)')
   deallocate(direct_ranks)

   if (active_platform%domain%rank < direct_count) then
      call MPI_Comm_rank(active_platform%domain%direct_comm, process_rank, ierr)
      call abort_on_mpi_error(ierr, 'MPI_Comm_rank(direct communicator)')
   else
      call MPI_Comm_rank(active_platform%domain%recip_comm, process_rank, ierr)
      call abort_on_mpi_error(ierr, 'MPI_Comm_rank(reciprocal communicator)')
   end if
end subroutine create_domain_subcommunicators

subroutine abort_on_mpi_error(error_code, operation)
   integer, intent(in) :: error_code
   character(*), intent(in) :: operation
   integer :: abort_error

   if (error_code == MPI_SUCCESS) return
   write(*,'(A,I0,2A)') 'MPI failure (', error_code, ') in ', trim(operation)
   call MPI_Abort(MPI_COMM_WORLD, error_code, abort_error)
   error stop 'MPI operation failed'
end subroutine abort_on_mpi_error

module subroutine polar_set_finish_failure_legacy(enabled)
   logical, intent(in) :: enabled

   inject_polar_finish_failure = enabled
end subroutine polar_set_finish_failure_legacy

subroutine report_polar_finish(status_value)
   type(status_t), intent(out), optional :: status_value

   if (present(status_value)) call clear_status(status_value)
   if (.not.inject_polar_finish_failure) return
   if (present(status_value)) call set_status(status_value, status_invalid_argument, &
      &'Injected active polarization-ticket finalization failure')
   inject_polar_finish_failure = .false.
end subroutine report_polar_finish

module subroutine gather_positions()
   implicit none
   integer, parameter :: Components = 3
   integer :: i, offset, offset_values, global_index, ierr
   integer :: send_offsets(nproc), receive_sizes(nproc), receive_offsets(nproc)
   integer :: global_indices(n)
   real(r_p) :: positions(Components*n)

   if (nproc.eq.1.or.nbig_recep.eq.nproc) return

!$acc wait
!$acc data create(global_indices,positions)
   send_offsets(1) = 0
   do i = 2, nproc
      send_offsets(i) = send_offsets(i-1) + domlen(i-1)
   end do
   offset = send_offsets(rank+1) + 1
   offset_values = Components*send_offsets(rank+1) + 1

!$acc parallel loop default(present)
   do i = 0, nloc-1
      global_index = glob(i+1)
      global_indices(offset+i) = global_index
      positions(offset_values+i*Components) = x(global_index)
      positions(offset_values+i*Components+1) = y(global_index)
      positions(offset_values+i*Components+2) = z(global_index)
   end do

   receive_sizes(:) = Components*domlen(1:nproc)
   receive_offsets(:) = Components*send_offsets(1:nproc)

!$acc host_data use_device(global_indices,positions)
   call MPI_Allgatherv(MPI_IN_PLACE,1,MPI_DATATYPE_NULL,global_indices,domlen, &
      &send_offsets,MPI_INT,get_domain_comm(),ierr)
   call MPI_Allgatherv(MPI_IN_PLACE,1,MPI_DATATYPE_NULL,positions,receive_sizes, &
      &receive_offsets,MPI_RPREC,get_domain_comm(),ierr)
!$acc end host_data

!$acc parallel loop default(present)
   do i = 0, n-1
      global_index = global_indices(i+1)
      x(global_index) = positions(Components*i+1)
      y(global_index) = positions(Components*i+2)
      z(global_index) = positions(Components*i+3)
   end do
!$acc end data
end subroutine gather_positions

module subroutine assign_direct()
   use atoms, only: x_primary => x, y_primary => y, z_primary => z
   use boxes
   use cell
   use inform, only: deb_Path
   use potent
   use timestat, only: timer_eneig, quiet_timers, timer_enter, timer_exit
   implicit none
   integer :: i, iproc, total_assigned
   integer :: local_rank, n_ar, n_br, n_cr
   real(t_p) :: xr, yr, zr, lena_box, lenb_box, lenc_box

   if (deb_Path) print*, '  assign_direct'
   call build_domain_delimiters

   if (nproc.eq.1.or.use_pmecore.and.ndir.eq.1) return
   call timer_enter(timer_eneig)

   if (use_pmecore) then
      local_rank = rank_bis
   else
      local_rank = rank
   end if

   lena_box = xcell/nxdd
   lenb_box = ycell/nydd
   lenc_box = zcell/nzdd

!$acc parallel loop async default(present)
   do i = 1, n
      glob(i) = 0
      repart(i) = -1
      repartrec(i) = -1
      loc(i) = 0
   end do

   domlen(:) = 0

!$acc parallel loop async default(present) copy(domlen)
   do i = 1, n
      xr = x_primary(i)
      yr = y_primary(i)
      zr = z_primary(i)
      if (monoclinic.or.triclinic) call parallel_ctf(xr,yr,zr)
      call parallel_image(xr,yr,zr)
      if ((xcell2-abs(xr)).lt.eps_cell) xr = xr-sign(4*eps_cell,xr)
      if ((ycell2-abs(yr)).lt.eps_cell) yr = yr-sign(4*eps_cell,yr)
      if ((zcell2-abs(zr)).lt.eps_cell) zr = zr-sign(4*eps_cell,zr)
      n_ar = int((xr+xcell2)/lena_box)
      n_br = int((yr+ycell2)/lenb_box)
      n_cr = int((zr+zcell2)/lenc_box)
      iproc = n_ar + n_br*nxdd + n_cr*nydd*nxdd
      repart(i) = iproc
!$acc atomic
      domlen(iproc+1) = domlen(iproc+1) + 1
   end do

!$acc wait
   if (use_pmecore.and.rank.le.ndir-1) then
      nloc = domlen(local_rank+1)
   else if (use_pmecore.and.rank.gt.ndir-1) then
      nloc = 0
   else
      nloc = domlen(local_rank+1)
   end if

   total_assigned = 0
   do i = 1, nproc
      total_assigned = total_assigned + domlen(i)
   end do
   if (total_assigned.ne.n) then
      write(0,'(A,I0,A,I0)') ' ERROR assign_direct: rank ', rank, &
         &' total atoms assigned ', total_assigned
      __TINKER_FATAL__
   end if

   call orderbuffer_gpu(.true.)
   if (associated(active_domain)) call refresh_domain_layout(active_domain)
   call timer_exit(timer_eneig,quiet_timers)
end subroutine assign_direct

module subroutine initialize_direct_communication()
   call orderbuffer(.true.)
end subroutine initialize_direct_communication

module subroutine assign_reciprocal()
   use atoms, only: x_primary => x, y_primary => y, z_primary => z
   use boxes, only: recip
   use chunks, only: grdoff
   use fft, only: kstart1, kend1, jstart1, jend1
   use inform, only: deb_Path
   use pme, only: pme_eps, nfft1, nfft2, nfft3, bsorder
   use potent, only: use_pmecore
   use tinheader, only: ti_p
   use timestat, only: timer_eneig, quiet_timers, timer_enter, timer_exit
   implicit none
   integer :: i, j, k, ifr, iproc, process_count
   real(t_p) :: xi, yi, zi, w, fractional

   if (use_pmecore) then
      process_count = ndir
   else
      process_count = nproc
   end if
   if (deb_Path) print*, 'assign_reciprocal'
   call timer_enter(timer_eneig)

!$acc parallel loop async default(present)
   do i = 1, n
      xi = x_primary(i)
      yi = y_primary(i)
      zi = z_primary(i)
      w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
      fractional = real(nfft3,t_p)*(w-anint(w)+0.5_ti_p)
      ifr = int(fractional-pme_eps)
      k = ifr-bsorder+grdoff
      if (k.lt.1) k = k+nfft3
      w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
      fractional = real(nfft2,t_p)*(w-anint(w)+0.5_ti_p)
      ifr = int(fractional-pme_eps)
      j = ifr-bsorder+grdoff
      if (j.lt.1) j = j+nfft2
!$acc loop seq
      do iproc = 0, process_count-1
         if ((k.ge.kstart1(iproc+1).and.k.le.kend1(iproc+1)).and. &
            &(j.ge.jstart1(iproc+1).and.j.le.jend1(iproc+1))) then
            repartrec(i) = iproc
         end if
      end do
   end do
   if (nproc.eq.1) then
!$acc wait
!$acc update self(repartrec(1:n))
      call orderbufferrec
   else
      call orderbufrec_gpu
   end if
   call timer_exit(timer_eneig,quiet_timers)
end subroutine assign_reciprocal

module subroutine polar_post_comm_legacy(ticket, nrhs, values, use_short_range, send_values, group_comm, channel)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   logical, intent(in) :: use_short_range
   real(t_p), intent(inout) :: send_values(3,nrhs,*)
   logical, intent(in) :: group_comm
   type(channel_t), intent(in) :: channel
   integer, allocatable :: receive_requests(:), send_requests(:)

   if (ticket%active .or. ticket%dip_active) error stop 'polar exchange ticket is already active'
   allocate(receive_requests(max(1,nproc)), send_requests(max(1,nproc)))
   receive_requests = MPI_REQUEST_NULL
   send_requests = MPI_REQUEST_NULL
   ticket%active = .true.
   ticket%use_short_range = use_short_range
   ticket%use_group = group_comm
   ticket%use_gpu = .not.use_short_range
   if (nproc.eq.1) then
      call move_alloc(receive_requests,ticket%recv_req)
      call move_alloc(send_requests,ticket%send_req)
      active_polar_ticket_count = active_polar_ticket_count + 1
      return
   end if
   if (ticket%use_group) then
      call post_group_comm(nrhs, values, receive_requests, send_requests)
   else if (ticket%use_short_range) then
      call commdirdirshort(nrhs, 1, send_values, receive_requests, send_requests, channel)
      call commdirdirshort(nrhs, 0, values, receive_requests, send_requests, channel)
   else
      call commdirdirgpu(nrhs, 0, values, receive_requests, send_requests)
      call commdirdirgpu(nrhs, 1, send_values, receive_requests, send_requests)
   end if
   call move_alloc(receive_requests,ticket%recv_req)
   call move_alloc(send_requests,ticket%send_req)
   active_polar_ticket_count = active_polar_ticket_count + 1
end subroutine polar_post_comm_legacy

module subroutine polar_finish_legacy(ticket, nrhs, values, status_value, channel)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   type(status_t), intent(out), optional :: status_value
   type(channel_t), intent(in) :: channel

   if (.not.ticket%active) error stop 'polar exchange ticket is not active'
   if (nproc.eq.1) then
      call reset_polar_ticket(ticket)
      call report_polar_finish(status_value)
      return
   end if
   if (ticket%use_group) then
      call finish_group_comm(ticket%recv_req, ticket%send_req)
   else if (ticket%use_short_range) then
      call commdirdirshort(nrhs, 2, values, ticket%recv_req, ticket%send_req, channel)
   else
      call commdirdirgpu(nrhs, 2, values, ticket%recv_req, ticket%send_req)
   end if
   call reset_polar_ticket(ticket)
   call report_polar_finish(status_value)
end subroutine polar_finish_legacy

module subroutine polar_post_dip_comm_legacy(ticket, nrhs, diprec, dip, &
   &buffer_send, buffer_recv, use_gpu)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
   real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
   logical, intent(in) :: use_gpu
   integer, allocatable :: receive_requests(:), send_requests(:)

   if (ticket%active .or. ticket%dip_active) error stop 'polar exchange ticket is already active'
   allocate(receive_requests(max(1,nproc)), send_requests(max(1,nproc)))
   receive_requests = MPI_REQUEST_NULL
   send_requests = MPI_REQUEST_NULL
   ticket%dip_active = .true.
   if (nproc.eq.1) then
      call move_alloc(receive_requests,ticket%dip_recv_req)
      call move_alloc(send_requests,ticket%dip_send_req)
      active_polar_ticket_count = active_polar_ticket_count + 1
      return
   end if
   if (use_gpu) then
      call commrecdirdipgpu(nrhs, 0, diprec, dip, buffer_send, buffer_recv, &
         &receive_requests, send_requests)
      call commrecdirdipgpu(nrhs, 1, diprec, dip, buffer_send, buffer_recv, &
         &receive_requests, send_requests)
   else
      call commrecdirdip(nrhs, 0, diprec, dip, buffer_send, buffer_recv, &
         &receive_requests, send_requests)
      call commrecdirdip(nrhs, 1, diprec, dip, buffer_send, buffer_recv, &
         &receive_requests, send_requests)
   end if
   call move_alloc(receive_requests,ticket%dip_recv_req)
   call move_alloc(send_requests,ticket%dip_send_req)
   active_polar_ticket_count = active_polar_ticket_count + 1
end subroutine polar_post_dip_comm_legacy

module subroutine polar_finish_dip_comm_legacy(ticket, nrhs, diprec, dip, &
   &buffer_send, buffer_recv, use_gpu, status_value)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
   real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
   logical, intent(in) :: use_gpu
   type(status_t), intent(out), optional :: status_value

   if (.not.ticket%dip_active) error stop 'polar exchange ticket is not active'
   if (nproc.eq.1) then
      call reset_polar_ticket(ticket)
      call report_polar_finish(status_value)
      return
   end if
   if (use_gpu) then
      call commrecdirdipgpu(nrhs, 2, diprec, dip, buffer_send, buffer_recv, &
         &ticket%dip_recv_req, ticket%dip_send_req)
   else
      call commrecdirdip(nrhs, 2, diprec, dip, buffer_send, buffer_recv, &
         &ticket%dip_recv_req, ticket%dip_send_req)
   end if
   call reset_polar_ticket(ticket)
   call report_polar_finish(status_value)
end subroutine polar_finish_dip_comm_legacy

module subroutine polar_post_receives_direct_legacy(ticket, nrhs, values, channel)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   type(channel_t), intent(in) :: channel

   if (ticket%active .or. ticket%dip_active .or. ticket%recdir_active) &
      &error stop 'polar exchange ticket is already active'
   allocate(ticket%recv_req(max(1,nproc)), ticket%send_req(max(1,nproc)))
   ticket%recv_req = MPI_REQUEST_NULL
   ticket%send_req = MPI_REQUEST_NULL
   ticket%active = .true.
   active_polar_ticket_count = active_polar_ticket_count + 1
   if (nproc > 1) call commdirdir(nrhs, 0, values, ticket%recv_req, ticket%send_req, channel)
end subroutine polar_post_receives_direct_legacy

module subroutine polar_post_sends_direct_legacy(ticket, nrhs, values, channel)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   type(channel_t), intent(in) :: channel

   if (.not.ticket%active) error stop 'polar exchange ticket is not active'
   if (nproc > 1) call commdirdir(nrhs, 1, values, ticket%recv_req, ticket%send_req, channel)
end subroutine polar_post_sends_direct_legacy

module subroutine polar_finish_direct_legacy(ticket, nrhs, values, status_value, channel)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   type(status_t), intent(out), optional :: status_value
   type(channel_t), intent(in) :: channel

   if (.not.ticket%active) error stop 'polar exchange ticket is not active'
   if (nproc > 1) call commdirdir(nrhs, 2, values, ticket%recv_req, ticket%send_req, channel)
   call reset_polar_ticket(ticket)
   call report_polar_finish(status_value)
end subroutine polar_finish_direct_legacy

module subroutine polar_post_receives_dip_legacy(ticket, nrhs, diprec, dip, &
   &buffer_send, buffer_recv)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
   real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)

   if (ticket%active .or. ticket%dip_active .or. ticket%recdir_active) &
      &error stop 'polar exchange ticket is already active'
   allocate(ticket%dip_recv_req(max(1,nproc)), ticket%dip_send_req(max(1,nproc)))
   ticket%dip_recv_req = MPI_REQUEST_NULL
   ticket%dip_send_req = MPI_REQUEST_NULL
   ticket%dip_active = .true.
   active_polar_ticket_count = active_polar_ticket_count + 1
   if (nproc > 1) call commrecdirdip(nrhs, 0, diprec, dip, buffer_send, buffer_recv, &
      &ticket%dip_recv_req, ticket%dip_send_req)
end subroutine polar_post_receives_dip_legacy

module subroutine polar_post_sends_dip_legacy(ticket, nrhs, diprec, dip, &
   &buffer_send, buffer_recv)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
   real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)

   if (ticket%active .or. ticket%recdir_active) &
      &error stop 'polar exchange ticket is already active'
   if (.not.ticket%dip_active) then
      allocate(ticket%dip_recv_req(max(1,nproc)), ticket%dip_send_req(max(1,nproc)))
      ticket%dip_recv_req = MPI_REQUEST_NULL
      ticket%dip_send_req = MPI_REQUEST_NULL
      ticket%dip_active = .true.
      active_polar_ticket_count = active_polar_ticket_count + 1
   end if
   if (nproc > 1) call commrecdirdip(nrhs, 1, diprec, dip, buffer_send, buffer_recv, &
      &ticket%dip_recv_req, ticket%dip_send_req)
end subroutine polar_post_sends_dip_legacy

module subroutine polar_finish_dip_legacy(ticket, nrhs, diprec, dip, &
   &buffer_send, buffer_recv, status_value)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
   real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
   type(status_t), intent(out), optional :: status_value

   if (.not.ticket%dip_active) error stop 'polar dipole ticket is not active'
   if (nproc > 1) call commrecdirdip(nrhs, 2, diprec, dip, buffer_send, buffer_recv, &
      &ticket%dip_recv_req, ticket%dip_send_req)
   call reset_polar_ticket(ticket)
   call report_polar_finish(status_value)
end subroutine polar_finish_dip_legacy

module subroutine reset_polar_ticket(ticket)
   type(polar_comm_t), intent(inout) :: ticket

   if (allocated(ticket%recv_req)) deallocate(ticket%recv_req)
   if (allocated(ticket%send_req)) deallocate(ticket%send_req)
   if (allocated(ticket%dip_recv_req)) deallocate(ticket%dip_recv_req)
   if (allocated(ticket%dip_send_req)) deallocate(ticket%dip_send_req)
   if (allocated(ticket%recdir_recv_req)) deallocate(ticket%recdir_recv_req)
   if (allocated(ticket%recdir_send_req)) deallocate(ticket%recdir_send_req)
   if (ticket%active .or. ticket%dip_active .or. ticket%recdir_active) then
      if (active_polar_ticket_count > 0) then
         active_polar_ticket_count = active_polar_ticket_count - 1
         call record_lifecycle_event(LifecycleTicket)
      end if
   end if
   ticket%active = .false.
   ticket%dip_active = .false.
   ticket%recdir_active = .false.
   ticket%recdir_send_phase = .false.
   ticket%recdir_receive_phase = .false.
   ticket%use_short_range = .false.
   ticket%use_group = .false.
   ticket%use_gpu = .false.
end subroutine reset_polar_ticket

module subroutine polar_post_recdir(ticket, nrhs, values, receive_buffer, send_buffer, channel, use_gpu, send_phase)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(in) :: values(3,nrhs,*)
   real(t_p), intent(inout) :: receive_buffer(3,nrhs,*), send_buffer(3,nrhs,*)
   type(channel_t), intent(in) :: channel
   logical, intent(in), optional :: use_gpu
   logical, intent(in), optional :: send_phase
   logical :: split_layout, gpu_path, send_path, receive_path
   integer :: i, j, k, l, ierr, peer, tag, lenbuf, ibufbeg

   gpu_path = .false.
   if (present(use_gpu)) gpu_path = use_gpu
   split_layout = use_pmecore .and. ndir > 0 .and. ndir < nproc
   send_path = .true.
   receive_path = .true.
   if (present(send_phase)) then
      send_path = send_phase
      receive_path = .not.send_phase
   else if (split_layout) then
      send_path = rank >= ndir
      receive_path = rank < ndir
   end if
   if (ticket%active .or. ticket%dip_active .or. ticket%recdir_active) &
      &error stop 'polar exchange ticket is already active'
   allocate(ticket%recdir_recv_req(max(1,nproc)), ticket%recdir_send_req(max(1,nproc)))
   ticket%recdir_recv_req = MPI_REQUEST_NULL
   ticket%recdir_send_req = MPI_REQUEST_NULL
   ticket%recdir_active = .true.
   ticket%recdir_send_phase = send_path
   ticket%recdir_receive_phase = receive_path
   active_polar_ticket_count = active_polar_ticket_count + 1
   if (nproc.eq.1) return

   call check_domain_layout()
   if (receive_path) then
      if (gpu_path) then
!$acc wait(channel%queue)
!$acc host_data use_device(receive_buffer)
         do i = 1, nrecdir_send1
            peer = precdir_send1(i)
            if (peer.ne.rank) then
               tag = nproc*rank + peer + 1
               call MPI_IRECV(receive_buffer(1,1,bufbeg1(peer+1)), 3*nrhs*buflen1(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_recv_req(i), ierr)
            end if
         end do
!$acc end host_data
      else
         do i = 1, nrecdir_send1
            peer = precdir_send1(i)
            if (peer.ne.rank) then
               tag = nproc*rank + peer + 1
               call MPI_IRECV(receive_buffer(1,1,bufbeg1(peer+1)), 3*nrhs*buflen1(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_recv_req(i), ierr)
            end if
         end do
      end if
   end if
   if (send_path) then
      if (gpu_path) then
         do i = 1, nrecdir_recep1
            peer = precdir_recep1(i)
            if (peer.ne.rank) then
               lenbuf = buflen2(peer+1)-1
               ibufbeg = bufbeg2(peer+1)
!$acc parallel loop gang collapse(3) async(channel%queue) default(present)
               do j = 0, lenbuf
                  do k = 1, nrhs
                     do l = 1, 3
                        send_buffer(l,k,ibufbeg+j) = values(l,k,polerecloc(buf2(ibufbeg+j)))
                     end do
                  end do
               end do
            end if
         end do
!$acc wait(channel%queue)
!$acc host_data use_device(send_buffer)
         do i = 1, nrecdir_recep1
            peer = precdir_recep1(i)
            if (peer.ne.rank) then
               tag = nproc*peer + rank + 1
               call MPI_ISEND(send_buffer(1,1,bufbeg2(peer+1)), 3*nrhs*buflen2(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_send_req(i), ierr)
            end if
         end do
!$acc end host_data
      else
         do i = 1, nrecdir_recep1
            peer = precdir_recep1(i)
            if (peer.ne.rank) then
               do j = 0, buflen2(peer+1)-1
                  call amove(3*nrhs, values(1,1,polerecloc(buf2(bufbeg2(peer+1)+j))), &
                     &send_buffer(1,1,bufbeg2(peer+1)+j))
               end do
               tag = nproc*peer + rank + 1
               call MPI_ISEND(send_buffer(1,1,bufbeg2(peer+1)), 3*nrhs*buflen2(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_send_req(i), ierr)
            end if
         end do
      end if
   end if
end subroutine polar_post_recdir

module subroutine polar_finish_recdir(ticket, nrhs, values, buffer, channel, use_gpu)
   type(polar_comm_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   real(t_p), intent(in) :: buffer(3,nrhs,*)
   type(channel_t), intent(in) :: channel
   logical, intent(in), optional :: use_gpu
   logical :: gpu_path
   integer :: i, j, k, l, ierr, peer, lenbuf, ibufbeg
   integer :: status(MPI_STATUS_SIZE)

   gpu_path = .false.
   if (present(use_gpu)) gpu_path = use_gpu
   if (.not.ticket%recdir_active) error stop 'polar exchange ticket is not active'
   if (nproc.eq.1) then
      call reset_polar_ticket(ticket)
      return
   end if
   call check_domain_layout()
   if (ticket%recdir_send_phase) then
      do i = 1, nrecdir_recep1
         peer = precdir_recep1(i)
         if (peer.ne.rank) call MPI_WAIT(ticket%recdir_send_req(i), status, ierr)
      end do
   end if
   if (ticket%recdir_receive_phase) then
      do i = 1, nrecdir_send1
         peer = precdir_send1(i)
         if (peer.ne.rank) then
            call MPI_WAIT(ticket%recdir_recv_req(i), status, ierr)
            lenbuf = buflen1(peer+1)-1
            ibufbeg = bufbeg1(peer+1)
            if (gpu_path) then
!$acc parallel loop gang collapse(3) async(channel%queue) default(present)
               do j = 0, lenbuf
                  do k = 1, nrhs
                     do l = 1, 3
                        values(l,k,poleloc(buf1(ibufbeg+j))) = buffer(l,k,ibufbeg+j)
                     end do
                  end do
               end do
            else
               do j = 0, lenbuf
                  call amove(3*nrhs, buffer(1,1,ibufbeg+j), values(1,1,poleloc(buf1(ibufbeg+j))))
               end do
            end if
         end if
      end do
   end if
   if (gpu_path) then
!$acc wait(channel%queue)
   end if
   call reset_polar_ticket(ticket)
end subroutine polar_finish_recdir

module subroutine polar_post_recdir_fields(ticket, values, receive_buffer, send_buffer, channel, use_gpu, send_phase)
   type(polar_comm_t), intent(inout) :: ticket
   real(t_p), intent(in) :: values(10,*)
   real(t_p), intent(inout) :: receive_buffer(10,*), send_buffer(10,*)
   type(channel_t), intent(in) :: channel
   logical, intent(in) :: use_gpu
   logical, intent(in), optional :: send_phase
   logical :: split_layout, send_path, receive_path
   integer :: i, j, k, ierr, peer, tag, lenbuf, ibufbeg

   split_layout = use_pmecore .and. ndir > 0 .and. ndir < nproc
   send_path = .true.
   receive_path = .true.
   if (present(send_phase)) then
      send_path = send_phase
      receive_path = .not.send_phase
   else if (split_layout) then
      send_path = rank >= ndir
      receive_path = rank < ndir
   end if
   if (ticket%active .or. ticket%dip_active .or. ticket%recdir_active) &
      &error stop 'polar exchange ticket is already active'
   allocate(ticket%recdir_recv_req(max(1,nproc)), ticket%recdir_send_req(max(1,nproc)))
   ticket%recdir_recv_req = MPI_REQUEST_NULL
   ticket%recdir_send_req = MPI_REQUEST_NULL
   ticket%recdir_active = .true.
   ticket%recdir_send_phase = send_path
   ticket%recdir_receive_phase = receive_path
   active_polar_ticket_count = active_polar_ticket_count + 1
   if (nproc.eq.1) return
   call check_domain_layout()
   if (receive_path) then
      if (use_gpu) then
!$acc wait(channel%queue)
!$acc host_data use_device(receive_buffer)
         do i = 1, nrecdir_send1
            peer = precdir_send1(i)
            if (peer.ne.rank) then
               tag = nproc*rank + peer + 1
               call MPI_IRECV(receive_buffer(1,bufbeg1(peer+1)), 10*buflen1(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_recv_req(i), ierr)
            end if
         end do
!$acc end host_data
      else
         do i = 1, nrecdir_send1
            peer = precdir_send1(i)
            if (peer.ne.rank) then
               tag = nproc*rank + peer + 1
               call MPI_IRECV(receive_buffer(1,bufbeg1(peer+1)), 10*buflen1(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_recv_req(i), ierr)
            end if
         end do
      end if
   end if
   if (send_path) then
      if (use_gpu) then
         do i = 1, nrecdir_recep1
            peer = precdir_recep1(i)
            if (peer.ne.rank) then
               lenbuf = buflen2(peer+1)-1
               ibufbeg = bufbeg2(peer+1)
!$acc parallel loop gang collapse(2) async(channel%queue) default(present)
               do j = 0, lenbuf
                  do k = 1, 10
                     send_buffer(k,ibufbeg+j) = values(k,polerecloc(buf2(ibufbeg+j)))
                  end do
               end do
            end if
         end do
!$acc wait(channel%queue)
!$acc host_data use_device(send_buffer)
         do i = 1, nrecdir_recep1
            peer = precdir_recep1(i)
            if (peer.ne.rank) then
               tag = nproc*peer + rank + 1
               call MPI_ISEND(send_buffer(1,bufbeg2(peer+1)), 10*buflen2(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_send_req(i), ierr)
            end if
         end do
!$acc end host_data
      else
         do i = 1, nrecdir_recep1
            peer = precdir_recep1(i)
            if (peer.ne.rank) then
               do j = 0, buflen2(peer+1)-1
                  call amove(10, values(1,polerecloc(buf2(bufbeg2(peer+1)+j))), &
                     &send_buffer(1,bufbeg2(peer+1)+j))
               end do
               tag = nproc*peer + rank + 1
               call MPI_ISEND(send_buffer(1,bufbeg2(peer+1)), 10*buflen2(peer+1), MPI_TPREC, peer, tag, &
                  &get_domain_comm(), ticket%recdir_send_req(i), ierr)
            end if
         end do
      end if
   end if
end subroutine polar_post_recdir_fields

module subroutine polar_finish_recdir_fields(ticket, values, buffer, channel, use_gpu)
   type(polar_comm_t), intent(inout) :: ticket
   real(t_p), intent(inout) :: values(10,*)
   real(t_p), intent(in) :: buffer(10,*)
   type(channel_t), intent(in) :: channel
   logical, intent(in) :: use_gpu
   integer :: i, j, k, ierr, peer, lenbuf, ibufbeg
   integer :: status(MPI_STATUS_SIZE)

   if (.not.ticket%recdir_active) error stop 'polar exchange ticket is not active'
   if (nproc.eq.1) then
      call reset_polar_ticket(ticket)
      return
   end if
   call check_domain_layout()
   if (ticket%recdir_send_phase) then
      do i = 1, nrecdir_recep1
         peer = precdir_recep1(i)
         if (peer.ne.rank) call MPI_WAIT(ticket%recdir_send_req(i), status, ierr)
      end do
   end if
   if (ticket%recdir_receive_phase) then
      do i = 1, nrecdir_send1
         peer = precdir_send1(i)
         if (peer.ne.rank) then
            call MPI_WAIT(ticket%recdir_recv_req(i), status, ierr)
            lenbuf = buflen1(peer+1)-1
            ibufbeg = bufbeg1(peer+1)
            if (use_gpu) then
!$acc parallel loop gang collapse(2) async(channel%queue) default(present)
               do j = 0, lenbuf
                  do k = 1, 10
                     values(k,poleloc(buf1(ibufbeg+j))) = buffer(k,ibufbeg+j)
                  end do
               end do
            else
               do j = 0, lenbuf
                  call amove(10, buffer(1,ibufbeg+j), values(1,poleloc(buf1(ibufbeg+j))))
               end do
            end if
         end if
      end do
   end if
   if (use_gpu) then
!$acc wait(channel%queue)
   end if
   call reset_polar_ticket(ticket)
end subroutine polar_finish_recdir_fields

subroutine post_group_comm(nrhs, values, recv_req, send_req)
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(3,nrhs,*)
   integer, intent(inout) :: recv_req(:), send_req(:)
   integer :: i, ierr, tag

   if (nproc.eq.1 .or. (use_pmecore .and. rank.gt.ndir-1)) return

!$acc host_data use_device(values)
   do i = 1, nproc
      if (i.ne.rank+1) then
         tag = nproc*rank + i
         call MPI_IRECV(values(1,1,bufbegpolegroup(i)), &
            &3*nrhs*domlenpolegroup(i),MPI_TPREC,i-1,tag,get_domain_comm(),recv_req(i),ierr)
      end if
   end do
   do i = 1, nproc
      if (i.ne.rank+1) then
         tag = nproc*(i-1) + rank + 1
         call MPI_ISEND(values,3*nrhs*npolelocgroup,MPI_TPREC,i-1,tag,get_domain_comm(),send_req(i),ierr)
      end if
   end do
!$acc end host_data
end subroutine post_group_comm

subroutine finish_group_comm(recv_req, send_req)
   integer, intent(inout) :: recv_req(:), send_req(:)
   integer :: i, ierr, status(MPI_STATUS_SIZE)

   if (nproc.eq.1 .or. (use_pmecore .and. rank.gt.ndir-1)) return
   do i = 1, nproc
      if (i.ne.rank+1) call MPI_WAIT(recv_req(i),status,ierr)
   end do
   do i = 1, nproc
      if (i.ne.rank+1) call MPI_WAIT(send_req(i),status,ierr)
   end do
end subroutine finish_group_comm

module function active_polar_tickets() result(count)
   integer :: count

   count = active_polar_ticket_count
end function active_polar_tickets

module subroutine domain_workspace_i(slot, required, values)
   integer, intent(in) :: slot, required
   integer, pointer, intent(out) :: values(:)

   call check_domain_layout()
!$acc wait
   select case (slot)
   case (1)
      if (.not.allocated(active_domain%pool_i1) .or. size(active_domain%pool_i1) < required) then
         if (allocated(active_domain%pool_i1)) then
!$acc wait
!$acc exit data delete(active_domain%pool_i1(1:size(active_domain%pool_i1)))
            deallocate(active_domain%pool_i1)
         end if
         allocate(active_domain%pool_i1(max(1,required)))
!$acc enter data create(active_domain%pool_i1(1:max(1,required)))
      end if
      values => active_domain%pool_i1
   case (2)
      if (.not.allocated(active_domain%pool_i2) .or. size(active_domain%pool_i2) < required) then
         if (allocated(active_domain%pool_i2)) then
!$acc wait
!$acc exit data delete(active_domain%pool_i2(1:size(active_domain%pool_i2)))
            deallocate(active_domain%pool_i2)
         end if
         allocate(active_domain%pool_i2(max(1,required)))
!$acc enter data create(active_domain%pool_i2(1:max(1,required)))
      end if
      values => active_domain%pool_i2
   case default
      error stop 'unknown integer communication workspace'
   end select
end subroutine domain_workspace_i

module subroutine domain_workspace_rp(slot, required, values)
   integer, intent(in) :: slot, required
   real(r_p), pointer, intent(out) :: values(:)

   call check_domain_layout()
!$acc wait
   select case (slot)
   case (1)
      if (.not.allocated(active_domain%pool_rp1) .or. size(active_domain%pool_rp1) < required) then
         if (allocated(active_domain%pool_rp1)) then
!$acc wait
!$acc exit data delete(active_domain%pool_rp1(1:size(active_domain%pool_rp1)))
            deallocate(active_domain%pool_rp1)
         end if
         allocate(active_domain%pool_rp1(max(1,required)))
!$acc enter data create(active_domain%pool_rp1(1:max(1,required)))
      end if
      values => active_domain%pool_rp1
   case (2)
      if (.not.allocated(active_domain%pool_rp2) .or. size(active_domain%pool_rp2) < required) then
         if (allocated(active_domain%pool_rp2)) then
!$acc wait
!$acc exit data delete(active_domain%pool_rp2(1:size(active_domain%pool_rp2)))
            deallocate(active_domain%pool_rp2)
         end if
         allocate(active_domain%pool_rp2(max(1,required)))
!$acc enter data create(active_domain%pool_rp2(1:max(1,required)))
      end if
      values => active_domain%pool_rp2
   case (3)
      if (.not.allocated(active_domain%pool_rp3) .or. size(active_domain%pool_rp3) < required) then
         if (allocated(active_domain%pool_rp3)) then
!$acc wait
!$acc exit data delete(active_domain%pool_rp3(1:size(active_domain%pool_rp3)))
            deallocate(active_domain%pool_rp3)
         end if
         allocate(active_domain%pool_rp3(max(1,required)))
!$acc enter data create(active_domain%pool_rp3(1:max(1,required)))
      end if
      values => active_domain%pool_rp3
   case default
      error stop 'unknown real communication workspace'
   end select
end subroutine domain_workspace_rp

#if defined(USE_DETERMINISTIC_REDUCTION)
module subroutine domain_workspace_md(required, values)
   integer, intent(in) :: required
   mdyn_rtyp, pointer, intent(out) :: values(:)

   call check_domain_layout()
!$acc wait
   if (.not.allocated(active_domain%work_md) .or. size(active_domain%work_md) < required) then
      if (allocated(active_domain%work_md)) then
!$acc exit data delete(active_domain%work_md(1:size(active_domain%work_md)))
         deallocate(active_domain%work_md)
      end if
      allocate(active_domain%work_md(max(1,required)))
!$acc enter data create(active_domain%work_md(1:max(1,required)))
   end if
   values => active_domain%work_md
end subroutine domain_workspace_md
#endif

module subroutine domain_workspace_tp2d(slot, n1, n2, values)
   integer, intent(in) :: slot, n1, n2
   real(t_p), pointer, intent(out) :: values(:,:)
   real(t_p), pointer :: raw(:)
   integer required

   call check_domain_layout()
!$acc wait
   required = max(1,n1)*max(1,n2)
   select case (slot)
   case (1)
      if (.not.allocated(active_domain%pool_t5) .or. size(active_domain%pool_t5) < required) then
         if (allocated(active_domain%pool_t5)) then
!$acc wait
!$acc exit data delete(active_domain%pool_t5(1:size(active_domain%pool_t5)))
            deallocate(active_domain%pool_t5)
         end if
         allocate(active_domain%pool_t5(required))
!$acc enter data create(active_domain%pool_t5(1:required))
      end if
      raw => active_domain%pool_t5
   case (2)
      if (.not.allocated(active_domain%pool_t6) .or. size(active_domain%pool_t6) < required) then
         if (allocated(active_domain%pool_t6)) then
!$acc wait
!$acc exit data delete(active_domain%pool_t6(1:size(active_domain%pool_t6)))
            deallocate(active_domain%pool_t6)
         end if
         allocate(active_domain%pool_t6(required))
!$acc enter data create(active_domain%pool_t6(1:required))
      end if
      raw => active_domain%pool_t6
   case default
      error stop 'unknown 2D real communication workspace'
   end select
   values(1:max(1,n1),1:max(1,n2)) => raw(1:required)
end subroutine domain_workspace_tp2d

module subroutine domain_workspace_tp3d(slot, n1, n2, n3, values)
   integer, intent(in) :: slot, n1, n2, n3
   real(t_p), pointer, intent(out) :: values(:,:,:)
   real(t_p), pointer :: raw(:)
   integer required

   call check_domain_layout()
!$acc wait
   required = max(1,n1)*max(1,n2)*max(1,n3)
   select case (slot)
   case (1)
      if (.not.allocated(active_domain%pool_t1) .or. size(active_domain%pool_t1) < required) then
         if (allocated(active_domain%pool_t1)) then
!$acc wait
!$acc exit data delete(active_domain%pool_t1(1:size(active_domain%pool_t1)))
            deallocate(active_domain%pool_t1)
         end if
         allocate(active_domain%pool_t1(required))
!$acc enter data create(active_domain%pool_t1(1:required))
      end if
      raw => active_domain%pool_t1
   case (2)
      if (.not.allocated(active_domain%pool_t2) .or. size(active_domain%pool_t2) < required) then
         if (allocated(active_domain%pool_t2)) then
!$acc wait
!$acc exit data delete(active_domain%pool_t2(1:size(active_domain%pool_t2)))
            deallocate(active_domain%pool_t2)
         end if
         allocate(active_domain%pool_t2(required))
!$acc enter data create(active_domain%pool_t2(1:required))
      end if
      raw => active_domain%pool_t2
   case (3)
      if (.not.allocated(active_domain%pool_t3) .or. size(active_domain%pool_t3) < required) then
         if (allocated(active_domain%pool_t3)) then
!$acc wait
!$acc exit data delete(active_domain%pool_t3(1:size(active_domain%pool_t3)))
            deallocate(active_domain%pool_t3)
         end if
         allocate(active_domain%pool_t3(required))
!$acc enter data create(active_domain%pool_t3(1:required))
      end if
      raw => active_domain%pool_t3
   case (4)
      if (.not.allocated(active_domain%pool_t4) .or. size(active_domain%pool_t4) < required) then
         if (allocated(active_domain%pool_t4)) then
!$acc wait
!$acc exit data delete(active_domain%pool_t4(1:size(active_domain%pool_t4)))
            deallocate(active_domain%pool_t4)
         end if
         allocate(active_domain%pool_t4(required))
!$acc enter data create(active_domain%pool_t4(1:required))
      end if
      raw => active_domain%pool_t4
   case default
      error stop 'unknown 3D real communication workspace'
   end select
   values(1:max(1,n1),1:max(1,n2),1:max(1,n3)) => raw(1:required)
end subroutine domain_workspace_tp3d

module subroutine capture_halo_layout(domain)
   type(domain_t), intent(inout) :: domain

   call release_halo_layout(domain)
   domain%nproc = nproc
   domain%nloc = nloc
   domain%nbig_send = nbig_send
   domain%nbig_recep = nbig_recep
   domain%nneig_send = nneig_send
   domain%nneig_recep = nneig_recep
   domain%nshort_send = nbigshort_send
   domain%nshort_recep = nbigshort_recep
   allocate(domain%domlen(nproc), domain%bufbeg(nproc))
   allocate(domain%pbig_send(max(1,nproc)), domain%pbig_recep(max(1,nproc)))
   allocate(domain%pneig_send(max(1,nproc)), domain%pneig_recep(max(1,nproc)))
   allocate(domain%pshort_send(max(1,nproc)), domain%pshort_recep(max(1,nproc)))
   allocate(domain%reqs(max(1,nproc)), domain%reqr(max(1,nproc)))
   domain%domlen = domlen
   domain%bufbeg = bufbeg
   domain%pbig_send = pbig_send
   domain%pbig_recep = pbig_recep
   domain%pneig_send = pneig_send
   domain%pneig_recep = pneig_recep
   domain%pshort_send = pbigshort_send
   domain%pshort_recep = pbigshort_recep
   domain%reqs = MPI_REQUEST_NULL
   domain%reqr = MPI_REQUEST_NULL
   domain%layout_ready = .true.
end subroutine capture_halo_layout

module subroutine refresh_domain_layout(domain)
   type(domain_t), intent(inout) :: domain

   call capture_halo_layout(domain)
end subroutine refresh_domain_layout

module subroutine release_halo_layout(domain, cleanup_status)
   type(domain_t), intent(inout) :: domain
   integer, intent(out), optional :: cleanup_status
   integer :: ierr, first_error

!$acc wait
   first_error = MPI_SUCCESS
   if (allocated(domain%reqs)) then
      call MPI_Waitall(size(domain%reqs), domain%reqs, MPI_STATUSES_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
   end if
   if (allocated(domain%reqr)) then
      call MPI_Waitall(size(domain%reqr), domain%reqr, MPI_STATUSES_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
   end if
   if (allocated(domain%recdir_recv_reqs)) then
      call MPI_Waitall(size(domain%recdir_recv_reqs), domain%recdir_recv_reqs, &
         &MPI_STATUSES_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
   end if
   if (allocated(domain%recdir_send_reqs)) then
      call MPI_Waitall(size(domain%recdir_send_reqs), domain%recdir_send_reqs, &
         &MPI_STATUSES_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
   end if
   if (allocated(domain%work_tp)) then
!$acc exit data delete(domain%work_tp(1:size(domain%work_tp)))
      deallocate(domain%work_tp)
   end if
   if (allocated(domain%work_rp)) then
!$acc exit data delete(domain%work_rp(1:size(domain%work_rp)))
      deallocate(domain%work_rp)
   end if
#if defined(USE_DETERMINISTIC_REDUCTION)
   if (allocated(domain%work_md)) then
!$acc exit data delete(domain%work_md(1:size(domain%work_md)))
      deallocate(domain%work_md)
   end if
#endif
   if (allocated(domain%recdir_receive_work)) then
!$acc exit data delete(domain%recdir_receive_work(1:size(domain%recdir_receive_work)))
      deallocate(domain%recdir_receive_work)
   end if
   if (allocated(domain%recdir_send_work)) then
!$acc exit data delete(domain%recdir_send_work(1:size(domain%recdir_send_work)))
      deallocate(domain%recdir_send_work)
   end if
   if (allocated(domain%recdir_recv_reqs)) deallocate(domain%recdir_recv_reqs)
   if (allocated(domain%recdir_send_reqs)) deallocate(domain%recdir_send_reqs)
   if (allocated(domain%pool_i1)) then
!$acc exit data delete(domain%pool_i1(1:size(domain%pool_i1)))
      deallocate(domain%pool_i1)
   end if
   if (allocated(domain%pool_i2)) then
!$acc exit data delete(domain%pool_i2(1:size(domain%pool_i2)))
      deallocate(domain%pool_i2)
   end if
   if (allocated(domain%pool_rp1)) then
!$acc exit data delete(domain%pool_rp1(1:size(domain%pool_rp1)))
      deallocate(domain%pool_rp1)
   end if
   if (allocated(domain%pool_rp2)) then
!$acc exit data delete(domain%pool_rp2(1:size(domain%pool_rp2)))
      deallocate(domain%pool_rp2)
   end if
   if (allocated(domain%pool_rp3)) then
!$acc exit data delete(domain%pool_rp3(1:size(domain%pool_rp3)))
      deallocate(domain%pool_rp3)
   end if
   if (allocated(domain%pool_t1)) then
!$acc exit data delete(domain%pool_t1(1:size(domain%pool_t1)))
      deallocate(domain%pool_t1)
   end if
   if (allocated(domain%pool_t2)) then
!$acc exit data delete(domain%pool_t2(1:size(domain%pool_t2)))
      deallocate(domain%pool_t2)
   end if
   if (allocated(domain%pool_t3)) then
!$acc exit data delete(domain%pool_t3(1:size(domain%pool_t3)))
      deallocate(domain%pool_t3)
   end if
   if (allocated(domain%pool_t4)) then
!$acc exit data delete(domain%pool_t4(1:size(domain%pool_t4)))
      deallocate(domain%pool_t4)
   end if
   if (allocated(domain%pool_t5)) then
!$acc exit data delete(domain%pool_t5(1:size(domain%pool_t5)))
      deallocate(domain%pool_t5)
   end if
   if (allocated(domain%pool_t6)) then
!$acc exit data delete(domain%pool_t6(1:size(domain%pool_t6)))
      deallocate(domain%pool_t6)
   end if
   if (allocated(domain%domlen)) deallocate(domain%domlen)
   if (allocated(domain%bufbeg)) deallocate(domain%bufbeg)
   if (allocated(domain%pbig_send)) deallocate(domain%pbig_send)
   if (allocated(domain%pbig_recep)) deallocate(domain%pbig_recep)
   if (allocated(domain%pneig_send)) deallocate(domain%pneig_send)
   if (allocated(domain%pneig_recep)) deallocate(domain%pneig_recep)
   if (allocated(domain%pshort_send)) deallocate(domain%pshort_send)
   if (allocated(domain%pshort_recep)) deallocate(domain%pshort_recep)
   if (allocated(domain%reqs)) deallocate(domain%reqs)
   if (allocated(domain%reqr)) deallocate(domain%reqr)
   domain%nproc = 0
   domain%nloc = 0
   domain%nshort_send = 0
   domain%nshort_recep = 0
   domain%layout_ready = .false.
   if (present(cleanup_status)) cleanup_status = first_error
end subroutine release_halo_layout

module subroutine extend_legacy_tp(values, components, halo, channel)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel

   call exchange_extend_tp(values, components, halo, channel)
end subroutine extend_legacy_tp

module subroutine accumulate_legacy_tp(values, components, halo, channel)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel

   call exchange_accumulate_tp(values, components, halo, channel)
end subroutine accumulate_legacy_tp

#if TINKER_MIXED_PREC
module subroutine extend_legacy_rp(values, components, halo, channel)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel

   call exchange_extend_rp(values, components, halo, channel)
end subroutine extend_legacy_rp

module subroutine accumulate_legacy_rp(values, components, halo, channel)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel

   call exchange_accumulate_rp(values, components, halo, channel)
end subroutine accumulate_legacy_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
module subroutine extend_legacy_md(values, components, halo, channel)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel

   call exchange_extend_md(values, components, halo, channel)
end subroutine extend_legacy_md

module subroutine accumulate_legacy_md(values, components, halo, channel)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel

   call exchange_accumulate_md(values, components, halo, channel)
end subroutine accumulate_legacy_md
#endif

module subroutine accumulate_rec_to_dir_legacy_tp(direct_values, reciprocal_values, components, channel)
   real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
   integer, intent(in) :: components
   type(channel_t), intent(in) :: channel

   call exchange_accumulate_rec_to_dir_tp(direct_values, reciprocal_values, components, channel)
end subroutine accumulate_rec_to_dir_legacy_tp

module subroutine replicate_legacy_rp(global_ids, owned_values, global_values, components)
   integer, intent(in) :: global_ids(:)
   real(r_p), intent(in) :: owned_values(:)
   real(r_p), intent(out) :: global_values(:)
   integer, intent(in) :: components
   integer :: ierr, i, local_count, max_global_id, value_count, total_count
   integer, allocatable :: counts(:), offsets(:)
   integer, pointer :: receive_ids(:)
   real(r_p), pointer :: received_values(:)

   call check_domain_layout()
   local_count = merge(active_domain%nloc, 0, .not.(use_pmecore .and. active_domain%rank > ndir-1))
   if (components <= 0 .or. size(global_ids) < local_count .or. &
      &size(owned_values) < components*local_count) then
      error stop 'invalid replication payload'
   end if
   allocate(counts(active_domain%nproc), offsets(active_domain%nproc))
   call MPI_Allgather(local_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, active_domain%communicator, ierr)
   offsets(1) = 0
   do i = 2, active_domain%nproc
      offsets(i) = offsets(i-1) + counts(i-1)
   end do
   total_count = sum(counts)
   call domain_workspace_i(1, total_count, receive_ids)
   call domain_workspace_rp(3, components*total_count, received_values)

!$acc data present_or_copyin(global_ids,owned_values) present_or_copyout(global_values)
!$acc wait
!$acc host_data use_device(global_ids,receive_ids)
   call MPI_Allgatherv(global_ids, local_count, MPI_INTEGER, receive_ids, counts, offsets, &
      &MPI_INTEGER, active_domain%communicator, ierr)
!$acc end host_data
   counts = components*counts
   offsets = components*offsets
!$acc host_data use_device(owned_values,received_values)
   call MPI_Allgatherv(owned_values, components*local_count, MPI_RPREC, received_values, counts, offsets, &
      &MPI_RPREC, active_domain%communicator, ierr)
!$acc end host_data

   max_global_id = 0
   if (total_count > 0) then
!$acc update self(receive_ids(1:total_count))
      max_global_id = maxval(receive_ids(1:total_count))
   end if
   if (size(global_values) < components*max_global_id) error stop 'replication result is too small'

!$acc parallel loop default(present)
   do i = 1, size(global_values)
      global_values(i) = real(0.0,r_p)
   end do
!$acc end parallel loop

   do i = 1, total_count
      if (receive_ids(i) <= 0 .or. receive_ids(i) > size(global_values)/components) then
         error stop 'replication received an invalid global identifier'
      end if
   end do
!$acc parallel loop default(present)
   do i = 1, total_count
      value_count = components*(receive_ids(i)-1)
      global_values(value_count+1:value_count+components) = &
         &received_values(components*(i-1)+1:components*i)
   end do
!$acc end parallel loop
!$acc end data
end subroutine replicate_legacy_rp

module subroutine replicate_full_rp(global_ids, global_values, components)
   integer, intent(in) :: global_ids(:)
   real(r_p), intent(inout) :: global_values(:)
   integer, intent(in) :: components
   integer :: i, component, atom_id, local_count
   real(r_p), pointer :: owned_values(:)

   call check_domain_layout()
   if (active_domain%nproc == 1) return
   local_count = merge(active_domain%nloc, 0, .not.(use_pmecore .and. active_domain%rank > ndir-1))
   if (components <= 0 .or. size(global_ids) < local_count) error stop 'invalid replication payload'
   call domain_workspace_rp(1, components*local_count, owned_values)

!$acc data present_or_copyin(global_ids) present_or_copy(global_values)
!$acc parallel loop collapse(2) async default(present) private(atom_id)
   do i = 1, local_count
      do component = 1, components
         atom_id = global_ids(i)
         if (atom_id > 0) then
            if (atom_id <= size(global_values)/components) then
               owned_values(components*(i-1)+component) = global_values(components*(atom_id-1)+component)
            else
               owned_values(components*(i-1)+component) = real(0.0,r_p)
            end if
         else
            owned_values(components*(i-1)+component) = real(0.0,r_p)
         end if
      end do
   end do
!$acc end parallel loop
!$acc wait
   call replicate_rp(global_ids, owned_values, global_values, components)
!$acc end data
end subroutine replicate_full_rp

#if TINKER_MIXED_PREC
module subroutine reduce_sum_legacy_tp_vector(values, visibility, root)
   real(t_p), intent(inout) :: values(:)
   integer, intent(in) :: visibility
   integer, intent(in), optional :: root
   integer :: ierr, root_rank, communicator, process_rank, process_count

   call check_reduction(visibility, root, root_rank, communicator, process_rank, process_count)
   if (process_count == 1) return
   if (visibility == reduction_all) then
      call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_TPREC, MPI_SUM, communicator, ierr)
   else if (process_rank == root_rank) then
      call MPI_Reduce(MPI_IN_PLACE, values, size(values), MPI_TPREC, MPI_SUM, root_rank, communicator, ierr)
   else
      call MPI_Reduce(values, values, size(values), MPI_TPREC, MPI_SUM, root_rank, communicator, ierr)
   end if
end subroutine reduce_sum_legacy_tp_vector
#endif

module subroutine reduce_sum_legacy_i_vector(values, visibility, root)
   integer, intent(inout) :: values(:)
   integer, intent(in) :: visibility
   integer, intent(in), optional :: root
   integer :: ierr, root_rank, communicator, process_rank, process_count

   call check_reduction(visibility, root, root_rank, communicator, process_rank, process_count)
   if (process_count == 1) return
   if (visibility == reduction_all) then
      call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_INTEGER, MPI_SUM, communicator, ierr)
   else if (process_rank == root_rank) then
      call MPI_Reduce(MPI_IN_PLACE, values, size(values), MPI_INTEGER, MPI_SUM, root_rank, communicator, ierr)
   else
      call MPI_Reduce(values, values, size(values), MPI_INTEGER, MPI_SUM, root_rank, communicator, ierr)
   end if
end subroutine reduce_sum_legacy_i_vector

module subroutine reduce_sum_legacy_rp_vector(values, visibility, root)
   real(r_p), intent(inout) :: values(:)
   integer, intent(in) :: visibility
   integer, intent(in), optional :: root
   integer :: ierr, root_rank, communicator, process_rank, process_count

   call check_reduction(visibility, root, root_rank, communicator, process_rank, process_count)
   if (process_count == 1) return
   if (visibility == reduction_all) then
      call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_RPREC, MPI_SUM, communicator, ierr)
   else if (process_rank == root_rank) then
      call MPI_Reduce(MPI_IN_PLACE, values, size(values), MPI_RPREC, MPI_SUM, root_rank, communicator, ierr)
   else
      call MPI_Reduce(values, values, size(values), MPI_RPREC, MPI_SUM, root_rank, communicator, ierr)
   end if
end subroutine reduce_sum_legacy_rp_vector

module subroutine reduce_sum_legacy_rp_matrix(values, visibility, root)
   real(r_p), intent(inout) :: values(:,:)
   integer, intent(in) :: visibility
   integer, intent(in), optional :: root
   integer :: ierr, root_rank, communicator, process_rank, process_count

   call check_reduction(visibility, root, root_rank, communicator, process_rank, process_count)
   if (process_count == 1) return
   if (visibility == reduction_all) then
      call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_RPREC, MPI_SUM, communicator, ierr)
   else if (process_rank == root_rank) then
      call MPI_Reduce(MPI_IN_PLACE, values, size(values), MPI_RPREC, MPI_SUM, root_rank, communicator, ierr)
   else
      call MPI_Reduce(values, values, size(values), MPI_RPREC, MPI_SUM, root_rank, communicator, ierr)
   end if
end subroutine reduce_sum_legacy_rp_matrix

#if defined(USE_DETERMINISTIC_REDUCTION)
module subroutine reduce_sum_legacy_md_vector(values, visibility, root)
   mdyn_rtyp, intent(inout) :: values(:)
   integer, intent(in) :: visibility
   integer, intent(in), optional :: root
   integer :: ierr, root_rank, communicator, process_rank, process_count

   call check_reduction(visibility, root, root_rank, communicator, process_rank, process_count)
   if (process_count == 1) return
   if (visibility == reduction_all) then
      call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_MDTYP, MPI_SUM, communicator, ierr)
   else if (process_rank == root_rank) then
      call MPI_Reduce(MPI_IN_PLACE, values, size(values), MPI_MDTYP, MPI_SUM, root_rank, communicator, ierr)
   else
      call MPI_Reduce(values, values, size(values), MPI_MDTYP, MPI_SUM, root_rank, communicator, ierr)
   end if
end subroutine reduce_sum_legacy_md_vector
#endif

subroutine exchange_extend_tp(values, components, halo, channel)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel
   integer :: i, ierr, peer, first, count, option, comm, request
   integer :: status(MPI_STATUS_SIZE)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   comm = active_domain%communicator
!$acc wait(channel%queue)
   do i = 1, receive_count(option)
      peer = receive_peer(option, i)
      first = components*(active_domain%bufbeg(peer+1)-1)+1
      count = components*active_domain%domlen(peer+1)
      request = active_domain%reqr(i)
!$acc host_data use_device(values)
      call MPI_Irecv(values(first), count, MPI_TPREC, peer, 0, comm, request, ierr)
!$acc end host_data
      active_domain%reqr(i) = request
   end do
   do i = 1, send_count(option)
      peer = send_peer(option, i)
      request = active_domain%reqs(i)
!$acc host_data use_device(values)
      call MPI_Isend(values, components*active_domain%nloc, MPI_TPREC, peer, 0, comm, request, ierr)
!$acc end host_data
      active_domain%reqs(i) = request
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
   do i = 1, receive_count(option); call MPI_Wait(active_domain%reqr(i), status, ierr); end do
end subroutine exchange_extend_tp

subroutine exchange_accumulate_tp(values, components, halo, channel)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel
   integer :: i, j, ierr, option, incoming, peer, comm, request, value_count, first, count
   integer :: status(MPI_STATUS_SIZE)
   real(t_p), pointer :: work(:)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   value_count = components*active_domain%nloc
   comm = active_domain%communicator
   incoming = receive_count(option)
   call ensure_work_tp(value_count*max(1,incoming))
   work => active_domain%work_tp
!$acc wait(channel%queue)
!$acc host_data use_device(work,values)
   do i = 1, incoming
      peer = receive_peer(option, i)
      request = active_domain%reqr(i)
      call MPI_Irecv(work(1+components*active_domain%nloc*(i-1)), &
         &value_count, MPI_TPREC, peer, 0, comm, request, ierr)
      active_domain%reqr(i) = request
   end do
   do i = 1, send_count(option)
      peer = send_peer(option, i)
      first = components*(active_domain%bufbeg(peer+1)-1)+1
      count = components*active_domain%domlen(peer+1)
      request = active_domain%reqs(i)
      call MPI_Isend(values(first), count, MPI_TPREC, peer, 0, comm, request, ierr)
      active_domain%reqs(i) = request
   end do
!$acc end host_data
   do i = 1, incoming; call MPI_Wait(active_domain%reqr(i), status, ierr); end do
!$acc parallel loop gang vector async(channel%queue) present(values,work(1:max(1,value_count*incoming)))
   do j = 1, value_count
      do i = 1, incoming
         values(j) = values(j) + work(j+components*active_domain%nloc*(i-1))
      end do
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
end subroutine exchange_accumulate_tp

subroutine exchange_accumulate_rec_to_dir_tp(direct_values, reciprocal_values, components, channel)
   real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
   integer, intent(in) :: components
   type(channel_t), intent(in) :: channel
   integer :: i, j, k, ierr, peer, first, count, dimn, rankloc
   integer :: status(MPI_STATUS_SIZE)
   real(t_p), pointer :: receive_work(:), send_work(:)

   call check_domain_layout()
   if (active_domain%nproc == 1) then
!$acc parallel loop async(channel%queue) default(present)
      do j = 1, components* nlocrec
         direct_values(j) = direct_values(j) + reciprocal_values(j)
      end do
      return
   end if

   dimn = components*active_domain%nloc
   call ensure_recdir_work(components*max(1,active_domain%nloc)*nrecdir_send1, &
      &components*max(1,nblocrecdir))
   receive_work => active_domain%recdir_receive_work
   send_work => active_domain%recdir_send_work
   rankloc = merge(rank_bis, active_domain%rank, use_pmecore)

!$acc wait(channel%queue)
!$acc host_data use_device(receive_work)
   do i = 1, nrecdir_send1
      peer = precdir_send1(i)
      if (peer /= active_domain%rank) then
         call MPI_Irecv(receive_work(1+(i-1)*dimn), dimn, MPI_TPREC, peer, 0, &
            &active_domain%communicator, active_domain%recdir_recv_reqs(i), ierr)
      end if
   end do
!$acc end host_data

   do i = 1, nrecdir_recep1
      peer = precdir_recep1(i)
      if (peer /= active_domain%rank) then
         first = active_domain%bufbeg(peer+1)-1
         count = active_domain%domlen(peer+1)
!$acc parallel loop collapse(2) async(channel%queue) &
!$acc      present(repartrec,locrec,glob,reciprocal_values,send_work)
         do j = 1, count
            do k = 1, components
               if (repartrec(glob(first+j)) == rankloc) then
                  send_work(k+(first+j-1)*components) = &
                     &reciprocal_values(k+(locrec(glob(first+j))-1)*components)
               else
                  send_work(k+(first+j-1)*components) = real(0.0,t_p)
               end if
            end do
         end do
      end if
   end do
!$acc wait(channel%queue)
!$acc host_data use_device(send_work)
   do i = 1, nrecdir_recep1
      peer = precdir_recep1(i)
      if (peer /= active_domain%rank) then
         first = 1 + (active_domain%bufbeg(peer+1)-1)*components
         count = components*active_domain%domlen(peer+1)
         call MPI_Isend(send_work(first), count, MPI_TPREC, peer, 0, &
            &active_domain%communicator, active_domain%recdir_send_reqs(i), ierr)
      end if
   end do
!$acc end host_data

   do i = 1, nrecdir_send1
      peer = precdir_send1(i)
      if (peer /= active_domain%rank) then
         call MPI_Wait(active_domain%recdir_recv_reqs(i), status, ierr)
!$acc parallel loop async(channel%queue) default(present)
         do j = 1, dimn
            direct_values(j) = direct_values(j) + receive_work(j+(i-1)*dimn)
         end do
      else
!$acc parallel loop collapse(2) async(channel%queue) default(present)
         do j = 1, nlocrec
            do k = 1, components
               if (repart(globrec(j)) == active_domain%rank) then
                  direct_values(k+(loc(globrec(j))-1)*components) = &
                     &direct_values(k+(loc(globrec(j))-1)*components) + reciprocal_values(k+(j-1)*components)
               end if
            end do
         end do
      end if
   end do

   do i = 1, nrecdir_recep1
      peer = precdir_recep1(i)
      if (peer /= active_domain%rank) then
         call MPI_Wait(active_domain%recdir_send_reqs(i), status, ierr)
      end if
   end do
end subroutine exchange_accumulate_rec_to_dir_tp

#if TINKER_MIXED_PREC
subroutine exchange_extend_rp(values, components, halo, channel)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel
   integer :: i, ierr, peer, first, count, option, comm, request
   integer :: status(MPI_STATUS_SIZE)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   comm = active_domain%communicator
!$acc wait(channel%queue)
   do i = 1, receive_count(option)
      peer = receive_peer(option, i)
      first = components*(active_domain%bufbeg(peer+1)-1)+1
      count = components*active_domain%domlen(peer+1)
      request = active_domain%reqr(i)
!$acc host_data use_device(values)
      call MPI_Irecv(values(first), count, MPI_RPREC, peer, 0, comm, request, ierr)
!$acc end host_data
      active_domain%reqr(i) = request
   end do
   do i = 1, send_count(option)
      peer = send_peer(option, i)
      request = active_domain%reqs(i)
!$acc host_data use_device(values)
      call MPI_Isend(values, components*active_domain%nloc, MPI_RPREC, peer, 0, comm, request, ierr)
!$acc end host_data
      active_domain%reqs(i) = request
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
   do i = 1, receive_count(option); call MPI_Wait(active_domain%reqr(i), status, ierr); end do
end subroutine exchange_extend_rp

subroutine exchange_accumulate_rp(values, components, halo, channel)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel
   integer :: i, j, ierr, option, incoming, peer, comm, request, value_count
   integer :: first, count
   integer :: status(MPI_STATUS_SIZE)
   real(r_p), pointer :: work(:)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   value_count = components*active_domain%nloc
   comm = active_domain%communicator
   incoming = receive_count(option)
   call ensure_work_rp(value_count*max(1,incoming))
   work => active_domain%work_rp
!$acc wait(channel%queue)
!$acc host_data use_device(work,values)
   do i = 1, incoming
      peer = receive_peer(option, i)
      request = active_domain%reqr(i)
      call MPI_Irecv(work(1+components*active_domain%nloc*(i-1)), &
         &value_count, MPI_RPREC, peer, 0, comm, request, ierr)
      active_domain%reqr(i) = request
   end do
   do i = 1, send_count(option)
      peer = send_peer(option, i)
      first = components*(active_domain%bufbeg(peer+1)-1)+1
      count = components*active_domain%domlen(peer+1)
      request = active_domain%reqs(i)
      call MPI_Isend(values(first), count, MPI_RPREC, peer, 0, comm, request, ierr)
      active_domain%reqs(i) = request
   end do
!$acc end host_data
   do i = 1, incoming; call MPI_Wait(active_domain%reqr(i), status, ierr); end do
!$acc parallel loop gang vector async(channel%queue) present(values,work(1:max(1,value_count*incoming)))
   do j = 1, value_count
      do i = 1, incoming
         values(j) = values(j) + work(j+components*active_domain%nloc*(i-1))
      end do
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
end subroutine exchange_accumulate_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
subroutine exchange_extend_md(values, components, halo, channel)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel
   integer :: i, ierr, peer, first, count, option, comm, request
   integer :: status(MPI_STATUS_SIZE)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   comm = active_domain%communicator
!$acc wait(channel%queue)
   do i = 1, receive_count(option)
      peer = receive_peer(option, i)
      first = components*(active_domain%bufbeg(peer+1)-1)+1
      count = components*active_domain%domlen(peer+1)
      request = active_domain%reqr(i)
!$acc host_data use_device(values)
      call MPI_Irecv(values(first), count, MPI_MDTYP, peer, 0, comm, request, ierr)
!$acc end host_data
      active_domain%reqr(i) = request
   end do
   do i = 1, send_count(option)
      peer = send_peer(option, i)
      request = active_domain%reqs(i)
!$acc host_data use_device(values)
      call MPI_Isend(values, components*active_domain%nloc, MPI_MDTYP, peer, 0, comm, request, ierr)
!$acc end host_data
      active_domain%reqs(i) = request
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
   do i = 1, receive_count(option); call MPI_Wait(active_domain%reqr(i), status, ierr); end do
end subroutine exchange_extend_md

subroutine exchange_accumulate_md(values, components, halo, channel)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   type(channel_t), intent(in) :: channel
   integer :: i, j, ierr, option, incoming, peer, comm, request, value_count, first, count
   integer :: status(MPI_STATUS_SIZE)
   mdyn_rtyp, pointer :: work(:)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   value_count = components*active_domain%nloc
   comm = active_domain%communicator
   incoming = receive_count(option)
   call ensure_work_md(value_count*max(1,incoming))
   work => active_domain%work_md
!$acc wait(channel%queue)
!$acc host_data use_device(work,values)
   do i = 1, incoming
      peer = receive_peer(option, i)
      request = active_domain%reqr(i)
      call MPI_Irecv(work(1+components*active_domain%nloc*(i-1)), &
         &value_count, MPI_MDTYP, peer, 0, comm, request, ierr)
      active_domain%reqr(i) = request
   end do
   do i = 1, send_count(option)
      peer = send_peer(option, i)
      first = components*(active_domain%bufbeg(peer+1)-1)+1
      count = components*active_domain%domlen(peer+1)
      request = active_domain%reqs(i)
      call MPI_Isend(values(first), count, MPI_MDTYP, peer, 0, comm, request, ierr)
      active_domain%reqs(i) = request
   end do
!$acc end host_data
   do i = 1, incoming; call MPI_Wait(active_domain%reqr(i), status, ierr); end do
!$acc parallel loop gang vector async(channel%queue) present(values,work(1:max(1,value_count*incoming)))
   do j = 1, value_count
      do i = 1, incoming
         values(j) = values(j) + work(j+components*active_domain%nloc*(i-1))
      end do
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
end subroutine exchange_accumulate_md
#endif

subroutine ensure_work_tp(required)
   integer, intent(in) :: required

!$acc wait
   if (allocated(active_domain%work_tp) .and. size(active_domain%work_tp) >= required) return
   if (allocated(active_domain%work_tp)) then
!$acc exit data delete(active_domain%work_tp(1:size(active_domain%work_tp)))
      deallocate(active_domain%work_tp)
   end if
   allocate(active_domain%work_tp(max(1,required)))
!$acc enter data create(active_domain%work_tp(1:max(1,required)))
end subroutine ensure_work_tp

subroutine ensure_recdir_work(required_receive, required_send)
   integer, intent(in) :: required_receive, required_send

!$acc wait
   if (.not.allocated(active_domain%recdir_receive_work) .or. &
      &size(active_domain%recdir_receive_work) < required_receive) then
      if (allocated(active_domain%recdir_receive_work)) then
!$acc exit data delete(active_domain%recdir_receive_work(1:size(active_domain%recdir_receive_work)))
         deallocate(active_domain%recdir_receive_work)
      end if
      allocate(active_domain%recdir_receive_work(max(1,required_receive)))
!$acc enter data create(active_domain%recdir_receive_work(1:max(1,required_receive)))
   end if
   if (.not.allocated(active_domain%recdir_send_work) .or. &
      &size(active_domain%recdir_send_work) < required_send) then
      if (allocated(active_domain%recdir_send_work)) then
!$acc exit data delete(active_domain%recdir_send_work(1:size(active_domain%recdir_send_work)))
         deallocate(active_domain%recdir_send_work)
      end if
      allocate(active_domain%recdir_send_work(max(1,required_send)))
!$acc enter data create(active_domain%recdir_send_work(1:max(1,required_send)))
   end if
   if (.not.allocated(active_domain%recdir_recv_reqs)) then
      allocate(active_domain%recdir_recv_reqs(max(1,active_domain%nproc)))
      active_domain%recdir_recv_reqs = MPI_REQUEST_NULL
   end if
   if (.not.allocated(active_domain%recdir_send_reqs)) then
      allocate(active_domain%recdir_send_reqs(max(1,active_domain%nproc)))
      active_domain%recdir_send_reqs = MPI_REQUEST_NULL
   end if
end subroutine ensure_recdir_work

#if TINKER_MIXED_PREC
subroutine ensure_work_rp(required)
   integer, intent(in) :: required

!$acc wait
   if (allocated(active_domain%work_rp) .and. size(active_domain%work_rp) >= required) return
   if (allocated(active_domain%work_rp)) then
!$acc exit data delete(active_domain%work_rp(1:size(active_domain%work_rp)))
      deallocate(active_domain%work_rp)
   end if
   allocate(active_domain%work_rp(max(1,required)))
!$acc enter data create(active_domain%work_rp(1:max(1,required)))
end subroutine ensure_work_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
subroutine ensure_work_md(required)
   integer, intent(in) :: required

!$acc wait
   if (allocated(active_domain%work_md) .and. size(active_domain%work_md) >= required) return
   if (allocated(active_domain%work_md)) then
!$acc exit data delete(active_domain%work_md(1:size(active_domain%work_md)))
      deallocate(active_domain%work_md)
   end if
   allocate(active_domain%work_md(max(1,required)))
!$acc enter data create(active_domain%work_md(1:max(1,required)))
end subroutine ensure_work_md
#endif

subroutine check_domain(halo)
   integer, intent(in) :: halo

   call check_domain_layout()
   if (halo /= halo_large .and. halo /= halo_neighbour .and. halo /= halo_short) &
      &error stop 'unknown halo neighbourhood'
end subroutine check_domain

module subroutine check_domain_layout()

   if (.not.associated(active_domain)) error stop 'parallel domain is not active'
   if (.not.active_domain%layout_ready) error stop 'parallel domain layout is not configured'
end subroutine check_domain_layout

subroutine check_reduction(visibility, root, root_rank, communicator, process_rank, process_count)
   integer, intent(in) :: visibility
   integer, intent(in), optional :: root
   integer, intent(out) :: root_rank, communicator, process_rank, process_count

   ! Unmigrated executables still initialize only domdec; keep their communicator
   ! behind this adapter while active domains use their owned topology metadata.
   if (associated(active_domain)) then
      communicator = active_domain%communicator
      process_rank = active_domain%rank
      process_count = active_domain%size
   else
      if (get_domain_comm() == MPI_COMM_NULL .or. nproc <= 0) error stop 'parallel domain is not active'
      communicator = get_domain_comm()
      process_rank = rank
      process_count = nproc
   end if
   if (visibility /= reduction_root .and. visibility /= reduction_all) error stop 'unknown reduction visibility'
   root_rank = 0
   if (present(root)) root_rank = root
   if (visibility == reduction_all .and. present(root)) error stop 'all-rank reduction has no root'
   if (visibility == reduction_root .and. (root_rank < 0 .or. root_rank >= process_count)) then
      error stop 'reduction root is outside the domain'
   end if
end subroutine check_reduction

logical function skip_domain_exchange()

   skip_domain_exchange = active_domain%nproc == 1 .or. &
      &(use_pmecore .and. active_domain%rank > ndir-1)
end function skip_domain_exchange

integer function halo_option(halo)
   integer, intent(in) :: halo

   select case (halo)
   case (halo_large)
      halo_option = 1
   case (halo_neighbour)
      halo_option = 2
   case (halo_short)
      halo_option = 3
   case default
      error stop 'unknown halo neighbourhood'
   end select
end function halo_option

integer function send_count(option)
   integer, intent(in) :: option

   select case (option)
   case (1)
      send_count = active_domain%nbig_send
   case (2)
      send_count = active_domain%nneig_send
   case (3)
      send_count = active_domain%nshort_send
   case default
      error stop 'unknown halo option'
   end select
end function send_count

integer function receive_count(option)
   integer, intent(in) :: option

   select case (option)
   case (1)
      receive_count = active_domain%nbig_recep
   case (2)
      receive_count = active_domain%nneig_recep
   case (3)
      receive_count = active_domain%nshort_recep
   case default
      error stop 'unknown halo option'
   end select
end function receive_count

integer function send_peer(option, index)
   integer, intent(in) :: option, index

   if (option == 1) then
      send_peer = active_domain%pbig_send(index)
   else if (option == 2) then
      send_peer = active_domain%pneig_send(index)
   else if (option == 3) then
      send_peer = active_domain%pshort_send(index)
   else
      error stop 'unknown halo option'
   end if
end function send_peer

integer function receive_peer(option, index)
   integer, intent(in) :: option, index

   if (option == 1) then
      receive_peer = active_domain%pbig_recep(index)
   else if (option == 2) then
      receive_peer = active_domain%pneig_recep(index)
   else if (option == 3) then
      receive_peer = active_domain%pshort_recep(index)
   else
      error stop 'unknown halo option'
   end if
end function receive_peer

! commdirdirgpu uses the reciprocal-work channel: its device dipoles are
! produced by the reciprocal induction path and the fixed stream contract is
! intentionally distinct from commfieldgpu's direct-field hand-off.
subroutine commdirdirgpu(nrhs,rule,mu,reqrec1,reqsend1)
   use domdec
   use iounit
   use inform  ,only: deb_Path
   use mdevice  ,only: reciprocal_work
   use mpole
   use mpi
   use sizes   ,only: tinkerdebug
   use timestat,only: timer_enter,timer_exit,timer_polsolvcomm&
      &,quiet_timers
   implicit none
   integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,tag0,i
   integer idomlen,ibufbeg,ipr
   integer, intent(inout) :: reqrec1(nproc), reqsend1(nproc)
   real(t_p) mu(3,nrhs,npolebloc)
   parameter(tag0=0)
1000 format(' illegal rule in commdirdirgpu.')
41 format(7x,'>> ',A20,   3x,'recv')
42 format(7x,   3x,A20,' >>','send')
43 format(7x,'<< ',A20,   3x,'wait')

   if ((n_recep1.eq.0.and.n_send1.eq.0) .or. npolebloc <= 0) return

   call timer_enter( timer_polsolvcomm )
   if (rule.eq.0) then
      if (deb_Path) write(*,41) 'commdirdirgpu       '
      do i = 1, n_recep1
         tag     = nproc*rank + p_recep1(i) + 1
         ibufbeg = bufbegpole(p_recep1(i)+1)
         idomlen = domlenpole(p_recep1(i)+1)
         if (idomlen <= 0) cycle
!$acc host_data use_device(mu)
         call MPI_IRECV(mu(1,1,ibufbeg),3*nrhs*idomlen&
            &,MPI_TPREC,p_recep1(i),tag0&
            &,get_domain_comm(),reqrec1(i),ierr)
!$acc end host_data
      end do
   else if (rule.eq.1) then
!$acc wait(reciprocal_work%queue)
      do i = 1, n_send1
         tag = nproc*p_send1(i) + rank + 1
         if (npoleloc <= 0) cycle
!$acc host_data use_device(mu)
         call MPI_ISEND(mu,3*nrhs*npoleloc&
            &,MPI_TPREC,p_send1(i),tag0&
            &,get_domain_comm(),reqsend1(i),ierr)
!$acc end host_data
      end do
      if (tinkerdebug.gt.0) call MPI_BARRIER(&
         &merge(get_direct_comm(),get_domain_comm(),use_pmecore),ierr)
      if (deb_Path) write(*,42) 'commdirdirgpu       '
   else if (rule.eq.2) then
      do i = 1, n_send1
         call MPI_WAIT(reqsend1(i),status,ierr)
      end do
      do i = 1, n_recep1
         call MPI_WAIT(reqrec1(i),status,ierr)
      end do
      if (tinkerdebug.gt.0) call MPI_BARRIER(&
         &merge(get_direct_comm(),get_domain_comm(),use_pmecore),ierr)
      if (deb_Path) write(*,43) 'commdirdirgpu       '
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if
   call timer_exit( timer_polsolvcomm,quiet_timers )
end

subroutine commdirdir(nrhs,rule,mu,reqrec,reqsend,channel)
   use domdec
   use iounit
   use inform  ,only: deb_Path
   use mpole
   use mpi
   use mdevice, only: channel_t
   use sizes   ,only: tinkerdebug
   use timestat,only:timer_enter,timer_exit,timer_polsolvcomm&
      &,quiet_timers
   implicit none
   integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
   integer :: reqrec(nproc),reqsend(nproc)
   real(t_p) mu(3,nrhs,npolebloc)
   type(channel_t), intent(in) :: channel
1000 format(' illegal rule in commdirdir.')

   call timer_enter( timer_polsolvcomm )
   if (rule.eq.0) then
      do i = 1, n_recep1
         tag = nproc*rank + p_recep1(i) + 1
         call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),&
            &3*nrhs*domlen(p_recep1(i)+1),MPI_TPREC,p_recep1(i),&
            &tag,get_domain_comm(),reqrec(i),ierr)
      end do
   else if (rule.eq.1) then
!$acc wait(channel%queue)
      do i = 1, n_send1
         tag = nproc*p_send1(i) + rank + 1
         call MPI_ISEND(mu,3*nrhs*npoleloc,&
            &MPI_TPREC,p_send1(i),tag,get_domain_comm(),&
            &reqsend(i),ierr)
      end do
   else if (rule.eq.2) then
      do i = 1, n_recep1
         call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_send1
         call MPI_WAIT(reqsend(i),status,ierr)
      end do
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if
   call timer_exit( timer_polsolvcomm,quiet_timers )
end

subroutine commdirdirshort(nrhs,rule,mu,reqrec1,reqsend1,channel)
   use domdec
   use iounit
   use inform  ,only: deb_Path
   use mdevice  ,only: channel_t
   use mpole
   use mpi
   use sizes
   use timestat,only:timer_enter,timer_exit,timer_polsolvcomm&
      &,quiet_timers
   implicit none
   integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
   integer ibufbeg
   integer, intent(inout) :: reqrec1(nproc),reqsend1(nproc)
   real(t_p)  mu(3,nrhs,npolebloc)
   type(channel_t), intent(in) :: channel
1000 format(' illegal rule in commdirdir.')
41 format(7x,'>> ',A20,   3x,'recv')
42 format(7x,   3x,A20,' >>','send')
43 format(7x,'<< ',A20,   3x,'wait')

   if ((n_recepshort1.eq.0.and.n_sendshort1.eq.0)) return

   call timer_enter( timer_polsolvcomm )
   if (rule.eq.0) then
      if (deb_Path) write(*,41) 'commdirdirshort     '
!$acc host_data use_device(mu)
      do i = 1, n_recepshort1
         tag = nproc*rank + p_recepshort1(i) + 1
         ibufbeg = bufbegpole(p_recepshort1(i)+1)
         call MPI_IRECV(mu(1,1,ibufbeg),&
            &3*nrhs*domlenpole(p_recepshort1(i)+1),MPI_TPREC,&
            &p_recepshort1(i),tag,get_domain_comm(),reqrec1(i),ierr)
      end do
!$acc end host_data
   else if (rule.eq.1) then
!$acc wait(channel%queue)
!$acc host_data use_device(mu)
      do i = 1, n_sendshort1
         tag = nproc*p_sendshort1(i) + rank + 1
         call MPI_ISEND(mu,3*nrhs*npoleloc,&
            &MPI_TPREC,p_sendshort1(i),tag,get_domain_comm(),&
            &reqsend1(i),ierr)
      end do
!$acc end host_data
      if (tinkerdebug.gt.0 .and. .not.use_pmecore) &
         &call MPI_BARRIER(get_host_comm(),ierr)
      if (deb_Path) write(*,42) 'commdirdirshort     '
   else if (rule.eq.2) then
      do i = 1, n_recepshort1
         call mpi_wait(reqrec1(i),status,ierr)
      end do
      do i = 1, n_sendshort1
         call mpi_wait(reqsend1(i),status,ierr)
      end do
      if (tinkerdebug.gt.0 .and. .not.use_pmecore) &
         &call MPI_BARRIER(get_host_comm(),ierr)
      if (deb_Path) write(*,43) 'commdirdirshort     '
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if
   call timer_exit( timer_polsolvcomm,quiet_timers )
end

module subroutine commfield(nrhs,ef)
   use domdec
   use mpi
   use mpole
   use potent
   use timestat,only:timer_enter,timer_exit,timer_polfieldcomm&
      &,quiet_timers
   implicit none
   integer, intent(in) :: nrhs
   integer i,j,k,tag,status(MPI_STATUS_SIZE)
   integer ierr
   integer, allocatable :: reqrec(:), reqsend(:)
   real(t_p), intent(inout) :: ef(3,nrhs,*)
   real(t_p), allocatable :: buffer(:,:,:,:)
!
   call timer_enter( timer_polfieldcomm )
!
   allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
   allocate (reqrec(nproc))
   allocate (reqsend(nproc))
!
   do i = 1, n_send1
      tag = nproc*rank + p_send1(i) + 1
      call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,&
         &p_send1(i),tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqrec(i),ierr)
   end do
!
   do i = 1, n_recep1
      tag = nproc*p_recep1(i) + rank + 1
      call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),&
         &3*nrhs*domlenpole(p_recep1(i)+1),MPI_TPREC,p_recep1(i),&
         &tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqsend(i),ierr)
   end do
!
   do i = 1, n_send1
      tag = nproc*rank + p_send1(i) + 1
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, n_recep1
      tag = nproc*p_recep1(i) + rank + 1
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
   do i = 1, n_send1
      do j = 1, npoleloc
         do k = 1, nrhs
            ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
            ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
            ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
         end do
      end do
   end do
!
   deallocate (buffer)
   deallocate (reqrec)
   deallocate (reqsend)
   call timer_exit( timer_polfieldcomm,quiet_timers )
end subroutine commfield

! commfieldgpu uses the direct-work channel: direct-field producers must be
! complete before their device buffer crosses the MPI hand-off.
module subroutine commfieldgpu(nrhs,ef)
   use domdec
   use mdevice,only: direct_work
   use mpi
   use mpole
   use potent
   use inform  ,only: deb_Path
   use sizes   ,only: tinkerdebug
   use timestat,only:timer_enter,timer_exit,timer_polfieldcomm&
      &,quiet_timers
   use tinMemory,only:prmem_request
   implicit none
   integer  ,intent(in)    :: nrhs
   real(t_p),intent(inout) :: ef(3,nrhs,*)
   integer i,i1,j,k,tag,status(MPI_STATUS_SIZE)
   integer ierr
   integer ibufbeg
   integer :: reqrec(nproc), reqsend(nproc)
   real(t_p), allocatable :: buffer(:,:,:,:)
   !real(t_p)  :: buffer(3,nrhs,max(npoleloc,1),n_send1)
!
41 format(7x,'>> ',A20,   3x)
42 format(7x,   3x,A20,' >>')
43 format(7x,'<< ',A20,   3x)

   if ( comm_opt%skip_direct ) return
   if (n_send1.gt.0.or.n_recep1.gt.0) then
      call timer_enter( timer_polfieldcomm )
      if (tinkerdebug.gt.0) call MPI_BARRIER(&
         &merge(get_direct_comm(),get_domain_comm(),use_pmecore),ierr)
      if (deb_Path) write(*,41) 'commfieldgpu        '
!        if (btest(tinkerdebug,0).and.nproc.eq.4) then
!22         format('_domlenpole ',I5,I7,3x,5I7,' _')
!           write(0,22) rank,npoleloc,domlenpole(1:4)
!        end if


      call prmem_request(buffer,3,nrhs,max(npoleloc,1),n_send1,async=.false.)
!
      do i = 1, n_send1
         tag = nproc*rank + p_send1(i) + 1
!$acc host_data use_device(buffer)
         call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,&
            &p_send1(i),tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqrec(i),ierr)
!$acc end host_data
      end do
!
!$acc wait(direct_work%queue)
      do i = 1, n_recep1
         tag = nproc*p_recep1(i) + rank + 1
         ibufbeg = bufbegpole(p_recep1(i)+1)
!$acc host_data use_device(ef)
         call MPI_ISEND(ef(1,1,ibufbeg),&
            &3*nrhs*domlenpole(p_recep1(i)+1),MPI_TPREC,&
            &p_recep1(i),tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqsend(i),ierr)
!$acc end host_data
      end do
!
      do i = 1, n_send1
         tag = nproc*rank + p_send1(i) + 1
         call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_recep1
         tag = nproc*p_recep1(i) + rank + 1
         call MPI_WAIT(reqsend(i),status,ierr)
      end do
!
!$acc parallel loop gang vector collapse(3) &
!$acc         present(ef,buffer) async(direct_work%queue)
      do j = 1, npoleloc
         do k = 1, nrhs
            do i1=1,3
               do i = 1, n_send1
                  ef(i1,k,j) = ef(i1,k,j) + buffer(i1,k,j,i)
               end do
            end do
         end do
      end do

!$acc wait(direct_work%queue)
      call prmem_request(buffer,0,0,0,0,async=.false.)
      if (tinkerdebug.gt.0) call MPI_BARRIER(&
         &merge(get_direct_comm(),get_domain_comm(),use_pmecore),ierr)
      if (deb_Path) write(*,43) 'commfieldgpu        '
      call timer_exit( timer_polfieldcomm,quiet_timers )
   end if
end subroutine commfieldgpu

module subroutine commfieldshort(nrhs,ef,channel)
   use domdec
   use mdevice,only: channel_t
   use mpi
   use mpole
   use potent
   use inform  ,only: deb_Path
   use sizes   ,only: tinkerdebug
   use timestat,only:timer_enter,timer_exit,timer_polfieldcomm&
      &,quiet_timers
   use tinMemory ,only: prmem_request
   implicit none
   integer, intent(in) :: nrhs
   integer i,j,k,l,tag,status(mpi_status_size)
   integer ierr,ibufbeg
   integer, allocatable :: reqrec(:), reqsend(:)
   real(t_p), intent(inout) :: ef(3,nrhs,*)
   type(channel_t), intent(in) :: channel
   real(t_p), allocatable :: buffer(:,:,:,:)
!
   if (comm_opt%skip_direct) return
   if (n_sendshort1.eq.0.and.n_recepshort1.eq.0) return
41 format(7x,'>> ',A20,   3x)
42 format(7x,   3x,A20,' >>')
43 format(7x,'<< ',A20,   3x)

   call timer_enter( timer_polfieldcomm )
   if (tinkerdebug.gt.0) call MPI_BARRIER(get_host_comm(),ierr)
   if (deb_Path) write(*,41) 'commfieldshort      '
!        if (btest(tinkerdebug,0).and.nproc.eq.4) then
!22         format('_',A6,I5,3x,5I7,'_ ')
!           write(0,22) 'npoleloc',rank,npoleloc
!        end if

!
   call prmem_request(buffer,3,nrhs,max(npoleloc,1),n_sendshort1,&
      &async=.true.)
   allocate (reqrec(nproc))
   allocate (reqsend(nproc))
!
!$acc host_data use_device(buffer)
   do i = 1, n_sendshort1
      tag = nproc*rank + p_sendshort1(i) + 1
      call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,&
         &p_sendshort1(i),tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqrec(i),ierr)
   end do
!$acc end host_data
!
!$acc wait(channel%queue)
!$acc host_data use_device(ef)
   do i = 1, n_recepshort1
      tag     = nproc*p_recepshort1(i) + rank + 1
      ibufbeg = bufbegpole(p_recepshort1(i)+1)
      call MPI_ISEND(ef(1,1,ibufbeg),&
         &3*nrhs*domlenpole(p_recepshort1(i)+1),MPI_TPREC,&
         &p_recepshort1(i),tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqsend(i),ierr)
   end do
!$acc end host_data
!
   do i = 1, n_sendshort1
      tag = nproc*rank + p_sendshort1(i) + 1
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, n_recepshort1
      tag = nproc*p_recepshort1(i) + rank + 1
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!$acc parallel loop gang vector collapse(3) &
!$acc         default(present) async(channel%queue)
   do j = 1, npoleloc
      do k = 1, nrhs
         do l = 1,3
            do i = 1, n_sendshort1
               ef(l,k,j) = ef(l,k,j) + buffer(l,k,j,i)
            end do
         end do
      end do
   end do
!
   deallocate (reqrec)
   deallocate (reqsend)
!$acc wait(channel%queue)
   call prmem_request(buffer,0,0,0,0,async=.false.)

   if (tinkerdebug.gt.0) call MPI_BARRIER(get_host_comm(),ierr)
   if (deb_Path) write(*,43) 'commfieldshort      '
   call timer_exit( timer_polfieldcomm,quiet_timers )
end subroutine commfieldshort

module subroutine commpot(pot,rule)
   use atmlst
   use atoms
   use charge
   use domdec
   use iounit
   use mpole
   use mpi
   implicit none
   integer i,tag,ierr,iproc
   integer iloc
   integer, intent(in) :: rule
   integer, allocatable :: reqrec(:),reqsend(:)
   integer :: n_recep, n_send
   integer, allocatable :: p_send(:),p_recep(:)
   real*8, allocatable :: buffer(:)
   real*8, allocatable :: buffers(:)
   real*8, intent(inout) :: pot(*)
   integer status(MPI_STATUS_SIZE)
1000 format(' illegal rule in commpot.')
!
   allocate (p_send(nproc))
   allocate (p_recep(nproc))
!
   if (rule.eq.0) then
      n_recep = nbig_recep
      n_send  = nbig_send
      p_send  = pbig_send
      p_recep = pbig_recep
   else if (rule.eq.1) then
      n_recep = nneig_recep
      n_send  = nneig_send
      p_send  = pneig_send
      p_recep = pneig_recep
   else
      if (rank.eq.0) write(iout,1000)
   end if

!
!     allocate some arrays
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
   allocate (buffer(max(1,nbloc)))
   allocate (buffers(max(1,nloc)))
!
!     MPI : begin reception in buffer
!
   do i = 1, n_recep
      tag = nproc*rank + p_recep(i) + 1
      call MPI_IRECV(buffer(bufbeg(p_recep(i)+1)),&
         &domlen(p_recep(i)+1),&
         &MPI_REAL8,p_recep(i),tag,get_domain_comm(),reqrec(i),ierr)
   end do
!
!     MPI : move in buffer
!
   do i = 1, nloc
      buffers(i) = pot(i)
   end do
!
!     MPI : begin sending
!
   do i = 1, n_send
      tag = nproc*p_send(i) + rank + 1
      call MPI_ISEND(buffers,nloc,&
         &MPI_REAL8,p_send(i),tag,get_domain_comm(),&
         &reqsend(i),ierr)
   end do
!
   do i = 1, n_recep
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, n_send
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!     MPI : move in global arrays
!
   do iproc = 1, n_recep
      do i = 1, domlen(p_recep(iproc)+1)
         iloc = bufbeg(p_recep(iproc)+1)+i-1
         pot(iloc) = buffer(iloc)
      end do
   end do
!
   deallocate (p_send)
   deallocate (p_recep)
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (buffer)
   deallocate (buffers)
   return
end subroutine commpot

module subroutine commpotsum(pot,rule)
   use atmlst
   use atoms
   use charge
   use domdec
   use iounit
   use mpole
   use mpi
   implicit none
   integer i,tag,ierr
   integer, intent(in) :: rule
   integer, allocatable :: reqrec(:),reqsend(:)
   integer :: n_recep, n_send
   integer, allocatable :: p_send(:),p_recep(:)
   real*8, allocatable :: buffer(:,:)
   real*8, allocatable :: buffers(:)
   real*8, intent(inout) :: pot(*)
   integer status(MPI_STATUS_SIZE)
1000 format(' illegal rule in commpot.')
!
   allocate (p_send(nproc))
   allocate (p_recep(nproc))
!
   if (rule.eq.0) then
      n_recep = nbig_recep
      n_send  = nbig_send
      p_send  = pbig_send
      p_recep = pbig_recep
   else if (rule.eq.1) then
      n_recep = nneig_recep
      n_send  = nneig_send
      p_send  = pneig_send
      p_recep = pneig_recep
   else
      if (rank.eq.0) write(iout,1000)
   end if

!
!     allocate some arrays
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
   allocate (buffer(max(1,nloc),n_send))
   allocate (buffers(max(1,nbloc)))
!
!     MPI : move in buffer
!
   buffers((nloc+1):nbloc) = pot((nloc+1):nbloc)
!
   do i = 1, n_send
      tag = nproc*rank + p_send(i) + 1
      call MPI_IRECV(buffer(1,i),nloc,MPI_REAL8,p_send(i),&
         &tag,get_domain_comm(),reqrec(i),ierr)
   end do
!
!     MPI : begin sending
!
   do i = 1, n_recep
      tag = nproc*p_recep(i) + rank + 1
      call MPI_ISEND(buffers(bufbeg(p_recep(i)+1)),&
         &domlen(p_recep(i)+1),&
         &MPI_REAL8,p_recep(i),tag,get_domain_comm(),&
         &reqsend(i),ierr)
   end do
!
   do i = 1, n_send
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, n_recep
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!     MPI : move in global arrays
!
   do i = 1, n_send
      pot(1:nloc) = pot(1:nloc) + buffer(1:nloc,i)
   end do
!
!
   deallocate (p_send)
   deallocate (p_recep)
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (buffer)
   deallocate (buffers)
   return
end subroutine commpotsum

module subroutine commpotrec(pot,potrec)
   use sizes
   use deriv
   use domdec
   use potent
   use mpi
   implicit none
   integer i,j,tag,ierr
   integer jloc,jglob,jlocrec
   integer, allocatable :: reqsend(:),reqrec(:)
   real*8, allocatable :: buffer(:,:)
   real*8, allocatable :: buffers(:)
   real*8, intent(inout) :: pot(*),potrec(*)
   integer status(MPI_STATUS_SIZE)
!
!     allocate some arrays
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
   allocate (buffer(max(1,nloc),nrecdir_send1))
   buffer = 0d0
   allocate (buffers(max(1,nblocrecdir)))
   buffers = 0d0
!
!     do rec dir communications
!
!     MPI : begin reception in buffer
!
   do i = 1, nrecdir_send1
      tag = nproc*rank + precdir_send1(i) + 1
      call MPI_IRECV(buffer(1,i),nloc,&
         &MPI_REAL8,precdir_send1(i),tag,get_domain_comm(),reqrec(i),ierr)
   end do
!
!     Move in array
!
   do i = 1, nrecdir_recep1
      do j = 1, domlen(precdir_recep1(i)+1)
         jloc = bufbeg(precdir_recep1(i)+1)+j-1
         jglob = glob(jloc)
         if (repartrec(jglob).eq.rank) then
            jlocrec = locrec(jglob)
            buffers(jloc) = potrec(jlocrec)
         end if
      end do
   end do
!
   do i = 1, nrecdir_recep1
      tag = nproc*precdir_recep1(i) + rank + 1
      call MPI_ISEND(buffers(bufbeg(precdir_recep1(i)+1)),&
         &domlen(precdir_recep1(i)+1),&
         &MPI_REAL8,precdir_recep1(i),tag,get_domain_comm(),reqsend(i),ierr)
   end do
!
   do i = 1, nrecdir_send1
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, nrecdir_recep1
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!     MPI : move in global arrays
!
   do i = 1, nrecdir_send1
      do j = 1, nloc
         pot(j) = pot(j) + buffer(j,i)
      end do
   end do

   return
end subroutine commpotrec

module subroutine commchgflx(pdelta)
   use atmlst
   use atoms
   use charge
   use domdec
   use mpole
   use mpi
   use potent
   implicit none
   integer i,tag,ierr,iproc
   integer iloc,iglob,j
   integer, allocatable :: reqrec(:),reqsend(:)
   real(t_p), allocatable :: buffer(:,:),buffer2(:)
   real(t_p), allocatable :: buffers(:)
   real(t_p), intent(inout) :: pdelta(*)
   integer status(MPI_STATUS_SIZE)
!
!     allocate some arrays
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
   allocate (buffer(max(1,nloc),nbig_send))
   allocate (buffers(max(1,nbloc)))
!
   do i = 1, nbig_recep
      do j = 1, domlen(pbig_recep(i)+1)
         iloc = bufbeg(pbig_recep(i)+1)+j-1
         iglob = glob(iloc)
         buffers(iloc) = pdelta(iglob)
      end do
   end do
!
!     MPI : begin reception in buffer
!
   do i = 1, nbig_send
      tag = nproc*rank + pbig_send(i) + 1
      call MPI_IRECV(buffer(1,i),nloc,MPI_TPREC,pbig_send(i),&
         &tag,get_domain_comm(),reqrec(i),ierr)
   end do
!
!     MPI : begin sending
!
   do i = 1, nbig_recep
      tag = nproc*pbig_recep(i) + rank + 1
      call MPI_ISEND(buffers(bufbeg(pbig_recep(i)+1)),&
         &domlen(pbig_recep(i)+1),&
         &MPI_TPREC,pbig_recep(i),tag,get_domain_comm(),&
         &reqsend(i),ierr)
   end do
!
   do i = 1, nbig_send
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, nbig_recep
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!     MPI : move in global arrays
!
   do iproc = 1, nbig_send
      do i = 1, nloc
         iglob = glob(i)
         pdelta(iglob) = pdelta(iglob) + buffer(i,iproc)
      end do
   end do
!
!     then send accumulated values to neighbors
!
   deallocate (buffer)
   allocate (buffer2(max(1,nbloc)))
!
!     MPI : begin reception in buffer
!
   do i = 1, nbig_recep
      tag = nproc*rank + pbig_recep(i) + 1
      call MPI_IRECV(buffer2(bufbeg(pbig_recep(i)+1)),&
         &domlen(pbig_recep(i)+1),&
         &MPI_TPREC,pbig_recep(i),tag,get_domain_comm(),reqrec(i),ierr)
   end do
!
!     MPI : move in buffer
!
   do i = 1, nloc
      iglob = glob(i)
      buffers(i) = pdelta(iglob)
   end do
!
!     MPI : begin sending
!
   do i = 1, nbig_send
      tag = nproc*pbig_send(i) + rank + 1
      call MPI_ISEND(buffers,nloc,&
         &MPI_TPREC,pbig_send(i),tag,get_domain_comm(),&
         &reqsend(i),ierr)
   end do
!
   do i = 1, nbig_recep
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, nbig_send
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!     MPI : move in global arrays
!
   do iproc = 1, nbig_recep
      do i = 1, domlen(pbig_recep(iproc)+1)
         iloc = bufbeg(pbig_recep(iproc)+1)+i-1
         iglob = glob(iloc)
         pdelta(iglob) = buffer2(iloc)
      end do
   end do
!
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (buffer2)
   deallocate (buffers)
   return
end subroutine commchgflx

subroutine commrecdirdip(nrhs,rule,diprec,dip,buffermpimu1,&
   &buffermpimu2,req2rec,req2send)
   use domdec
   use iounit
   use mpole
   use mpi
   use timestat,only:timer_enter,timer_exit,timer_polsolvcomm&
      &,quiet_timers
   implicit none
   integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
   integer :: req2rec(nproc),req2send(nproc)
   real(t_p) dip(3,nrhs,max(1,npoleloc))
   real(t_p) diprec(3,nrhs,max(1,npolerecloc))
   real(t_p) buffermpimu1(3,nrhs,*)
   real(t_p) buffermpimu2(3,nrhs,*)
1000 format(' illegal rule in commrecdirdip.')

   call timer_enter( timer_polsolvcomm )
   if (rule.eq.0) then
      do i = 1, nrecdir_recep1
         proc = precdir_recep1(i)
         if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),3*nrhs&
               &*buflen2(proc+1),MPI_TPREC,proc,tag,get_domain_comm(),&
               &req2rec(i),ierr)
         end if
      end do
   else if (rule.eq.1) then
      do i = 1, nrecdir_send1
         proc = precdir_send1(i)
         if (proc.ne.rank) then
            do j = 0, buflen1(proc+1)-1
               call amove(3*nrhs,dip(1,1,poleloc(&
                  &buf1(bufbeg1(proc+1)+j))),&
                  &buffermpimu1(1,1,bufbeg1(proc+1)+j))
            end do
         end if
      end do
      do i = 1, nrecdir_send1
         proc = precdir_send1(i)
         if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),3*nrhs&
               &*buflen1(proc+1),MPI_TPREC,proc,tag,get_domain_comm(),&
               &req2send(i),ierr)
         end if
      end do
   else if (rule.eq.2) then
      do i = 1, nrecdir_send1
         proc = precdir_send1(i)
         if (proc.ne.rank) call MPI_WAIT(req2send(i),status,ierr)
      end do
      do i = 1, nrecdir_recep1
         proc = precdir_recep1(i)
         if (proc.ne.rank) then
            call MPI_WAIT(req2rec(i),status,ierr)
            do j = 0, buflen2(proc+1)-1
               call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),&
                  &diprec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
            end do
         end if
      end do
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if
   call timer_exit( timer_polsolvcomm,quiet_timers )
end

subroutine commrecdirdipgpu(nrhs,rule,diprec,dip,buffermpimu1,&
   &buffermpimu2,req2rec,req2send)
   use domdec
   use iounit
   use inform  ,only:deb_path
   use mpole
   use mpi
   use timestat,only:timer_enter,timer_exit,timer_polsolvcomm&
      &,quiet_timers
   use sizes   ,only:tinkerdebug
   implicit none
   integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,k,l,proc
   integer lenbuf,ibufbeg
   integer :: req2rec(nproc),req2send(nproc)
   real(t_p) dip(3,nrhs,max(1,npoleloc))
   real(t_p) diprec(3,nrhs,max(1,npolerecloc))
   real(t_p) buffermpimu1(3,nrhs,*)
   real(t_p) buffermpimu2(3,nrhs,*)
1000 format(' illegal rule in commrecdirdipgpu.')
41 format(7x,'>> ',A20,   3x,'recv')
42 format(7x,   3x,A20,' >>','send')
43 format(7x,'<< ',A20,   3x,'wait')

   if (nrecdir_send1.eq.nrecdir_recep1.and.nrecdir_recep1.le.1)&
      &return
   call timer_enter( timer_polsolvcomm )
   if (rule.eq.0) then
      do i = 1, nrecdir_recep1
         proc = precdir_recep1(i)
         if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
!$acc host_data use_device(buffermpimu2)
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),&
               &3*nrhs*buflen2(proc+1),MPI_TPREC,proc,tag,&
               &get_domain_comm(),req2rec(i),ierr)
!$acc end host_data
         end if
      end do
   else if (rule.eq.1) then
      do i = 1, nrecdir_send1
         proc = precdir_send1(i)
         if (proc.ne.rank) then
            lenbuf = buflen1(proc+1)-1
            ibufbeg = bufbeg1(proc+1)
!$acc parallel loop gang collapse(3) async &
!$acc         default(present)
            do j = 0,lenbuf
               do k = 1, 2
                  do l = 1, 3
                     buffermpimu1(l,k,ibufbeg+j)&
                        &= dip(l,k,poleloc(buf1(ibufbeg+j)))
                  end do
               end do
            end do
         end if
      end do
!$acc wait
      do i = 1, nrecdir_send1
         proc = precdir_send1(i)
         if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
!$acc host_data use_device(buffermpimu1)
            call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),&
               &3*nrhs*buflen1(proc+1),MPI_TPREC,proc,tag,&
               &get_domain_comm(),req2send(i),ierr)
!$acc end host_data
         end if
      end do
      if (tinkerdebug.gt.0 .and. .not.use_pmecore) &
         &call MPI_BARRIER(get_host_comm(),ierr)
      if (deb_Path) write(*,42) 'commrecdirdipgpu    '
   else if (rule.eq.2) then
      do i = 1, nrecdir_send1
         proc = precdir_send1(i)
         if (proc.ne.rank) call MPI_WAIT(req2send(i),status,ierr)
      end do
      do i = 1, nrecdir_recep1
         proc = precdir_recep1(i)
         if (proc.ne.rank) then
            call MPI_WAIT(req2rec(i),status,ierr)
            lenbuf = buflen2(proc+1)-1
            ibufbeg = bufbeg2(proc+1)
!$acc parallel loop gang collapse(3) async &
!$acc         default(present)
            do j = 0, lenbuf
               do k = 1, 2
                  do l = 1, 3
                     diprec(l,k,polerecloc(buf2(ibufbeg+j)))&
                        &= buffermpimu2(l,k,ibufbeg+j)
                  end do
               end do
            end do
         end if
      end do
!$acc wait
      if (tinkerdebug.gt.0 .and. .not.use_pmecore) &
         &call MPI_BARRIER(get_host_comm(),ierr)
      if (deb_Path) write(*,43) 'commrecdirdipgpu    '
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if
   call timer_exit( timer_polsolvcomm,quiet_timers )
end

!
!     subroutine commposrec: communicate positions for reciprocal space
!
module subroutine commposrec
   use atomsMirror
   use domdec
   use inform   ,only: deb_Path
   use mpi
   use potent   ,only: use_pmecore,use_mlpot_only
   use pme      ,only: GridDecomp1d
   use timestat ,only: timer_enter,timer_exit,timer_commstep&
      &,quiet_timers
   use tinMemory,only: prmem_requestm
   implicit none
   integer i,iproc,iprec
   integer iglob,iloc
   integer tag,ierr,ibufbeg,idomlen
   integer rankloc,nprocloc
   integer status(MPI_STATUS_SIZE)
   real(r_p),pointer :: buffer(:,:),buffers(:,:),pool1(:),pool2(:)
   integer, allocatable :: reqsend(:),reqrec(:)
!
   if (nproc.eq.1.or.use_mlpot_only) return
   if (GridDecomp1d.and.Bdecomp1d.and..not.use_pmecore) return

   allocate(reqsend(max(1,nproc)),reqrec(max(1,nproc)))
   reqsend = MPI_REQUEST_NULL
   reqrec = MPI_REQUEST_NULL

   call timer_enter( timer_commstep )
   if (deb_Path) write(*,*) '   >> commposrec'
!
   if (use_pmecore) then
      nprocloc = nrec
      rankloc  = rank_bis
   else
      nprocloc = nproc
      rankloc  = rank
   end if
!
!     allocate some arrays
!
   call domain_workspace_rp(1,3*max(nblocrecdir,nlocrec2),pool1)
   call domain_workspace_rp(2,3*max(nloc,nlocrec),pool2)
   buffer (1:3,1:nblocrecdir)=>pool1(1:3*nblocrecdir)
   buffers(1:3,1:nloc)       =>pool2(1:3*nloc)
!
!     first do dir rec communications
!
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(glob,buffer,buffers,x,y,z)
   do i = 1, nrecdir_recep2
      if (precdir_recep2(i).ne.rank) then
         tag = nproc*rank + precdir_recep2(i) + 1
         call MPI_IRECV(buffer(1,bufbeg(precdir_recep2(i)+1)),&
            &3*domlen(precdir_recep2(i)+1),MPI_RPREC,&
            &precdir_recep2(i),tag,get_domain_comm(),reqrec(i),ierr)
      end if
   end do
!
!     MPI : move in buffer
!
!$acc parallel loop async deviceptr(glob,x,y,z,buffers)
   do i = 1, nloc
      iglob        = glob(i)
      buffers(1,i) = x(iglob)
      buffers(2,i) = y(iglob)
      buffers(3,i) = z(iglob)
   end do
!
!     MPI : begin sending
!
!$acc wait
   do i = 1, nrecdir_send2
      if (precdir_send2(i).ne.rank) then
         tag = nproc*precdir_send2(i) + rank + 1
         call MPI_ISEND(buffers,3*nloc,MPI_RPREC,precdir_send2(i),&
            &tag,get_domain_comm(),reqsend(i),ierr)
      end if
   end do

   do i = 1, nrecdir_recep2
      if (precdir_recep2(i).ne.rank) then
         call MPI_Wait(reqrec(i),status,ierr)
      end if
   end do
!
!     MPI : move in global arrays
!
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         call MPI_WAIT(reqrec(iproc),status,ierr)
         ibufbeg = bufbeg(precdir_recep2(iproc)+1)
         idomlen = domlen(precdir_recep2(iproc)+1)
!$acc parallel loop async &
!$acc         deviceptr(glob,x,y,z,buffer)
         do i = 1, idomlen
            iloc     = ibufbeg+i-1
            iglob    = glob(iloc)
            x(iglob) = buffer(1,iloc)
            y(iglob) = buffer(2,iloc)
            z(iglob) = buffer(3,iloc)
         end do
      end if
   end do
!$acc end host_data
!
   do i = 1, nrecdir_send2
      if (precdir_send2(i).ne.rank) then
         call MPI_WAIT(reqsend(i),status,ierr)
      end if
   end do
!
   nullify(buffer)
   nullify(buffers)
!$acc wait
!
!     then do rec rec communications
!
   buffer (1:3,1:nlocrec2) =>pool1(1:3*nlocrec2)
   buffers(1:3,1:nlocrec)  =>pool2(1:3*nlocrec)
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(globrec,buffer,buffers,x,y,z)
   do i = 1, nrec_recep
      tag = nprocloc*rankloc + prec_recep(i) + 1
      call MPI_IRECV(buffer(1,bufbegrec(prec_recep(i)+1))&
         &,3*domlenrec(prec_recep(i)+1),MPI_RPREC&
         &,prec_recep(i),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(i),ierr)
   end do

!     MPI : move in buffer
!$acc parallel loop deviceptr(globrec,buffers,x,y,z)
   do i = 1, nlocrec
      iglob        = globrec(i)
      buffers(1,i) = x(iglob)
      buffers(2,i) = y(iglob)
      buffers(3,i) = z(iglob)
   end do
!
!     MPI : begin sending
!
   do i = 1, nrec_send
      tag = nprocloc*prec_send(i) + rankloc+1
      call MPI_ISEND(buffers,3*nlocrec,MPI_RPREC,prec_send(i)&
         &,tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(i),ierr)
   end do

   ! Wait for receiving requests
   do i = 1, nrec_recep
      call MPI_WAIT(reqrec(i),status,ierr)
   end do

!     !MPI : move in global arrays
   do iproc = 1, nrec_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
      ibufbeg = bufbegrec(prec_recep(iproc)+1)-1
      idomlen = domlenrec(prec_recep(iproc)+1)
!$acc parallel loop async deviceptr(globrec,buffer,x,y,z)
      do i = 1, idomlen
         iloc     = ibufbeg+i
         iglob    = globrec(iloc)
         x(iglob) = buffer(1,iloc)
         y(iglob) = buffer(2,iloc)
         z(iglob) = buffer(3,iloc)
      end do
   end do
!$acc end host_data
   call reCast_position

   ! Wait for sending requests
   do i = 1, nrec_send
      call MPI_WAIT(reqsend(i),status,ierr)
   end do

   nullify(buffer)
   nullify(buffers)
   deallocate(reqsend,reqrec)

   if (deb_Path) write(*,*) '   << commposrec'
   call timer_exit( timer_commstep,quiet_timers )
end subroutine commposrec
!
!     subroutine commpos : communicate positions to neighboring processes
!
module subroutine commpos
   use atomsMirror
   use domdec
   use inform   ,only: deb_Path
   use mdevice  ,only: direct_work
#ifdef _OPENACC
   use mdevice  ,only: after,reciprocal_work
   use status   ,only: failed,status_t
#endif
   use potent   ,only: use_pmecore
   use timestat ,only: timer_enter,timer_exit,timer_commstep&
      &,quiet_timers
   implicit none
   integer i,iproc,iloc,idomlen,iglob
   real(r_p), pointer :: pool1(:)
#ifdef _OPENACC
   type(status_t) :: dep_status
#endif
!
!
!     communicate positions
!
   if (nproc.eq.1.or.(use_pmecore).and.(rank.gt.ndir-1)) return
   call timer_enter( timer_commstep )
   if (deb_Path) write(*,*) '   >> commpos'

   call domain_workspace_rp(1,3*nbloc,pool1)
!$acc parallel loop async(direct_work%queue) default(present)
   do i = 1, nloc
      iglob = glob(i)
      pool1(3*i-2) = x(iglob)
      pool1(3*i-1) = y(iglob)
      pool1(3*i) = z(iglob)
   end do
   call extend(pool1,3,halo_large,channel=direct_work)
   do iproc = 1, nbig_recep
      iloc = bufbeg(pbig_recep(iproc)+1)
      idomlen = domlen(pbig_recep(iproc)+1)
!$acc parallel loop async(direct_work%queue)
      do i = 1, idomlen
         iglob  = glob(iloc+i-1)
         x(iglob) = pool1(3*(iloc+i-1)-2)
         y(iglob) = pool1(3*(iloc+i-1)-1)
         z(iglob) = pool1(3*(iloc+i-1))
      end do
   end do
#ifdef _OPENACC
   call after(reciprocal_work,direct_work,dep_status)
   if (failed(dep_status)) __TINKER_FATAL__
   call reCast_position(reciprocal_work%queue)
#else
   call reCast_position
#endif

   if (deb_Path) write(*,*) '   << commpos'
   call timer_exit( timer_commstep,quiet_timers )
end subroutine commpos
!


!     subroutine commposrespa : deal with communications of positions between two time steps
!     with respa integrator
!
module subroutine commposrespa(fast)
   use atomsMirror
   use domdec
   use inform   ,only: deb_Path
   use mdevice  ,only: direct_work
#ifdef _OPENACC
   use mdevice  ,only: after,reciprocal_work
   use status   ,only: failed,status_t
#endif
   use potent   ,only: use_pmecore
   use timestat ,only: timer_enter,timer_exit,timer_commstep&
      &,quiet_timers
   use sizes    ,only: tinkerdebug
   implicit none
   integer i,iproc,iloc,idomlen,iglob
   real(r_p), pointer :: pool1(:)
#ifdef _OPENACC
   type(status_t) :: dep_status
#endif
   logical, intent(in) :: fast
!
   if (ndir.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

   if (deb_Path) then
      write(*,*) '   >> commposrespa',fast
   end if
   call timer_enter( timer_commstep )
!
   call domain_workspace_rp(1,3*nbloc,pool1)
!$acc parallel loop async(direct_work%queue) default(present)
   do i = 1, nloc
      iglob = glob(i)
      pool1(3*i-2) = x(iglob)
      pool1(3*i-1) = y(iglob)
      pool1(3*i) = z(iglob)
   end do
   call extend(pool1,3,merge(halo_neighbour,halo_large,fast),channel=direct_work)
   if (fast) then
      do iproc = 1, nneig_recep
         iloc = bufbeg(pneig_recep(iproc)+1)
         idomlen = domlen(pneig_recep(iproc)+1)
!$acc parallel loop async(direct_work%queue)
         do i = 1, idomlen
            iglob = glob(iloc+i-1)
            x(iglob) = pool1(3*(iloc+i-1)-2)
            y(iglob) = pool1(3*(iloc+i-1)-1)
            z(iglob) = pool1(3*(iloc+i-1))
         end do
      end do
   else
      do iproc = 1, nbig_recep
         iloc = bufbeg(pbig_recep(iproc)+1)
         idomlen = domlen(pbig_recep(iproc)+1)
!$acc parallel loop async(direct_work%queue)
         do i = 1, idomlen
            iglob = glob(iloc+i-1)
            x(iglob) = pool1(3*(iloc+i-1)-2)
            y(iglob) = pool1(3*(iloc+i-1)-1)
            z(iglob) = pool1(3*(iloc+i-1))
         end do
      end do
   end if
#ifdef _OPENACC
   call after(reciprocal_work,direct_work,dep_status)
   if (failed(dep_status)) __TINKER_FATAL__
   call reCast_position(reciprocal_work%queue)
#else
   call reCast_position
#endif

   if (deb_Path) then
      write(*,*) '   << commposrespa'
   end if
   call timer_exit( timer_commstep,quiet_timers )
end subroutine commposrespa


!

!     subroutine orderbuffer : get the indexes of the particules in the neighboring processes
!
module subroutine orderbuffer(init)
   use sizes
   use atoms  ,only: n
   use inform ,only: deb_Path
   use cutoff
   use domdec
   use mpi
   implicit none
   integer, allocatable :: counttemp(:),ind1temp(:)
   integer count
   integer i,iproc,iproc1,tag,ierr

   integer, allocatable :: reqsend(:),reqrec(:)
   integer status(MPI_STATUS_SIZE)
   logical, intent(in) :: init
!
   allocate (counttemp(nproc))
   allocate (ind1temp(n))
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
!
   if (tinkerdebug.gt.0) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   if (deb_Path) write(*,*) '   >> orderbuffer',init
!
!      ind1temp = 0
   counttemp = 0
!
!     get the domlen values of the coupled real space processes
!
   do iproc = 1, nbig_recep
      tag = nproc*rank + pbig_recep(iproc) + 1
      call MPI_IRECV(domlen(pbig_recep(iproc)+1),1,MPI_INT,&
         &pbig_recep(iproc),tag,get_domain_comm(),reqrec(iproc),ierr)
   end do
   do iproc = 1, nbig_send
      tag = nproc*pbig_send(iproc) + rank + 1
      call MPI_ISEND(domlen(rank+1),1,MPI_INT,pbig_send(iproc),tag,&
         &get_domain_comm(),reqsend(iproc),ierr)
   end do
   do iproc = 1, nbig_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
   do iproc = 1, nbig_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!
!     get the domlen values of the coupled real space processes (dir-recip neighbors but not present in the previous list)
!
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         tag = nproc*rank + precdir_recep2(iproc) + 1
         call MPI_IRECV(domlen(precdir_recep2(iproc)+1),1,MPI_INT,&
            &precdir_recep2(iproc),tag,get_domain_comm(),reqrec(iproc),ierr)
      end if
   end do
   do iproc = 1, nrecdir_send2
      if (precdir_send2(iproc).ne.rank) then
         tag = nproc*precdir_send2(iproc) + rank + 1
         call MPI_ISEND(domlen(rank+1),1,MPI_INT,precdir_send2(iproc),&
            &tag,get_domain_comm(),reqsend(iproc),ierr)
      end if
   end do
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         call MPI_WAIT(reqrec(iproc),status,ierr)
      end if
   end do
   do iproc = 1, nrecdir_send2
      if (precdir_send2(iproc).ne.rank) then
         call MPI_WAIT(reqsend(iproc),status,ierr)
      end if
   end do
!
   if (init) then
      do i = 1, n
!
!       check if atom is in the domain or in one of the neighboring domains
!
         counttemp(repart(i)+1) = counttemp(repart(i)+1) + 1
         ind1temp(i) = counttemp(repart(i)+1)
         if (repart(i).eq.rank) then
            glob(ind1temp(i)) = i
            loc(i) = ind1temp(i)
         end if
      end do
   end if
!
   bufbeg(rank+1) = 1
   count = domlen(rank+1)
   do iproc = 1, nbig_recep
      if (domlen(pbig_recep(iproc)+1).ne.0) then
         bufbeg(pbig_recep(iproc)+1) = count + 1
      else
         bufbeg(pbig_recep(iproc)+1) = 1
      end if
      count = count + domlen(pbig_recep(iproc)+1)
   end do
   nbloc = count
!     nvdwblocloop is nvdwbloc if nvdwbloc is a multiple of 16, or the first one greater
   nblocloop = merge( nbloc,&
      &(int(nbloc/16)+1)*16,&
      &(mod(nbloc,16).eq.0))

   nblocrecdir = count
!!$acc update device(nbloc,nblocloop) async
   do iproc = 1, nrecdir_recep1
      if (precdir_recep1(iproc).eq.rank) goto 10
      do iproc1 = 1, nbig_recep
         if (pbig_recep(iproc1).eq.precdir_recep1(iproc)) goto 10
      end do
      if (domlen(precdir_recep1(iproc)+1).ne.0) then
         bufbeg(precdir_recep1(iproc)+1) = count + 1
      else
         bufbeg(precdir_recep1(iproc)+1) = 1
      end if
      count = count + domlen(precdir_recep1(iproc)+1)
10    continue
   end do
   nblocrecdir = count
!
!     send and receive the indexes of the neighboring processes
!
   do iproc = 1, nbig_recep
      tag = nproc*rank + pbig_recep(iproc) + 1
      call MPI_IRECV(glob(bufbeg(pbig_recep(iproc)+1)),&
         &domlen(pbig_recep(iproc)+1),MPI_INT,pbig_recep(iproc),tag,&
         &get_domain_comm(),reqrec(iproc),ierr)
   end do
   do iproc = 1, nbig_send
      tag = nproc*pbig_send(iproc) + rank + 1
      call MPI_ISEND(glob,domlen(rank+1),MPI_INT,pbig_send(iproc),&
         &tag,get_domain_comm(),reqsend(iproc),ierr)
   end do
   do iproc = 1, nbig_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
   do iproc = 1, nbig_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!
   do iproc = 1, nbig_recep
      do i = 1, domlen(pbig_recep(iproc)+1)
         loc(glob(bufbeg(pbig_recep(iproc)+1)+i-1)) =&
            &bufbeg(pbig_recep(iproc)+1)+i-1
         repart(glob(bufbeg(pbig_recep(iproc)+1)+i-1)) =&
            &pbig_recep(iproc)
      end do
   end do
!
!     also get the indexes of recdir neighboring processes
!
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         tag = nproc*rank + precdir_recep2(iproc) + 1
         call MPI_IRECV(glob(bufbeg(precdir_recep2(iproc)+1)),&
            &domlen(precdir_recep2(iproc)+1),MPI_INT,&
            &precdir_recep2(iproc),tag,get_domain_comm(),reqrec(iproc),ierr)
      end if
   end do
   do iproc = 1, nrecdir_send2
      if (precdir_send2(iproc).ne.rank) then
         tag = nproc*precdir_send2(iproc) + rank + 1
         call MPI_ISEND(glob,domlen(rank+1),MPI_INT,&
            &precdir_send2(iproc),tag,get_domain_comm(),reqsend(iproc),ierr)
      end if
   end do
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         call MPI_WAIT(reqrec(iproc),status,ierr)
      end if
   end do
   do iproc = 1, nrecdir_send2
      if (precdir_send2(iproc).ne.rank) then
         call MPI_WAIT(reqsend(iproc),status,ierr)
      end if
   end do
!
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         do i = 1, domlen(precdir_recep2(iproc)+1)
            loc(glob(bufbeg(precdir_recep2(iproc)+1)+i-1)) =&
               &bufbeg(precdir_recep2(iproc)+1)+i-1
            repart(glob(bufbeg(precdir_recep2(iproc)+1)+i-1)) =&
               &precdir_recep2(iproc)
         end do
      end if
   end do
!$acc update device(repart,loc,glob) async
!
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (ind1temp)
   deallocate (counttemp)

   if (tinkerdebug.gt.0) call MPI_BARRIER(get_host_comm(),ierr)
   if (deb_Path) write(*,*) '   << orderbuffer'
end subroutine orderbuffer

!
!     Orderbuffer with openacc directives
!
module subroutine orderbuffer_gpu(init)
   use sizes
   use atoms  , only : n
   use cutoff
   use domdec
   use inform , only : deb_Path
   use mpi
   use sizes  , only : tinkerdebug
   implicit none
   integer :: counttemp(nproc)
   integer count
   integer i,iproc,iproc1,tag,ierr
   integer iloc,iglob
   integer idomlen,iprec,ibeg
   integer :: reqsend(nproc),reqrec(nproc)
   integer status(MPI_STATUS_SIZE)
   integer cap,cpt
   logical, intent(in) :: init

   ! This approach MPI_Allgather is fine when total number of prcocess
   ! doesn't exceed a certain amount
   ! Consider using classical method with more than 64 MPI process

   if (tinkerdebug.gt.0) call MPI_BARRIER(get_host_comm(),ierr)
   if (deb_Path) write(*,*) '   >> orderbuffer_gpu',init
   ! Get domlen values of all processes
   call MPI_IAllgather( MPI_IN_PLACE, 1, MPI_DATATYPE_NULL,&
      &domlen      , 1, MPI_INT,&
      &get_domain_comm(),reqsend(1), ierr )
   call MPI_WAIT(reqsend(1),status,ierr)

   if (init) then
      cpt = 0
!$acc parallel loop present(repart,glob,loc) copyin(cpt)
      do i = 1, n
         !check if atom is in the domain or in one of the neighboring domains
         if (repart(i).eq.rank) then
!$acc atomic capture
            cpt       = cpt + 1
            cap       = cpt
!$acc end atomic
            glob(cap) = i
            loc(i)    = cap
         end if
      end do
   end if

   bufbeg(rank+1) = 1
   count = domlen(rank+1)
   do iproc = 1, nbig_recep
      if (domlen(pbig_recep(iproc)+1).ne.0) then
         bufbeg(pbig_recep(iproc)+1) = count + 1
      else
         bufbeg(pbig_recep(iproc)+1) = 1
      end if
      count = count + domlen(pbig_recep(iproc)+1)
   end do
   nbloc = count
!     nvdwblocloop is nvdwbloc if nvdwbloc is a multiple of 16, or the first one greater
   nblocloop = merge( nbloc,&
      &(int(nbloc/16)+1)*16,&
      &(mod(nbloc,16).eq.0))

   nblocrecdir = count

   do iproc = 1, nrecdir_recep1
      if (precdir_recep1(iproc).eq.rank) goto 10
      do iproc1 = 1, nbig_recep
         if (pbig_recep(iproc1).eq.precdir_recep1(iproc)) goto 10
      end do
      if (domlen(precdir_recep1(iproc)+1).ne.0) then
         bufbeg(precdir_recep1(iproc)+1) = count + 1
      else
         bufbeg(precdir_recep1(iproc)+1) = 1
      end if
      count = count + domlen(precdir_recep1(iproc)+1)
10    continue
   end do
   nblocrecdir = count

   !send and receive the indexes of the neighboring processes
!$acc host_data use_device(glob)
   do iproc = 1, nbig_recep
      tag = nproc*rank + pbig_recep(iproc) + 1
      call MPI_IRECV(glob(bufbeg(pbig_recep(iproc)+1))&
         &,domlen(pbig_recep(iproc)+1),MPI_INT&
         &,pbig_recep(iproc),tag,get_domain_comm(),reqrec(iproc),ierr)
   end do
!$acc wait
   do iproc = 1, nbig_send
      tag = nproc*pbig_send(iproc) + rank + 1
      call MPI_ISEND(glob,domlen(rank+1),MPI_INT&
         &,pbig_send(iproc),tag,get_domain_comm(),reqsend(iproc),ierr)
   end do
!$acc end host_data
   do iproc = 1, nbig_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
   do iproc = 1, nbig_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!
!     also get the indexes of recdir neighboring processes
!
!$acc host_data use_device(glob)
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         tag = nproc*rank + precdir_recep2(iproc) + 1
         call MPI_IRECV(glob(bufbeg(precdir_recep2(iproc)+1))&
            &,domlen(precdir_recep2(iproc)+1),MPI_INT&
            &,precdir_recep2(iproc),tag,get_domain_comm(),reqrec(iproc),ierr)
      end if
   end do
   do iproc = 1, nrecdir_send2
      if (precdir_send2(iproc).ne.rank) then
         tag = nproc*precdir_send2(iproc) + rank + 1
         call MPI_ISEND(glob,domlen(rank+1),MPI_INT&
            &,precdir_send2(iproc),tag,get_domain_comm(),reqsend(iproc),ierr)
      end if
   end do
!$acc end host_data
   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         call MPI_WAIT(reqrec(iproc),status,ierr)
      end if
   end do
   do iproc = 1, nrecdir_send2
      if (precdir_send2(iproc).ne.rank) then
         call MPI_WAIT(reqsend(iproc),status,ierr)
      end if
   end do

   !Update loc & repart

!$acc data present(repart,glob,loc)
   do iproc = 1, nbig_recep
      iprec   = pbig_recep(iproc)
      idomlen = domlen(iprec+1)
      ibeg    = bufbeg(iprec+1)
!$acc parallel loop async
      do i = 1, idomlen
         iglob = glob(ibeg+i-1)
         repart(iglob) = iprec
         loc   (iglob) = ibeg+i-1
      end do
   end do

   do iproc = 1, nrecdir_recep2
      if (precdir_recep2(iproc).ne.rank) then
         iprec   = precdir_recep2(iproc)
         idomlen = domlen(iprec+1)
         ibeg    = bufbeg(iprec+1)
!$acc parallel loop async
         do i = 1, idomlen
            iglob = glob(ibeg+i-1)
            repart(iglob) = iprec
            loc   (iglob) = ibeg+i-1
         end do
      end if
   end do
!$acc end data

   if (tinkerdebug.gt.0) call MPI_BARRIER(get_host_comm(),ierr)
   if (deb_Path) write(*,*) '   << orderbuffer_gpu'
end subroutine orderbuffer_gpu
!
!
!     subroutine orderbufferrec : get the indexes of the particules in the reciprocal processes
!
module subroutine orderbufferrec
   use atoms  , only: n
   use domdec
   use inform , only: deb_Path
   use potent
   use sizes  , only: tinkerdebug
   use mpi
   implicit none
   integer counttemp(nproc),reqsend(nproc),reqrec(nproc)
   integer,allocatable::ind1temp(:)
   integer count,rankloc,nprocloc
   integer i,iproc
   integer iloc,iglob,ierr,tag,status(MPI_STATUS_SIZE)
   integer iprec,ind1
!
   if (use_mlpot_only) return
!
   if (tinkerdebug.gt.0) call MPI_Barrier(MPI_COMM_WORLD,i)
   if (deb_Path) write(*,*) '   >> orderbufferrec'
!     allocate (ind1temp(n))
!
   if (use_pmecore) then
      nprocloc = nrec
      rankloc  = rank_bis
   else
      nprocloc = nproc
      rankloc  = rank
   end if
!
!      ind1temp = 0
   counttemp = 0
   domlenrec = 0
!$acc wait
!
   if ((.not.(use_pmecore)).or.((use_pmecore).and.(rank.gt.ndir-1)))&
      &then
      do iproc = 1, nrecdir_recep
         do i = 1, domlen(precdir_recep(iproc)+1)
            iloc  = bufbeg(precdir_recep(iproc)+1)+i-1
            iglob = glob(iloc)
            domlenrec(repartrec(iglob)+1) =&
               &domlenrec(repartrec(iglob)+1) + 1
            ind1  = domlenrec(repartrec(iglob)+1)
!
!       check if atom is in the domain or in one of the neighboring domains
!
            if (repartrec(iglob).eq.rankloc) then
               globrec(ind1) = iglob
               locrec(iglob) = ind1
            end if
         end do
      end do
   end if
!
   if (use_pmecore) then
      nlocrec = domlenrec(rank_bis+1)
   else
      nlocrec = domlenrec(rank+1)
   end if
!
!     get the domlen values of the coupled reciprocal space processes
!
   do iproc = 1, nrec_recep
      tag = nprocloc*rankloc + prec_recep(iproc) + 1
      call MPI_IRECV(domlenrec(prec_recep(iproc)+1),1,MPI_INT,&
         &prec_recep(iproc),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(iproc),ierr)
   end do
   do iproc = 1, nrec_send
      tag = nprocloc*prec_send(iproc) + rankloc + 1
      call MPI_ISEND(nlocrec,1,MPI_INT,prec_send(iproc),tag,&
         &merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(iproc),ierr)
   end do
   do iproc = 1, nrec_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
   do iproc = 1, nrec_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!
   bufbegrec(rankloc+1) = 1
   count = nlocrec

   do iproc = 1, nrec_recep
      if (domlenrec(prec_recep(iproc)+1).ne.0) then
         bufbegrec(prec_recep(iproc)+1) = count + 1
      else
         bufbegrec(prec_recep(iproc)+1) = 1
      end if
      count = count + domlenrec(prec_recep(iproc)+1)
10    continue
   end do
   nlocrec2 = count
!
!     send and receive the indexes of the neighboring reciprocal processes
!
   do iproc = 1, nrec_recep
      tag = nprocloc*rankloc + prec_recep(iproc) + 1
      call MPI_IRECV(globrec(bufbegrec(prec_recep(iproc)+1)),&
         &domlenrec(prec_recep(iproc)+1),MPI_INT,prec_recep(iproc),tag,&
         &merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(iproc),ierr)
   end do
   do iproc = 1, nrec_send
      tag = nprocloc*prec_send(iproc) + rankloc + 1
      call MPI_ISEND(globrec,domlenrec(rankloc+1),MPI_INT,&
         &prec_send(iproc),&
         &tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(iproc),ierr)
   end do
   do iproc = 1, nrec_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
   do iproc = 1, nrec_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!
   do iproc = 1, nrec_recep
      do i = 1, domlenrec(prec_recep(iproc)+1)
         locrec(globrec(bufbegrec(prec_recep(iproc)+1)+i-1)) =&
            &bufbegrec(prec_recep(iproc)+1)+i-1
         repartrec(globrec(bufbegrec(prec_recep(iproc)+1)+i-1)) =&
            &prec_recep(iproc)
      end do
   end do
!$acc update device(repartrec)
!
!     count = 0
!     bufbegrec(precdir_recep(1)+1) = 1
!     if (nrecdir_recep.gt.0) then
!       count = domlen(precdir_recep(1)+1)
!     end if
!     do iproc = 2, nrecdir_recep
!       if (domlen(precdir_recep(iproc)+1).ne.0) then
!         bufbegrec(precdir_recep(iproc)+1) = count + 1
!       else
!         bufbegrec(precdir_recep(iproc)+1) = 1
!       end if
!       count = count + domlen(precdir_recep(iproc)+1)
!     end do
!     nblocrec = count
!     do iproc = 1, nrecdir_recep
!       do i = 1, domlen(precdir_recep(iproc)+1)
!         globrec1(bufbegrec(precdir_recep(iproc)+1)+i-1) =
!    $       glob(bufbeg(precdir_recep(iproc)+1)+i-1)
!         locrec1(globrec1(bufbegrec(precdir_recep(iproc)+1)+i-1)) =
!    $      bufbegrec(precdir_recep(iproc)+1)+i-1
!       end do
!     end do
!$acc update device(globrec,locrec)
!
!     deallocate (ind1temp)

   if (tinkerdebug.gt.0) call MPI_Barrier(MPI_COMM_WORLD,i)
   if (deb_Path) write(*,*) '   << orderbufferrec'
end subroutine orderbufferrec

module subroutine orderbufrec_gpu
   use atoms    ,only : n
   use domdec
   use inform   ,only : deb_Path
   use potent
   use pme      ,only : GridDecomp1d
   use mpi
   use sizes    ,only : tinkerdebug
   implicit none
!     integer, allocatable :: counttemp(:),ind1temp(:)
   integer count,rankloc,nprocloc
   integer i,iproc
   integer idomlen,ibeg,ibeg1,iprec,irep,ind1
   integer iloc,iglob
   integer ierr,tag,status(MPI_STATUS_SIZE)
   integer, allocatable :: reqsend(:),reqrec(:)
!
   if (use_mlpot_only) return

   allocate(reqsend(max(1,nproc)),reqrec(max(1,nproc)))
   reqsend = MPI_REQUEST_NULL
   reqrec = MPI_REQUEST_NULL
!
   if (tinkerdebug.gt.0) call MPI_Barrier(MPI_COMM_WORLD,i)
   if (deb_Path) write(*,*) '   >> orderbufrec_gpu'
!
   if (use_pmecore) then
      nprocloc = nrec
      rankloc = rank_bis
   else
      nprocloc = nproc
      rankloc = rank
   end if
!

   ! Debug check
   !if (tinkerdebug.gt.0) then
!$acc parallel loop async
   do i = 1,n
      locrec(i) = 0
   end do
   !end if

!$acc parallel loop async
   do i = 1,nproc
      domlenrec(i)=0
   end do
!
   if ((.not.(use_pmecore)).or.((use_pmecore).and.(rank.gt.ndir-1)))&
      &then
      do iproc = 1, nrecdir_recep
         iprec   = precdir_recep(iproc)
         idomlen = domlen(iprec+1)
         ibeg    = bufbeg(iprec+1)
!$acc parallel loop async
         do i = 1, idomlen
            iloc  = ibeg+i-1
            iglob = glob(iloc)
            irep  = repartrec(iglob)+1
!$acc atomic capture
            domlenrec(irep) = domlenrec(irep) + 1
            ind1            = domlenrec(irep)
!$acc end atomic
!
!       check if atom is in the domain or in one of the neighboring domains
!
            if (irep.eq.rankloc+1) then
               globrec(ind1) = iglob
               locrec(iglob) = ind1
            end if
         end do
      end do
   end if
!$acc update host(domlenrec) async
!
!$acc wait
   ! Update nlocrec
   if (use_pmecore) then
      nlocrec = domlenrec(rank_bis+1)
   else
      nlocrec = domlenrec(rank+1)
   end if
!
!     get the domlen values of the coupled reciprocal space processes
!
   do iproc = 1, nrec_recep
      tag = nprocloc*rankloc + prec_recep(iproc) + 1
      call MPI_IRECV(domlenrec(prec_recep(iproc)+1),1,MPI_INT&
         &,prec_recep(iproc),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(iproc),ierr)
   end do
   do iproc = 1, nrec_send
      tag = nprocloc*prec_send(iproc) + rankloc + 1
      call MPI_ISEND(nlocrec,1,MPI_INT,prec_send(iproc),tag&
         &,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(iproc),ierr)
   end do
!
   do iproc = 1, nrec_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
!$acc update device(domlenrec) async
   do iproc = 1, nrec_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!
   bufbegrec(rankloc+1) = 1
   count = nlocrec

   do iproc = 1, nrec_recep
      if (domlenrec(prec_recep(iproc)+1).ne.0) then
         bufbegrec(prec_recep(iproc)+1) = count + 1
      else
         bufbegrec(prec_recep(iproc)+1) = 1
      end if
      count = count + domlenrec(prec_recep(iproc)+1)
10    continue
   end do
   nlocrec2 = count
!
!     send and receive the indexes of the neighboring reciprocal processes
!
!$acc host_data use_device(globrec)
   do iproc = 1, nrec_recep
      tag = nprocloc*rankloc + prec_recep(iproc) + 1
      call MPI_IRECV(globrec(bufbegrec(prec_recep(iproc)+1))&
         &,domlenrec(prec_recep(iproc)+1),MPI_INT&
         &,prec_recep(iproc),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(iproc),ierr)
   end do
   do iproc = 1, nrec_send
      tag = nprocloc*prec_send(iproc) + rankloc + 1
      call MPI_ISEND(globrec,domlenrec(rankloc+1),MPI_INT&
         &,prec_send(iproc),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(iproc),ierr)
   end do
   do iproc = 1, nrec_recep
      call MPI_WAIT(reqrec(iproc),status,ierr)
   end do
   do iproc = 1, nrec_send
      call MPI_WAIT(reqsend(iproc),status,ierr)
   end do
!$acc end host_data
!
!$acc data present(globrec,locrec,repartrec)
   do iproc = 1, nrec_recep
      iprec  = prec_recep(iproc)
      ibeg   = bufbegrec(iprec+1)-1
      idomlen= domlenrec(iprec+1)
!$acc parallel loop async
      do i = 1, idomlen
         iglob            = globrec(ibeg+i)
         locrec(iglob)    = ibeg+i
         repartrec(iglob) = iprec
      end do
   end do
!$acc end data
   !if (tinkerdebug.gt.0) call MPI_Barrier(MPI_COMM_WORLD,i)
   !if (deb_Path) write(0,*) __FILE__,__LINE__

!     call MPI_BARRIER(hostcomm,ierr)
!46   format(I3,'domlen',4I8,'bufbeg',2I8)
!     print 46, rank,domlen,domlenrec,bufbegrec

   !Debug Stuff
   if (tinkerdebug.gt.0) then
      if (.not.use_pmecore .or. rank.ge.ndir) then
         call MPI_AllReduce(nlocrec,reqsend(1),1,MPI_INT&
            &,MPI_SUM,merge(get_recip_comm(),get_domain_comm(),use_pmecore),ierr)
         if (reqsend(1).ne.n.and.&
            &((use_pmecore.and.rank.eq.ndir).or.&
            &(.not.use_pmecore.and.rank.eq.0)) ) then
48          format('orderbufrec_gpu::Error',/&
               &,' nlocrec reduction',I8,'.ne. natoms',I9)
            print 48, reqsend(1),n
         end if
      end if
      call MPI_BARRIER(get_host_comm(),ierr)

      count=0
!$acc wait
!$acc parallel loop present(globrec,locrec)
      do i = 1,nlocrec2
         iglob = globrec(i)
         if (locrec(iglob).ne.i) then
            print*,'error locrec',rank,i,iglob,locrec(iglob)
            count = count + 1
         end if
      end do

   end if

   if (tinkerdebug.gt.0) call MPI_Barrier(MPI_COMM_WORLD,i)
   if (deb_Path) write(*,*) '   << orderbufrec_gpu'
   deallocate(reqsend,reqrec)
end subroutine orderbufrec_gpu

   module subroutine exchange_polymer_positions(polymer,send_centroid,fast)
      use atomsMirror
      use domdec
      use parallel, only: get_domain_comm
      use inform   ,only: deb_Path
      use freeze
      use mpi
      use potent
      use timestat ,only: timer_enter,timer_exit,timer_commstep&
         &,quiet_timers
      use parallel ,only: domain_workspace_rp
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      logical, intent(in) :: fast
      integer i,k,j,tag,ierr,iproc
      integer iloc,iglob,idomlen,ibufbeg
      real(r_p),pointer:: buffer(:,:,:),buffers(:,:,:)
      real(r_p),pointer:: buffer_rp1(:),buffer_rp2(:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer,dimension(nproc):: ptemp_recep,ptemp_send
      integer status(MPI_STATUS_SIZE)
      integer ntemp_recep,ntemp_send
      integer nbeadsproc,ibeadbeg,ibeadend,nbeadssend
      integer cshift

      if (nproc.eq.1.or.(use_pmecore).and.(rank.gt.ndir-1)) return

      if (fast) then
         ptemp_recep = pneig_recep
         ptemp_send  = pneig_send
         ntemp_recep = nneig_recep
         ntemp_send  = nneig_send
      else
         ptemp_recep = pbig_recep
         ptemp_send  = pbig_send
         ntemp_recep = nbig_recep
         ntemp_send  = nbig_send
      end if

      if(send_centroid) then
         cshift=1
      else
         cshift=0
      endif

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1
      nbeadssend = nbeadsproc+cshift
!
!
!     communicate positions
!
      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commpos'

      !Fetch memory from pool
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
!$acc wait
      call domain_workspace_rp(1,3*nbloc*nbeadssend,buffer_rp1)
      call domain_workspace_rp(2,3*nloc*nbeadssend,buffer_rp2)
      buffer (1:nbeadssend,1:3,1:nbloc)&
         &=>buffer_rp1(1:3*nbloc*nbeadssend)
      buffers(1:nbeadssend,1:3,1:nloc)&
         &=>buffer_rp2(1:3*nloc*nbeadssend)
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(buffer)
      do i = 1, ntemp_recep
         tag = nproc*rank + ptemp_recep(i) + 1
         call MPI_IRECV(buffer(1,1,bufbeg(ptemp_recep(i)+1)),&
            &3*domlen(ptemp_recep(i)+1)*nbeadssend,&
            &MPI_RPREC,ptemp_recep(i),tag,get_domain_comm(),reqrec(i),ierr)
      end do
!$acc end host_data
!
!     MPI : move in buffer
!
!$acc parallel loop async default(present) collapse(3)
      do i = 1, nloc; do j=1,3; do k=1,nbeadsproc
               iglob = glob(i)
               if(k==cshift) then
                  buffers(1,j,i) = polymer%eigpos(j,iglob,1)
               endif
               buffers(k+cshift,j,i) = polymer%pos(j,iglob,k+ibeadbeg-1)
            end do ; enddo ; enddo
!
!     MPI : begin sending
!
!$acc wait
!$acc host_data use_device(buffers)
      do i = 1, ntemp_send
         tag = nproc*ptemp_send(i) + rank + 1
         call MPI_ISEND(buffers,3*nloc*nbeadssend&
            &,MPI_RPREC,&
            &ptemp_send(i),tag,get_domain_comm(),&
            &reqsend(i),ierr)
      end do
!$acc end host_data
!
      do i = 1, ntemp_recep
         call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, ntemp_send
         call MPI_WAIT(reqsend(i),status,ierr)
      end do
!
!     MPI : move in global arrays
!
      do iproc = 1, ntemp_recep
         idomlen = domlen(ptemp_recep(iproc)+1)
         ibufbeg = bufbeg(ptemp_recep(iproc)+1)
!$acc parallel loop async default(present)  collapse(2)
         do i = 1, idomlen; do j=1,3
               iloc   = ibufbeg+i-1
               iglob  = glob(iloc)
               if(cshift==1) then
                  polymer%eigpos(j,iglob,1)=buffer(1,j,iloc)
               endif
!$acc loop seq
               do k = 1, nbeadsproc
                  polymer%pos(j,iglob,k+ibeadbeg-1)&
                     &= buffer(k+cshift,j,iloc)
               enddo
            end do; enddo
      end do
!
      nullify(buffer)
      nullify(buffers)
      deallocate (reqsend)
      deallocate (reqrec)

      if (deb_Path) write(*,*) '   << commpos'
      call timer_exit( timer_commstep,quiet_timers )
   end subroutine exchange_polymer_positions

   module subroutine exchange_polymer_recip_positions(polymer,send_centroid)
      use atomsMirror
      use domdec
      use parallel, only: get_domain_comm, get_recip_comm
      use inform   ,only: deb_Path
      use mpi
      use potent   ,only: use_pmecore
      use pme      ,only: GridDecomp1d
      use timestat ,only: timer_enter,timer_exit,timer_commstep&
         &,quiet_timers
      use parallel ,only: domain_workspace_rp
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      integer i,iproc,iprec
      integer iglob,iloc,j,k
      integer tag,ierr,ibufbeg,idomlen
      integer rankloc,nprocloc
      integer status(MPI_STATUS_SIZE)
      real(r_p),pointer :: buffer(:,:,:),buffers(:,:,:)
      real(r_p),pointer :: buffer_rp1(:),buffer_rp2(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer nbeadsproc,ibeadbeg,ibeadend,nbeadssend
      integer cshift

      if (nproc.eq.1) return
      if(send_centroid) then
         cshift=1
      else
         cshift=0
      endif

!
      if (GridDecomp1d.and.Bdecomp1d.and..not.use_pmecore) return

      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commposrecpi'

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1
      nbeadssend = nbeadsproc+cshift
      allocate(reqsend(max(1,nproc)),reqrec(max(1,nproc)))
      reqsend = MPI_REQUEST_NULL
      reqrec = MPI_REQUEST_NULL
!
      if (use_pmecore) then
         nprocloc = nrec
         rankloc  = rank_bis
      else
         nprocloc = nproc
         rankloc  = rank
      end if
!
!     allocate some arrays
!
!$acc wait
      call domain_workspace_rp(1,3*max(nblocrecdir,nlocrec2)*nbeadssend,buffer_rp1)
      call domain_workspace_rp(2,3*max(nloc,nlocrec)*nbeadssend,buffer_rp2)
      buffer (1:nbeadssend,1:3,1:nblocrecdir)&
         &=>buffer_rp1(1:3*nblocrecdir*nbeadssend)
      buffers(1:nbeadssend,1:3,1:nloc)&
         &=>buffer_rp2(1:3*nloc*nbeadssend)
!
!     first do dir rec communications
!
!
!     MPI : begin reception in buffer
!
!$acc data present(x,y,z,glob,globrec)
!$acc data present(buffer,buffers)

!$acc host_data use_device(buffer,buffers)
      do i = 1, nrecdir_recep2
         if (precdir_recep2(i).ne.rank) then
            tag = nproc*rank + precdir_recep2(i) + 1
            call MPI_IRECV(buffer(1,1,bufbeg(precdir_recep2(i)+1)),&
               &3*domlen(precdir_recep2(i)+1)*nbeadssend,MPI_RPREC,&
               &precdir_recep2(i),tag,get_domain_comm(),reqrec(i),ierr)
         end if
      end do
!
!     MPI : move in buffer
!
!$acc parallel loop async default(present) collapse(3)
      do i = 1, nloc; do j=1,3; do k=1,nbeadsproc
               iglob = glob(i)
               if(k==cshift) then
                  buffers(1,j,i) = polymer%eigpos(j,iglob,1)
               endif
               buffers(k+cshift,j,i) = polymer%pos(j,iglob,k+ibeadbeg-1)
            end do ; enddo ; enddo
!
!     MPI : begin sending
!
!$acc wait
      do i = 1, nrecdir_send2
         if (precdir_send2(i).ne.rank) then
            tag = nproc*precdir_send2(i) + rank + 1
            call MPI_ISEND(buffers,3*nloc*nbeadssend&
               &,MPI_RPREC,precdir_send2(i),&
               &tag,get_domain_comm(),reqsend(i),ierr)
         end if
      end do
!$acc end host_data
!
!     MPI : move in global arrays
!
      do iproc = 1, nrecdir_recep2
         if (precdir_recep2(iproc).ne.rank) then
            call MPI_WAIT(reqrec(iproc),status,ierr)
            ibufbeg = bufbeg(precdir_recep2(iproc)+1)
            idomlen = domlen(precdir_recep2(iproc)+1)
!$acc parallel loop async collapse(3) default(present)
            do i = 1, idomlen; do j=1,3; do k=1,nbeadsproc
                     iloc     = ibufbeg+i-1
                     iglob    = glob(iloc)
                     if(k==cshift) then
                        polymer%eigpos(j,iglob,1) = buffer(1,j,iloc)
                     endif
                     polymer%pos(j,iglob,k+ibeadbeg-1) = buffer(k+cshift,j,iloc)
                  enddo; enddo; enddo
         end if
      end do
!
      do i = 1, nrecdir_send2
         if (precdir_send2(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
         end if
      end do
!
!$acc end data
!$acc wait
      nullify(buffer)
      nullify(buffers)
!
!     then do rec rec communications
!
      buffer (1:nbeadssend,1:3,1:nlocrec2)&
         &=>buffer_rp1(1:3*nlocrec2*nbeadssend)
      buffers(1:nbeadssend,1:3,1:nlocrec)&
         &=>buffer_rp2(1:3*nlocrec*nbeadssend)

!$acc data present(buffer,buffers)
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(buffer,buffers)
      do i = 1, nrec_recep
         tag = nprocloc*rankloc + prec_recep(i) + 1
         call MPI_IRECV(buffer(1,1,bufbegrec(prec_recep(i)+1))&
            &,3*domlenrec(prec_recep(i)+1)*nbeadssend,MPI_RPREC&
            &,prec_recep(i),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(i),ierr)
      end do
!
!     MPI : move in buffer
!
!$acc parallel loop async default(present) collapse(3)
      do i = 1, nlocrec; do j=1,3; do k=1,nbeadsproc
               iglob = globrec(i)
               if(k==cshift) then
                  buffers(1,j,i) = polymer%eigpos(j,iglob,1)
               endif
               buffers(k+cshift,j,i) = polymer%pos(j,iglob,k+ibeadbeg-1)
            end do ; enddo ; enddo
!
!     MPI : begin sending
!
!$acc wait
      do i = 1, nrec_send
         tag = nprocloc*prec_send(i) + rankloc+1
         call MPI_ISEND(buffers,3*nlocrec*nbeadssend&
            &,MPI_RPREC,prec_send(i)&
            &,tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(i),ierr)
      end do
!$acc end host_data
!
!     MPI : move in global arrays
!
      do iproc = 1, nrec_recep
         call MPI_WAIT(reqrec(iproc),status,ierr)
         ibufbeg = bufbegrec(prec_recep(iproc)+1)
         idomlen = domlenrec(prec_recep(iproc)+1)
!$acc parallel loop async collapse(3) default(present)
         do i = 1, idomlen; do j=1,3; do k=1,nbeadsproc
                  iloc     = ibufbeg+i-1
                  iglob    = globrec(iloc)
                  if(k==cshift) then
                     polymer%eigpos(j,iglob,1) = buffer(1,j,iloc)
                  endif
                  polymer%pos(j,iglob,k+ibeadbeg-1) = buffer(k+cshift,j,iloc)
               enddo; enddo; enddo
      end do
      ! Wait for sending requests
      do i = 1, nrec_send
         call MPI_WAIT(reqsend(i),status,ierr)
      end do

!$acc end data
!$acc end data
      nullify(buffer)
      nullify(buffers)
      deallocate(reqsend,reqrec)

      if (deb_Path) write(*,*) '   << commposrecpi'
      call timer_exit( timer_commstep,quiet_timers )
   end subroutine exchange_polymer_recip_positions

   module subroutine exchange_polymer_gradient_state(polymer,send_centroid,send_vel)
      use beads
      use mpi
      use domdec
      use parallel, only: get_replica_comm
      use inform, only: deb_path
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      logical, intent(in), optional :: send_vel
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real(r_p) :: sqrtnu
      integer :: nbeadslocmax, nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      integer :: cshift,nelt
      integer :: send_vel_
      real(r_p), allocatable :: bufferpi_s(:,:,:)
      real(r_p), allocatable :: bufferpi_r(:,:,:)

      if(nproc_polymer==1) return
      reqs(:)=MPI_REQUEST_NULL

      if(deb_path) then
         write(*,*) '>> exchange_polymer_gradient_state'
      endif
      if(send_centroid) then
         cshift=1
      else
         cshift=0
      endif

      send_vel_=.true.
      if(present(send_vel)) send_vel_=send_vel
      if(allocated(polymer%vel) .and. send_vel_) then
         nelt=2
      else
         nelt=1
      endif

      nu=polymer%nbeads
      sqrtnu=sqrt(real(nu))

      if(pi_comm_report) then
         write(pi_comm_unit,*) "  entering exchange_polymer_gradient_state"
         FLUSH(pi_comm_unit)
      endif

      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(bufferpi_s(nelt*3,maxloc*(nbeadslocmax+cshift)&
         &,nproc_polymer))
      allocate(bufferpi_r(nelt*3,maxloc*(nbeadslocmax+cshift)&
         &,nproc_polymer))

!$acc enter data create(bufferpi_r,bufferpi_s) async

      !PUT DATA INTO BUFFER ARRAYS

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)

!$acc parallel loop collapse(2) async &
!$acc present(bufferpi_s)
      DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
            iproc=polymer%iproc(k)
            if(rank_polymer==iproc-1) CYCLE
            ibeadbeg=polymer%ibead_beg(iproc)
            i=glob(ilocbeg+iloc)
            if(k-ibeadbeg==0 .and. send_centroid) then
               kk=iloc +1
!$acc loop seq
               DO j=1,3
                  bufferpi_s(j,kk,iproc) = polymer%eigpos(j,i,1)
               ENDDO
               if(nelt==2) then
!$acc loop seq
                  DO j=1,3
                     bufferpi_s(3+j,kk,iproc) = polymer%eigvel(j,i,1)
                  ENDDO
               endif
            endif
            kk=(k-ibeadbeg + cshift)*nlocpi+ iloc +1
!$acc loop seq
            DO j=1,3
               bufferpi_s(j,kk,iproc)=polymer%pos(j,i,k)
            ENDDO
            if(nelt==2) then
!$acc loop seq
               DO j=1,3
                  bufferpi_s(3+j,kk,iproc)=polymer%vel(j,i,k)
               ENDDO
            endif
         ENDDO; ENDDO !; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
         if(rank_polymer==iproc-1) CYCLE
         ii=ii+1
         ibeadbeg = polymer%ibead_beg(iproc)
         ibeadend = polymer%ibead_end(iproc)
         if(pi_comm_report) then
            write(pi_comm_unit,*) "  send"&
               &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)&
               &,"ats to",iproc-1
         endif
         call MPI_ISEND(bufferpi_s(1,1,iproc)&
            &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)&
            &,MPI_RPREC&
            &,iproc-1, mptag_data, get_replica_comm(), reqs(ii),ierr)
      ENDDO

      ibeadbeg = polymer%ibead_beg(rank_polymer+1)
      ibeadend = polymer%ibead_end(rank_polymer+1)
      DO iproc=1,nproc_polymer
         if(rank_polymer==iproc-1) CYCLE
         ii=ii+1
         ilocbeg=ilocpi_beg(iproc)
         ilocend=ilocpi_end(iproc)
         if(pi_comm_report) then
            write(pi_comm_unit,*) "  receive"&
               &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)&
               &,"ats from",iproc-1
         endif
         call MPI_IRECV(bufferpi_r(1,1,iproc)&
            &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)&
            &,MPI_RPREC,iproc-1,mptag_data, get_replica_comm()&
            &,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      !WAIT FOR COMMUNICATIONS
      if(pi_comm_report) then
         write(pi_comm_unit,*) "  Wait for data"
         FLUSH(pi_comm_unit)
      endif
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)


      !PUT DATA INTO POLYMER ARRAYS

!$acc parallel loop collapse(2) async &
!$acc present(bufferpi_r,ilocpi_beg,ilocpi_end)
      DO k=1-cshift,ibeadend-ibeadbeg+1; DO iloc=1,nloc !; DO j=1,3
            iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
            if(rank_polymer==iproc-1) CYCLE
            ilocbeg=ilocpi_beg(iproc)
            ilocend=ilocpi_end(iproc)
            i=glob(iloc)
            if(k==0) then
               kk=iloc-ilocbeg + 1
!$acc loop seq
               DO j=1,3
                  polymer%eigpos(j,i,1)=bufferpi_r(j,kk,iproc)
               ENDDO
               if(nelt==2) then
!$acc loop seq
                  DO j=1,3
                     polymer%eigvel(j,i,1)=bufferpi_r(3+j,kk,iproc)
                  ENDDO
               endif
            else
               ibead=k+ibeadbeg-1
               kk=(k+cshift-1)*(ilocend-ilocbeg+1)+ iloc-ilocbeg +1
!$acc loop seq
               DO j=1,3
                  polymer%pos(j,i,ibead)=bufferpi_r(j,kk,iproc)
               ENDDO
               if(nelt==2) then
!$acc loop seq
                  DO j=1,3
                     polymer%vel(j,i,ibead)=bufferpi_r(3+j,kk,iproc)
                  ENDDO
               endif
            endif
         ENDDO; ENDDO !; ENDDO

!!$acc wait
!      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r) async
      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
         write(pi_comm_unit,*) "  exiting exchange_polymer_gradient_state"
         FLUSH(pi_comm_unit)
      endif

      if(deb_path) then
         write(*,*) '<< exchange_polymer_gradient_state'
      endif

   end subroutine exchange_polymer_gradient_state
   module subroutine exchange_polymer_normal_modes(polymer&
      &,array1,array2,array3,array4)
      use beads
      use mpi
      use domdec
      use parallel, only: get_replica_comm
      use inform, only: deb_path
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer
      real(r_p), intent(inout) :: array1(:,:,:)
      real(r_p), intent(inout), optional :: array2(:,:,:)
      real(r_p), intent(inout), optional :: array3(:,:,:)
      real(r_p), intent(inout), optional :: array4(:,:,:)
      logical :: send_pos_, send_vel_,slow_
      logical :: correctorder
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real(r_p) :: sqrtnu
      integer :: nbeadslocmax,nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg,nelt,ielt
      real(r_p), allocatable :: bufferpi_s(:,:,:)
      real(r_p), allocatable :: bufferpi_r(:,:,:)

      if(nproc_polymer==1) return

      if(deb_path) then
!$acc wait
         write(*,*) '>> exchange_polymer_normal_modes'
      endif

      reqs(:)=MPI_REQUEST_NULL

      nu=polymer%nbeads
      sqrtnu=sqrt(real(nu,r_p))

      if(pi_comm_report) then
         write(pi_comm_unit,*) "  entering exchange_polymer_normal_modes"
         FLUSH(pi_comm_unit)
      endif

      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      nelt=1
      if(present(array2)) then
         nelt=nelt+1
      endif
      if(present(array3)) then
         if(.not. present(array2)) then
            write(*,*) 'Error in exchange_polymer_normal_modes: array3'&
               &,' is present but array2 is not'
            call fatal
         endif
         nelt=nelt+1
      endif
      if(present(array4)) then
         if(.not. present(array3)) then
            write(*,*) 'Error in exchange_polymer_normal_modes: array4'&
               &,' is present but array3 is not'
            call fatal
         endif
         nelt=nelt+1
      endif
      allocate(bufferpi_s(nelt*3,maxloc*nbeadslocmax,nproc_polymer))
      allocate(bufferpi_r(nelt*3,maxloc*nbeadslocmax,nproc_polymer))

!      allocate(bufferpi_s(3*3,maxloc*nbeadslocmax,nproc_polymer))
!      allocate(bufferpi_r(3*3,maxloc*nbeadslocmax,nproc_polymer))

!$acc enter data create(bufferpi_r,bufferpi_s) async

      !PUT DATA INTO BUFFER ARRAYS
      ibeadbeg = polymer%ibead_beg(rank_polymer+1)
      ibeadend = polymer%ibead_end(rank_polymer+1)
!$acc parallel loop collapse(2) async &
!$acc present(bufferpi_s,ilocpi_beg,ilocpi_end)
      DO k=0,ibeadend-ibeadbeg; DO iloc=1,nloc !; DO j=1,3
            iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
            if(rank_polymer==iproc-1) cycle
            ilocbeg=ilocpi_beg(iproc)
            ilocend=ilocpi_end(iproc)
            i=glob(iloc)
            ibead=k+ibeadbeg
            kk=k*(ilocend-ilocbeg+1)+ iloc-ilocbeg + 1
!$acc loop seq
            do j=1,3
               bufferpi_s(  j,kk,iproc)=array1(j,i,ibead)
            enddo
            if(nelt==1) cycle
!$acc loop seq
            do j=1,3
               bufferpi_s(3+j,kk,iproc)=array2(j,i,ibead)
            enddo
            if(nelt==2) cycle
!$acc loop seq
            do j=1,3
               bufferpi_s(6+j,kk,iproc)=array3(j,i,ibead)
            enddo
            if(nelt==3) cycle
!$acc loop seq
            do j=1,3
               bufferpi_s(9+j,kk,iproc)=array4(j,i,ibead)
            enddo
         ENDDO; ENDDO !; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
         if(rank_polymer==iproc-1) CYCLE
         ii=ii+1
         ilocbeg=ilocpi_beg(iproc)
         ilocend=ilocpi_end(iproc)
         if(pi_comm_report) then
            write(pi_comm_unit,*) "  send"&
               &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
               &,"ats to",iproc-1
         endif
         call MPI_ISEND(bufferpi_s(1,1,iproc)&
            &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
            &,MPI_RPREC&
            &,iproc-1, mptag_data, get_replica_comm(), reqs(ii),ierr)
      ENDDO

      !RECEIVE DATA
      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      DO iproc=1,nproc_polymer
         if(rank_polymer==iproc-1) CYCLE
         ii=ii+1
         ibeadbeg = polymer%ibead_beg(iproc)
         ibeadend = polymer%ibead_end(iproc)
         if(pi_comm_report) then
            write(pi_comm_unit,*) "  receive"&
               &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
               &,"ats to",iproc-1
         endif
         call MPI_IRECV(bufferpi_r(1,1,iproc)&
            &,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
            &,MPI_RPREC,iproc-1,mptag_data, get_replica_comm()&
            &,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      !WAIT FOR COMMUNICATIONS
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      !PUT DATA INTO POLYMER ARRAYS

!$acc parallel loop collapse(2) async &
!$acc present(polymer,bufferpi_r,polymer%ibead_beg,polymer%iproc)
      DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
            iproc=polymer%iproc(k)
            if(rank_polymer==iproc-1) cycle
            ibeadbeg=polymer%ibead_beg(iproc)
            i=glob(ilocbeg+iloc)
            kk=(k-ibeadbeg)*nlocpi + iloc + 1
!$acc loop seq
            do j=1,3
               array1(j,i,k)=bufferpi_r(j,kk,iproc)
            enddo
            if(nelt==1) cycle
!$acc loop seq
            do j=1,3
               array2(j,i,k)=bufferpi_r(3+j,kk,iproc)
            enddo
            if(nelt==2) cycle
!$acc loop seq
            do j=1,3
               array3(j,i,k)=bufferpi_r(6+j,kk,iproc)
            enddo
            if(nelt==3) cycle
!$acc loop seq
            do j=1,3
               array4(j,i,k)=bufferpi_r(9+j,kk,iproc)
            enddo
         ENDDO; ENDDO !; ENDDO

!!$acc wait
!      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r) async
      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
         write(pi_comm_unit,*) "  exiting exchange_polymer_normal_modes"
         FLUSH(pi_comm_unit)
      endif

      if(deb_path) then
         write(*,*) '<< exchange_polymer_normal_modes'
      endif

   end subroutine exchange_polymer_normal_modes
   module subroutine exchange_polymer_recip_forces(polymer,derivs,temprec)
      use deriv
      use domdec
      use parallel, only: get_domain_comm, get_recip_comm
      use mpi
      use potent ,pa=>PotentialAll
      use pme    ,only: GridDecomp1d
      use inform ,only: deb_Path
      use sizes
      use timestat ,only: timer_enter,timer_exit,quiet_timers&
         &,timer_dirreccomm,timer_fmanage
      use parallel ,only: domain_workspace_rp
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      integer :: nbeadsproc
      real(r_p), intent(inout) :: derivs(:,:,:)
      real(r_p), intent(inout) :: temprec(:,:,:)
      integer i,j,k,tag,ierr,iproc,iglob,iloc,ilocrec,ibufbeg
      integer ibeadbeg,ibeadend
      integer jloc,jglob,jlocrec,ii
      integer sz1,sz2,sz3
      integer rankloc,nprocloc
      integer reqsend(nproc),reqrec(nproc)
      integer status(MPI_STATUS_SIZE)
      real(r_p),pointer:: buffer(:,:,:,:)
      real(r_p),pointer:: buffers(:,:,:),temprecbuf(:,:,:)
      real(r_p),pointer:: buffer_rp1(:),buffer_rp2(:),buffer_rp3(:)

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc=ibeadend-ibeadbeg+1

      !case sequential
      if (nproc.eq.1) then
!$acc parallel loop collapse(3) &
!$acc         default(present) async
         do k= 1,nbeadsproc; do i = 1, nloc
               do j = 1, 3
                  derivs(j,i,k) = derivs(j,i,k) + temprec(j,i,k)
               end do
            end do; end do

         return
      end if

      if (GridDecomp1d.and.Bdecomp1d.and..not.use_pmecore) then
!$acc parallel loop collapse(3) &
!$acc         default(present) async
         do k= 1,nbeadsproc; do i = 1, nbloc
               do j = 1, 3
                  ilocrec = locrec(glob(i))
                  if (ilocrec.eq.0) cycle
                  derivs(j,i,k) = derivs(j,i,k) + temprec(j,ilocrec,k)
               end do
            end do; end do

         return
      end if

      if (deb_Path) write(*,*) '   >> exchange_polymer_recip_forces'
      call timer_enter( timer_dirreccomm )
!
      if (use_pmecore) then
         nprocloc = nrec
         rankloc  = rank_bis
      else
         nprocloc = nproc
         rankloc  = rank
      end if
!
!     allocate some arrays
!
!$acc wait
      sz1=3*max(max(1,nloc)*nrecdir_send1,nlocrec*nrec_send)*nbeadsproc
      sz2=3*max(1,nblocrecdir)*nbeadsproc
      sz3=3*max(1,nlocrec2)*nbeadsproc
      call domain_workspace_rp(1,sz1,buffer_rp1)
      call domain_workspace_rp(2,sz2,buffer_rp2)
      call domain_workspace_rp(3,sz3,buffer_rp3)

      sz1 = 3*max(1,nlocrec)*nrec_send
      buffer(1:nbeadsproc,1:3,1:max(1,nlocrec),1:nrec_send)&
         &=> buffer_rp1(1:sz1)
      buffers(1:nbeadsproc,1:3,1:max(1,nblocrecdir))&
         &=> buffer_rp2(1:sz2)
      temprecbuf(1:nbeadsproc,1:3,1:max(1,nlocrec2))&
         &=> buffer_rp3(1:sz3)
!$acc enter data attach(buffer,buffers)
!
!     first do rec rec communications
!
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(buffer)
      do i = 1, nrec_send
         tag = nprocloc*rankloc + prec_send(i) + 1
         call MPI_IRECV(buffer(1,1,1,i),3*nlocrec*nbeadsproc,MPI_RPREC&
            &,prec_send(i),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqrec(i),ierr)
      end do
!$acc end host_data
!
!     MPI : move in buffer
!
!$acc data present(temprec,temprecbuf,buffers,derivs)

!$acc parallel loop collapse(3) &
!$acc         default(present) async
      do k= 1,nbeadsproc; do i = 1, nlocrec2
            do j = 1, 3
               temprecbuf(k,j,i) = temprec(j,i,k)
            end do
         end do; end do
!
!$acc wait
!$acc host_data use_device(temprec)
      do i = 1, nrec_recep
         tag = nprocloc*prec_recep(i) + rankloc + 1
         call MPI_ISEND(temprecbuf(1,1,bufbegrec(prec_recep(i)+1))&
            &,3*domlenrec(prec_recep(i)+1)*nbeadsproc&
            &,MPI_RPREC&
            &,prec_recep(i),tag,merge(get_recip_comm(),get_domain_comm(),use_pmecore),reqsend(i),ierr)
      end do
!$acc end host_data
!
      do i = 1, nrec_send
         call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_recep
         call MPI_WAIT(reqsend(i),status,ierr)
      end do
!
!     mpi : move in global arrays
!
!$acc parallel loop collapse(3) present(buffer)
      do k=1,nbeadsproc; do i = 1,nlocrec
            do j = 1,3
!$acc loop seq
               do ii = 1,nrec_send
                  temprec(j,i,k) = temprec(j,i,k) + buffer(k,j,i,ii)
               end do
            end do
         end do; enddo

!$acc exit data detach(buffer)
      sz1 = 3*max(1,nloc)*nrecdir_send1*nbeadsproc
      nullify(buffer)
      buffer(1:nbeadsproc,1:3,1:max(1,nloc),1:nrecdir_send1)&
         &=> buffer_rp1(1:sz1)
!$acc enter data attach(buffer)
!
!     then do rec dir communications
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(buffer)
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
            tag = nproc*rank + precdir_send1(i) + 1
            call MPI_IRECV(buffer(1,1,1,i),3*nloc*nbeadsproc,MPI_RPREC&
               &,precdir_send1(i),tag,get_domain_comm(),reqrec(i),ierr)
         end if
      end do
!$acc end host_data
!
!     Move in array
!
      do ii = 1, nrecdir_recep1
         if (precdir_recep1(ii).ne.rank) then
            ibufbeg =bufbeg(precdir_recep1(ii)+1)-1
!$acc parallel loop collapse(3) &
!$acc    present(repartrec,locrec,glob)
            do k=1,nbeadsproc; do i = 1, domlen(precdir_recep1(ii)+1)
                  do j = 1,3
                     iloc = ibufbeg+i
                     iglob = glob(iloc)
                     if (repartrec(iglob).eq.rankloc) then
                        ilocrec = locrec(iglob)
                        buffers(k,j,iloc) = temprec(j,ilocrec,k)
                     else
                        buffers(k,j,iloc) = 0
                     end if
                  end do
               end do; enddo
         end if
      end do
!
!$acc host_data use_device(buffers)
      do i = 1, nrecdir_recep1
         if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buffers(1,1,bufbeg(precdir_recep1(i)+1))&
               &,3*domlen(precdir_recep1(i)+1)*nbeadsproc,MPI_RPREC&
               &,precdir_recep1(i),tag,get_domain_comm(),reqsend(i),ierr)
         end if
      end do
!$acc end host_data
!
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
         end if
      end do
!
!     MPI : move in global arrays
!
      do ii = 1, nrecdir_send1
         if (precdir_send1(ii).ne.rank) then
!$acc parallel loop collapse(3) async present(buffer)
            do k=1,nbeadsproc; do i = 1,nloc
                  do j = 1,3
                     derivs(j,i,k) = derivs(j,i,k) + buffer(k,j,i,ii)
                  end do
               end do; enddo
         else
!$acc parallel loop collapse(3) async
            do k=1,nbeadsproc; do i = 1, nlocrec
                  do j = 1,3
                     iglob = globrec(i)
                     iloc  = loc(iglob)
                     if (repart(iglob).eq.rank) then
                        derivs(j,iloc,k) = derivs(j,iloc,k) + temprec(j,i,k)
                     end if
                  end do
               end do; enddo
         end if
      end do

      do i = 1, nrecdir_recep1
         if (precdir_recep1(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
         end if
      end do

!$acc wait
!$acc end data
!$acc exit data detach(buffer,buffers,temprecbuf)
      nullify(buffer)
      nullify(buffers)
      nullify(temprecbuf)

      call timer_exit( timer_dirreccomm,quiet_timers )
      if (deb_Path) write(*,*) '   << exchange_polymer_recip_forces'


   end subroutine exchange_polymer_recip_forces
   module subroutine exchange_polymer_contracted_forces(polymer_ctr,derivs_ctr)
      use beads
      use mpi
      use domdec
      use parallel, only: get_replica_comm
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer_ctr
      real(r_p), intent(inout) :: derivs_ctr(:,:,:)
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real(r_p) :: sqrtnu
      integer :: nbeadslocmax,nlocproc,nbeadsproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      real(r_p), allocatable :: bufferpi_s(:,:,:)
      real(r_p), allocatable :: bufferpi_r(:,:,:)

      nu=polymer_ctr%nbeads

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      ibeadbeg = polymer_ctr%ibead_beg(rank_polymer+1)
      ibeadend = polymer_ctr%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

!$acc parallel loop collapse(3) async &
!$acc present(polymer_ctr,glob,derivs_ctr)
      DO k=1,nbeadsproc; DO iloc=ilocbeg,ilocend ; DO j=1,3
               i=glob(iloc)
               polymer_ctr%forces_slow(j,i,k+ibeadbeg-1)&
                  &=-derivs_ctr(j,iloc,k)
            ENDDO; ENDDO ; ENDDO

      if(nproc_polymer==1)   return

      reqs(:)=MPI_REQUEST_NULL

      if(pi_comm_report) then
         write(pi_comm_unit,*) "  entering exchange_polymer_contracted_forces"
         FLUSH(pi_comm_unit)
      endif

      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer_ctr%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(bufferpi_s(3,maxloc*nbeadslocmax,nproc_polymer))
      allocate(bufferpi_r(3,maxloc*nbeadslocmax,nproc_polymer))

!$acc enter data create(bufferpi_r,bufferpi_s) async

      !PUT DATA INTO BUFFER ARRAYS

!$acc parallel loop collapse(2) async &
!$acc present(polymer_ctr,bufferpi_s,ilocpi_beg,ilocpi_end)
      DO k=1,ibeadend-ibeadbeg+1; DO iloc=1,nloc !; DO j=1,3
            iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
            if(rank_polymer/=iproc-1) then
               ilocbeg=ilocpi_beg(iproc)
               ilocend=ilocpi_end(iproc)
               kk=(k-1)*(ilocend-ilocbeg+1)+ iloc-ilocbeg + 1
!$acc loop seq
               do j=1,3
                  bufferpi_s(j,kk,iproc)=derivs_ctr(j,iloc,k)
               enddo
            endif
         ENDDO; ENDDO !; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
         if(rank_polymer==iproc-1) CYCLE
         ii=ii+1
         ilocbeg=ilocpi_beg(iproc)
         ilocend=ilocpi_end(iproc)
         if(pi_comm_report) then
            write(pi_comm_unit,*) "  send"&
               &,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
               &,"ats to",iproc-1
         endif
         call MPI_ISEND(bufferpi_s(1,1,iproc)&
            &,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
            &,MPI_RPREC&
            &,iproc-1, mptag_data, get_replica_comm(), reqs(ii),ierr)
      ENDDO

      !RECEIVE DATA
      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      DO iproc=1,nproc_polymer
         if(rank_polymer==iproc-1) CYCLE
         ii=ii+1
         ibeadbeg = polymer_ctr%ibead_beg(iproc)
         ibeadend = polymer_ctr%ibead_end(iproc)
         if(pi_comm_report) then
            write(pi_comm_unit,*) "  receive"&
               &,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
               &,"ats to",iproc-1
         endif
         call MPI_IRECV(bufferpi_r(1,1,iproc)&
            &,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)&
            &,MPI_RPREC,iproc-1,mptag_data, get_replica_comm()&
            &,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      !WAIT FOR COMMUNICATIONS
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      !PUT DATA INTO POLYMER ARRAYS

!$acc parallel loop collapse(2) async &
!$acc present(polymer_ctr,bufferpi_r,polymer_ctr%ibead_beg,polymer_ctr%iproc)
      DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
            iproc=polymer_ctr%iproc(k)
            if(rank_polymer/=iproc-1) then
               ibeadbeg=polymer_ctr%ibead_beg(iproc)
               i=glob(ilocbeg+iloc)
               kk=(k-ibeadbeg)*nlocpi + iloc + 1
!$acc loop seq
               do j=1,3
                  polymer_ctr%forces_slow(j,i,k)=-bufferpi_r(j,kk,iproc)
               enddo
            endif
         ENDDO; ENDDO !; ENDDO

!!$acc wait
!      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r) async
      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
         write(pi_comm_unit,*) "  exiting exchange_polymer_contracted_forces"
         FLUSH(pi_comm_unit)
      endif

   end subroutine exchange_polymer_contracted_forces

end submodule parallel_communication
