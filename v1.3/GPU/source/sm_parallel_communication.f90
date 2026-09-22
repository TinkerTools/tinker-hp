!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (parallel) parallel_communication
   use domdec
   use potent, only: use_pmecore
   use utilgpu, only: def_queue
   use mpi
   implicit none

   type(domain_t), pointer, save :: active_domain => null()
   integer, save :: active_polar_ticket_count = 0

contains

module subroutine activate_domain(platform)
   type(parallel_t), target, intent(inout) :: platform

   active_domain => platform%domain
end subroutine activate_domain

module subroutine deactivate_domain(platform)
   type(parallel_t), target, intent(inout) :: platform

   if (associated(active_domain, platform%domain)) nullify(active_domain)
end subroutine deactivate_domain

module subroutine migrate_atoms_legacy()

   if (associated(active_domain)) call check_domain_layout()
   call reassign
   if (associated(active_domain)) call refresh_domain_layout(active_domain)
end subroutine migrate_atoms_legacy

module subroutine polar_post_comm_legacy(ticket, nrhs, values)
   type(polar_ticket_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(:,:,:)

   if (ticket%active) error stop 'polar exchange ticket is already active'
   allocate(ticket%receive_requests(max(1,nproc)), ticket%send_requests(max(1,nproc)))
   ticket%receive_requests = MPI_REQUEST_NULL
   ticket%send_requests = MPI_REQUEST_NULL
   ticket%active = .true.
   call commdirdirgpu(nrhs, 0, values, ticket%receive_requests, ticket%send_requests)
   call commdirdirgpu(nrhs, 1, values, ticket%receive_requests, ticket%send_requests)
   active_polar_ticket_count = active_polar_ticket_count + 1
end subroutine polar_post_comm_legacy

module subroutine polar_finish_legacy(ticket, nrhs, values)
   type(polar_ticket_t), intent(inout) :: ticket
   integer, intent(in) :: nrhs
   real(t_p), intent(inout) :: values(:,:,:)

   if (.not.ticket%active) error stop 'polar exchange ticket is not active'
   call commdirdirgpu(nrhs, 2, values, ticket%receive_requests, ticket%send_requests)
   call reset_polar_ticket(ticket)
end subroutine polar_finish_legacy

module subroutine reset_polar_ticket(ticket)
   type(polar_ticket_t), intent(inout) :: ticket

   if (allocated(ticket%receive_requests)) deallocate(ticket%receive_requests)
   if (allocated(ticket%send_requests)) deallocate(ticket%send_requests)
   if (ticket%active .and. active_polar_ticket_count > 0) &
      &active_polar_ticket_count = active_polar_ticket_count - 1
   ticket%active = .false.
end subroutine reset_polar_ticket

module function active_polar_tickets() result(count)
   integer :: count

   count = active_polar_ticket_count
end function active_polar_tickets

module subroutine capture_halo_layout(domain)
   type(domain_t), intent(inout) :: domain

   call release_halo_layout(domain)
   domain%nproc = nproc
   domain%nloc = nloc
   domain%nbig_send = nbig_send
   domain%nbig_recep = nbig_recep
   domain%nneig_send = nneig_send
   domain%nneig_recep = nneig_recep
   allocate(domain%domlen(nproc), domain%bufbeg(nproc))
   allocate(domain%pbig_send(max(1,nproc)), domain%pbig_recep(max(1,nproc)))
   allocate(domain%pneig_send(max(1,nproc)), domain%pneig_recep(max(1,nproc)))
   allocate(domain%reqs(max(1,nproc)), domain%reqr(max(1,nproc)))
   domain%domlen = domlen
   domain%bufbeg = bufbeg
   domain%pbig_send = pbig_send
   domain%pbig_recep = pbig_recep
   domain%pneig_send = pneig_send
   domain%pneig_recep = pneig_recep
   domain%reqs = MPI_REQUEST_NULL
   domain%reqr = MPI_REQUEST_NULL
   domain%layout_ready = .true.
end subroutine capture_halo_layout

module subroutine refresh_domain_layout(domain)
   type(domain_t), intent(inout) :: domain

   call capture_halo_layout(domain)
end subroutine refresh_domain_layout

module subroutine release_halo_layout(domain)
   type(domain_t), intent(inout) :: domain

!$acc wait
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
   if (allocated(domain%domlen)) deallocate(domain%domlen, domain%bufbeg)
   if (allocated(domain%pbig_send)) deallocate(domain%pbig_send, domain%pbig_recep)
   if (allocated(domain%pneig_send)) deallocate(domain%pneig_send, domain%pneig_recep)
   if (allocated(domain%reqs)) deallocate(domain%reqs, domain%reqr)
   domain%nproc = 0
   domain%nloc = 0
   domain%layout_ready = .false.
end subroutine release_halo_layout

module subroutine extend_legacy_tp(values, components, halo)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo

   call exchange_extend_tp(values, components, halo)
end subroutine extend_legacy_tp

module subroutine accumulate_legacy_tp(values, components, halo)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo

   call exchange_accumulate_tp(values, components, halo)
end subroutine accumulate_legacy_tp

#if TINKER_MIXED_PREC
module subroutine extend_legacy_rp(values, components, halo)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo

   call exchange_extend_rp(values, components, halo)
end subroutine extend_legacy_rp

module subroutine accumulate_legacy_rp(values, components, halo)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo

   call exchange_accumulate_rp(values, components, halo)
end subroutine accumulate_legacy_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
module subroutine extend_legacy_md(values, components, halo)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo

   call exchange_extend_md(values, components, halo)
end subroutine extend_legacy_md

module subroutine accumulate_legacy_md(values, components, halo)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo

   call exchange_accumulate_md(values, components, halo)
end subroutine accumulate_legacy_md
#endif

module subroutine accumulate_rec_to_dir_legacy_tp(direct_values, reciprocal_values, components)
   use utilcomm, only: commDDrd, ucComm
   real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
   integer, intent(in) :: components

   ! The R/D family is still a temporary compatibility adapter.  It is not
   ! part of the D/D state copied by this first migration tranche.
   call commDDrd(direct_values, reciprocal_values, components, ucComm)
end subroutine accumulate_rec_to_dir_legacy_tp

module subroutine replicate_legacy_rp(global_ids, owned_values, global_values, components)
   integer, intent(in) :: global_ids(:)
   real(r_p), intent(in) :: owned_values(:)
   real(r_p), intent(out) :: global_values(:)
   integer, intent(in) :: components
   integer :: ierr, i, local_count, max_global_id, value_count, total_count
   integer, allocatable :: counts(:), offsets(:), received_ids(:)
   real(r_p), allocatable :: received_values(:)

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
   allocate(received_ids(total_count), received_values(components*total_count))
   call MPI_Allgatherv(global_ids, local_count, MPI_INTEGER, received_ids, counts, offsets, &
      &MPI_INTEGER, active_domain%communicator, ierr)
   counts = components*counts
   offsets = components*offsets
   call MPI_Allgatherv(owned_values, components*local_count, MPI_RPREC, received_values, counts, offsets, &
      &MPI_RPREC, active_domain%communicator, ierr)
   max_global_id = 0
   if (total_count > 0) max_global_id = maxval(received_ids)
   if (size(global_values) < components*max_global_id) error stop 'replication result is too small'
   global_values = real(0.0,r_p)
   do i = 1, total_count
      value_count = components*(received_ids(i)-1)
      if (received_ids(i) <= 0 .or. value_count+components > size(global_values)) then
         error stop 'replication received an invalid global identifier'
      end if
      global_values(value_count+1:value_count+components) = received_values(components*(i-1)+1:components*i)
   end do
end subroutine replicate_legacy_rp

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

subroutine exchange_extend_tp(values, components, halo)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   integer :: i, ierr, peer, first, count, option, comm, request
   integer :: status(MPI_STATUS_SIZE)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   comm = active_domain%communicator
!$acc wait(def_queue)
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

subroutine exchange_accumulate_tp(values, components, halo)
   real(t_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
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
!$acc wait(def_queue)
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
!$acc parallel loop gang vector async(def_queue) present(values,work(1:max(1,value_count*incoming)))
   do j = 1, value_count
      do i = 1, incoming
         values(j) = values(j) + work(j+components*active_domain%nloc*(i-1))
      end do
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
end subroutine exchange_accumulate_tp

#if TINKER_MIXED_PREC
subroutine exchange_extend_rp(values, components, halo)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   integer :: i, ierr, peer, first, count, option, comm, request
   integer :: status(MPI_STATUS_SIZE)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   comm = active_domain%communicator
!$acc wait(def_queue)
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

subroutine exchange_accumulate_rp(values, components, halo)
   real(r_p), intent(inout) :: values(*)
   integer, intent(in) :: components, halo
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
!$acc wait(def_queue)
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
!$acc parallel loop gang vector async(def_queue) present(values,work(1:max(1,value_count*incoming)))
   do j = 1, value_count
      do i = 1, incoming
         values(j) = values(j) + work(j+components*active_domain%nloc*(i-1))
      end do
   end do
   do i = 1, send_count(option); call MPI_Wait(active_domain%reqs(i), status, ierr); end do
end subroutine exchange_accumulate_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
subroutine exchange_extend_md(values, components, halo)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo
   integer :: i, ierr, peer, first, count, option, comm, request
   integer :: status(MPI_STATUS_SIZE)

   call check_domain(halo)
   if (skip_domain_exchange()) return
   option = halo_option(halo)
   comm = active_domain%communicator
!$acc wait(def_queue)
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

subroutine exchange_accumulate_md(values, components, halo)
   mdyn_rtyp, intent(inout) :: values(*)
   integer, intent(in) :: components, halo
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
!$acc wait(def_queue)
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
!$acc parallel loop gang vector async(def_queue) present(values,work(1:max(1,value_count*incoming)))
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

   if (allocated(active_domain%work_tp) .and. size(active_domain%work_tp) >= required) return
!$acc wait
   if (allocated(active_domain%work_tp)) then
!$acc exit data delete(active_domain%work_tp(1:size(active_domain%work_tp)))
      deallocate(active_domain%work_tp)
   end if
   allocate(active_domain%work_tp(max(1,required)))
!$acc enter data create(active_domain%work_tp(1:max(1,required)))
end subroutine ensure_work_tp

#if TINKER_MIXED_PREC
subroutine ensure_work_rp(required)
   integer, intent(in) :: required

   if (allocated(active_domain%work_rp) .and. size(active_domain%work_rp) >= required) return
!$acc wait
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

   if (allocated(active_domain%work_md) .and. size(active_domain%work_md) >= required) return
!$acc wait
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
   if (halo /= halo_large .and. halo /= halo_neighbour) error stop 'unknown halo neighbourhood'
end subroutine check_domain

subroutine check_domain_layout()

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
      if (COMM_TINKER == MPI_COMM_NULL .or. nproc <= 0) error stop 'parallel domain is not active'
      communicator = COMM_TINKER
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

   halo_option = merge(1, 2, halo == halo_large)
end function halo_option

integer function send_count(option)
   integer, intent(in) :: option

   send_count = merge(active_domain%nbig_send, active_domain%nneig_send, option == 1)
end function send_count

integer function receive_count(option)
   integer, intent(in) :: option

   receive_count = merge(active_domain%nbig_recep, active_domain%nneig_recep, option == 1)
end function receive_count

integer function send_peer(option, index)
   integer, intent(in) :: option, index

   if (option == 1) then
      send_peer = active_domain%pbig_send(index)
   else
      send_peer = active_domain%pneig_send(index)
   end if
end function send_peer

integer function receive_peer(option, index)
   integer, intent(in) :: option, index

   if (option == 1) then
      receive_peer = active_domain%pbig_recep(index)
   else
      receive_peer = active_domain%pneig_recep(index)
   end if
end function receive_peer

end submodule parallel_communication
