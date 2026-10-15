!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
module parallel
   use tinTypes , only: POLYMER_COMM_TYPE
   use mdevice_channel, only: channel_t
   use mpi      , only: MPI_COMM_NULL
   use status_types, only: status_t
   implicit none

   private :: reject_parallel_copy, reject_domain_copy, reject_polar_ticket_copy, &
      &capture_halo_layout, refresh_domain_layout, release_halo_layout, check_domain_layout, activate_domain, &
      &deactivate_domain, migrate_atoms_legacy, polar_post_comm_legacy, &
      &polar_finish_legacy, &
      &active_domain, active_platform, &
      &polar_post_dip_comm_legacy, polar_finish_dip_comm_legacy, &
      &polar_post_receives_direct_legacy, polar_post_sends_direct_legacy, &
      &polar_finish_direct_legacy, polar_post_receives_dip_legacy, &
      &polar_post_sends_dip_legacy, polar_finish_dip_legacy, &
      &polar_post_receives_direct, polar_post_sends_direct, polar_finish_direct, &
      &polar_post_receives_dip, polar_post_sends_dip, polar_finish_dip, &
      &polar_set_finish_failure_legacy, &
      &reset_polar_ticket, replicate_rp, replicate_full_rp, replicate_legacy_rp, &
      &orderbuffer, orderbuffer_gpu, orderbufferrec, orderbufrec_gpu, &
      &extend_tp, accumulate_tp, accumulate_rec_to_dir_tp, &
      &reduce_sum_i_vector, reduce_sum_rp_vector, reduce_sum_rp_matrix, &
      &reduce_sum_legacy_i_vector, reduce_sum_legacy_rp_vector, reduce_sum_legacy_rp_matrix, &
      &extend_legacy_tp, accumulate_legacy_tp, accumulate_rec_to_dir_legacy_tp, &
      &parallel_ctf, parallel_image

#if TINKER_MIXED_PREC
   private :: extend_rp, accumulate_rp, reduce_sum_tp_vector, &
      &extend_legacy_rp, accumulate_legacy_rp, reduce_sum_legacy_tp_vector
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
   private :: extend_md, accumulate_md, reduce_sum_md_vector, &
      &extend_legacy_md, accumulate_legacy_md, reduce_sum_legacy_md_vector
#endif

   !
   ! TYPES
   !
   type :: topo_cfg_t
      logical :: replicas_enabled = .false.
      integer :: replica_count = 0
   end type topo_cfg_t

   type :: domain_t
      ! Noncopyable owner; reset_topology releases and restores its state explicitly.
      integer :: communicator = MPI_COMM_NULL
      integer :: direct_comm = MPI_COMM_NULL
      integer :: recip_comm = MPI_COMM_NULL
      integer :: size = 0
      integer :: rank = 0
      integer :: nproc = 0, nloc = 0
      integer :: nbig_send = 0, nbig_recep = 0
      integer :: nneig_send = 0, nneig_recep = 0
      integer :: nshort_send = 0, nshort_recep = 0
      integer, allocatable :: domlen(:), bufbeg(:)
      integer, allocatable :: pbig_send(:), pbig_recep(:)
      integer, allocatable :: pneig_send(:), pneig_recep(:)
      integer, allocatable :: pshort_send(:), pshort_recep(:)
      integer, private, allocatable :: reqs(:), reqr(:)
      real(t_p), allocatable :: work_tp(:)
      real(r_p), allocatable :: work_rp(:)
      real(t_p), allocatable :: recdir_receive_work(:), recdir_send_work(:)
      integer, private, allocatable :: recdir_recv_reqs(:), recdir_send_reqs(:)
      ! Communication-family workspaces.  Scientific arrays are borrowed by
      ! the operation; these pools belong to the active domain and are
      ! released with its layout.
      integer, private, allocatable :: pool_i1(:), pool_i2(:)
      real(r_p), private, allocatable :: pool_rp1(:), pool_rp2(:), pool_rp3(:)
      real(t_p), private, allocatable :: pool_t1(:), pool_t2(:), pool_t3(:), pool_t4(:)
      real(t_p), private, allocatable :: pool_t5(:), pool_t6(:)
#if defined(USE_DETERMINISTIC_REDUCTION)
      mdyn_rtyp, allocatable :: work_md(:)
#endif
      logical :: layout_ready = .false.
   end type domain_t

   type :: comm_opt_t
      ! Orchestration policy for optional direct-space halo communication.
      logical :: skip_direct = .false.
      logical :: pole_ordered = .false.
   end type comm_opt_t

   type :: polar_comm_t
      private
      ! Noncopyable owner of in-flight MPI requests; polar_finish releases it.
      integer, allocatable :: recv_req(:), send_req(:)
      integer, allocatable :: dip_recv_req(:), dip_send_req(:)
      integer, allocatable :: recdir_recv_req(:), recdir_send_req(:)
      logical :: active = .false.
      logical :: dip_active = .false.
      logical :: recdir_active = .false.
      logical :: recdir_send_phase = .false.
      logical :: recdir_receive_phase = .false.
      logical :: use_short_range = .false.
      logical :: use_group = .false.
      logical :: use_gpu = .false.
   end type polar_comm_t

   type :: host_t
      integer :: communicator = MPI_COMM_NULL
      integer :: rank = 0
   end type host_t

   type :: replica_t
      integer :: communicator = MPI_COMM_NULL
      integer :: size = 0
      integer :: rank = 0
      integer :: index = 0
   end type replica_t

   type :: root_exch_t
      integer :: communicator = MPI_COMM_NULL
   end type root_exch_t

   type :: parallel_t
      ! Noncopyable MPI topology owner. Components are loaned through
      ! the handle queries and legacy adapter; reset_topology releases every loan.
      integer :: world_size = 0
      integer :: world_rank = 0
      type(domain_t) :: domain
      type(host_t) :: host
      type(replica_t) :: replica
      type(root_exch_t) :: root_exchange
   end type parallel_t

   !
   !  Globals
   !
   integer, parameter :: topology_classical = 1, topology_polymer = 2
   integer, parameter :: halo_large = 1, halo_neighbour = 2, halo_short = 3
   integer, parameter :: reduction_root = 1, reduction_all = 2

   type(domain_t), pointer, save :: active_domain => null()
   type(parallel_t), pointer, save :: active_platform => null()
   type(comm_opt_t) :: comm_opt


   !
   ! INTERFACES
   !
   interface polar_post_comm
      module procedure polar_post_comm_default, polar_post_comm_short, polar_post_comm_send
   end interface polar_post_comm

   interface polar_post_dip_comm
      module procedure polar_post_dip_comm_default, polar_post_dip_comm_gpu
   end interface polar_post_dip_comm

   interface polar_post_receives
      module procedure polar_post_receives_direct, polar_post_receives_dip
   end interface polar_post_receives

   interface polar_post_sends
      module procedure polar_post_sends_direct, polar_post_sends_dip
   end interface polar_post_sends

   interface polar_finish
      module procedure polar_finish_direct, polar_finish_dip
   end interface polar_finish

   interface polar_finish_dip_comm
      module procedure polar_finish_dip_comm_default, polar_finish_dip_comm_gpu
   end interface polar_finish_dip_comm

   interface extend
      module procedure extend_tp
#if TINKER_MIXED_PREC
      module procedure extend_rp
#endif
#if defined(USE_DETERMINISTIC_REDUCTION)
      module procedure extend_md
#endif
   end interface extend

   interface accumulate
      module procedure accumulate_tp, accumulate_rec_to_dir_tp
#if TINKER_MIXED_PREC
      module procedure accumulate_rp
#endif
#if defined(USE_DETERMINISTIC_REDUCTION)
      module procedure accumulate_md
#endif
   end interface accumulate

   interface replicate
      module procedure replicate_rp
      module procedure replicate_full_rp
   end interface replicate

   interface reduce_sum
#if TINKER_MIXED_PREC
      module procedure reduce_sum_tp_vector
#endif
      module procedure reduce_sum_i_vector
      module procedure reduce_sum_rp_vector, reduce_sum_rp_matrix
#if defined(USE_DETERMINISTIC_REDUCTION)
      module procedure reduce_sum_md_vector
#endif
   end interface reduce_sum

   interface
      module subroutine construct_topology(platform, topo_config, topology_kind, status_value)
         type(parallel_t), target, intent(inout) :: platform
         type(topo_cfg_t), intent(in) :: topo_config
         integer, intent(in) :: topology_kind
         type(status_t), intent(out) :: status_value
      end subroutine construct_topology

      module subroutine configure_layout(platform, key_discovered, atom_count, status_value)
         type(parallel_t), target, intent(inout) :: platform
         logical, intent(in) :: key_discovered
         integer, intent(in) :: atom_count
         type(status_t), intent(out) :: status_value
      end subroutine configure_layout

      module subroutine reset_topology(platform)
         type(parallel_t), target, intent(inout) :: platform
      end subroutine reset_topology

      module subroutine release(platform, status_value)
         type(parallel_t), target, intent(inout) :: platform
         type(status_t), intent(inout) :: status_value
      end subroutine release

      module subroutine reject_parallel_copy(lhs, rhs)
         type(parallel_t), intent(out) :: lhs
         type(parallel_t), intent(in) :: rhs
      end subroutine reject_parallel_copy

      module subroutine reject_domain_copy(lhs, rhs)
         type(domain_t), intent(out) :: lhs
         type(domain_t), intent(in) :: rhs
      end subroutine reject_domain_copy

      module subroutine reject_polar_ticket_copy(lhs, rhs)
         type(polar_comm_t), intent(out) :: lhs
         type(polar_comm_t), intent(in) :: rhs
      end subroutine reject_polar_ticket_copy

      module subroutine capture_halo_layout(domain)
         type(domain_t), intent(inout) :: domain
      end subroutine capture_halo_layout

      module subroutine refresh_domain_layout(domain)
         type(domain_t), intent(inout) :: domain
      end subroutine refresh_domain_layout

      module subroutine release_halo_layout(domain, cleanup_status)
         type(domain_t), intent(inout) :: domain
         integer, intent(out), optional :: cleanup_status
      end subroutine release_halo_layout

      module subroutine migrate_polymer_atoms(polymer, polymer_ctr)
         type(POLYMER_COMM_TYPE), intent(inout) :: polymer
         type(POLYMER_COMM_TYPE), intent(inout), optional :: polymer_ctr
      end subroutine migrate_polymer_atoms

      module subroutine exchange_polymer_positions(polymer, send_centroid, fast)
         type(POLYMER_COMM_TYPE), intent(inout) :: polymer
         logical, intent(in) :: send_centroid, fast
      end subroutine exchange_polymer_positions

      module subroutine exchange_polymer_recip_positions(polymer, send_centroid)
         type(POLYMER_COMM_TYPE), intent(inout) :: polymer
         logical, intent(in) :: send_centroid
      end subroutine exchange_polymer_recip_positions

      module subroutine exchange_polymer_gradient_state(polymer,send_centroid,send_vel)
         type(POLYMER_COMM_TYPE), intent(inout) :: polymer
         logical, intent(in) :: send_centroid
         logical, intent(in), optional :: send_vel
      end subroutine exchange_polymer_gradient_state

      module subroutine exchange_polymer_normal_modes(polymer,array1,array2,array3,array4)
         type(POLYMER_COMM_TYPE), intent(inout) :: polymer
         real(r_p), intent(inout) :: array1(:,:,:)
         real(r_p), intent(inout), optional :: array2(:,:,:)
         real(r_p), intent(inout), optional :: array3(:,:,:)
         real(r_p), intent(inout), optional :: array4(:,:,:)
      end subroutine exchange_polymer_normal_modes

      module subroutine exchange_polymer_recip_forces(polymer,derivs,temprec)
         type(POLYMER_COMM_TYPE), intent(in) :: polymer
         real(r_p), intent(inout) :: derivs(:,:,:)
         real(r_p), intent(inout) :: temprec(:,:,:)
      end subroutine exchange_polymer_recip_forces

      module subroutine exchange_polymer_contracted_forces(polymer_ctr,derivs_ctr)
         type(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
         real(r_p), intent(inout) :: derivs_ctr(:,:,:)
      end subroutine exchange_polymer_contracted_forces

      module subroutine release_migration_workspace()
      end subroutine release_migration_workspace

      module subroutine domain_workspace_i(slot, required, values)
         integer, intent(in) :: slot, required
         integer, pointer, intent(out) :: values(:)
      end subroutine domain_workspace_i

      module subroutine domain_workspace_rp(slot, required, values)
         integer, intent(in) :: slot, required
         real(r_p), pointer, intent(out) :: values(:)
      end subroutine domain_workspace_rp

#if defined(USE_DETERMINISTIC_REDUCTION)
      module subroutine domain_workspace_md(required, values)
         integer, intent(in) :: required
         mdyn_rtyp, pointer, intent(out) :: values(:)
      end subroutine domain_workspace_md
#endif

      module subroutine domain_workspace_tp2d(slot, n1, n2, values)
         integer, intent(in) :: slot, n1, n2
         real(t_p), pointer, intent(out) :: values(:,:)
      end subroutine domain_workspace_tp2d

      module subroutine domain_workspace_tp3d(slot, n1, n2, n3, values)
         integer, intent(in) :: slot, n1, n2, n3
         real(t_p), pointer, intent(out) :: values(:,:,:)
      end subroutine domain_workspace_tp3d

      module subroutine activate_domain(platform)
         type(parallel_t), target, intent(inout) :: platform
      end subroutine activate_domain

      module subroutine deactivate_domain(platform)
         type(parallel_t), target, intent(inout) :: platform
      end subroutine deactivate_domain

      module subroutine create_domain_subcommunicators(direct_count, process_rank)
         integer, intent(in) :: direct_count
         integer, intent(out) :: process_rank
      end subroutine create_domain_subcommunicators

      module subroutine check_domain_layout()
      end subroutine check_domain_layout

      module subroutine migrate_atoms_legacy()
      end subroutine migrate_atoms_legacy

      module subroutine polar_post_comm_legacy(ticket, nrhs, values, use_short_range, send_values, group_comm, &
         &channel)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(3,nrhs,*)
         logical, intent(in) :: use_short_range
         real(t_p), intent(inout) :: send_values(3,nrhs,*)
         logical, intent(in) :: group_comm
         type(channel_t), intent(in) :: channel
      end subroutine polar_post_comm_legacy

      module subroutine polar_finish_legacy(ticket, nrhs, values, status_value, channel)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(3,nrhs,*)
         type(status_t), intent(out), optional :: status_value
         type(channel_t), intent(in) :: channel
      end subroutine polar_finish_legacy

      module subroutine polar_post_dip_comm_legacy(ticket, nrhs, diprec, dip, &
         &buffer_send, buffer_recv, use_gpu)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: diprec(3,nrhs,*)
         real(t_p), intent(inout) :: dip(3,nrhs,*)
         real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
         logical, intent(in) :: use_gpu
      end subroutine polar_post_dip_comm_legacy

      module subroutine polar_finish_dip_comm_legacy(ticket, nrhs, diprec, dip, &
         &buffer_send, buffer_recv, use_gpu, status_value)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: diprec(3,nrhs,*)
         real(t_p), intent(inout) :: dip(3,nrhs,*)
         real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
         logical, intent(in) :: use_gpu
         type(status_t), intent(out), optional :: status_value
      end subroutine polar_finish_dip_comm_legacy

      module subroutine polar_post_receives_direct_legacy(ticket, nrhs, values, channel)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(3,nrhs,*)
         type(channel_t), intent(in) :: channel
      end subroutine polar_post_receives_direct_legacy

      module subroutine polar_post_sends_direct_legacy(ticket, nrhs, values, channel)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(3,nrhs,*)
         type(channel_t), intent(in) :: channel
      end subroutine polar_post_sends_direct_legacy

      module subroutine polar_finish_direct_legacy(ticket, nrhs, values, status_value, channel)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(3,nrhs,*)
         type(status_t), intent(out), optional :: status_value
         type(channel_t), intent(in) :: channel
      end subroutine polar_finish_direct_legacy

      module subroutine polar_post_receives_dip_legacy(ticket, nrhs, diprec, dip, &
         &buffer_send, buffer_recv)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
         real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      end subroutine polar_post_receives_dip_legacy

      module subroutine polar_post_sends_dip_legacy(ticket, nrhs, diprec, dip, &
         &buffer_send, buffer_recv)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
         real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      end subroutine polar_post_sends_dip_legacy

      module subroutine polar_finish_dip_legacy(ticket, nrhs, diprec, dip, &
         &buffer_send, buffer_recv, status_value)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
         real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
         type(status_t), intent(out), optional :: status_value
      end subroutine polar_finish_dip_legacy

      module subroutine polar_set_finish_failure(enabled)
         logical, intent(in) :: enabled
      end subroutine polar_set_finish_failure

      module subroutine polar_set_finish_failure_legacy(enabled)
         logical, intent(in) :: enabled
      end subroutine polar_set_finish_failure_legacy

      module subroutine polar_post_recdir(ticket, nrhs, values, receive_buffer, send_buffer, channel, use_gpu, send_phase)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(in) :: values(3,nrhs,*)
         real(t_p), intent(inout) :: receive_buffer(3,nrhs,*), send_buffer(3,nrhs,*)
         type(channel_t), intent(in) :: channel
         logical, intent(in), optional :: use_gpu
         logical, intent(in), optional :: send_phase
      end subroutine polar_post_recdir

      module subroutine polar_finish_recdir(ticket, nrhs, values, buffer, channel, use_gpu)
         type(polar_comm_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(3,nrhs,*)
         real(t_p), intent(in) :: buffer(3,nrhs,*)
         type(channel_t), intent(in) :: channel
         logical, intent(in), optional :: use_gpu
      end subroutine polar_finish_recdir

      module subroutine polar_post_recdir_fields(ticket, values, receive_buffer, send_buffer, channel, use_gpu, send_phase)
         type(polar_comm_t), intent(inout) :: ticket
         real(t_p), intent(in) :: values(10,*)
         real(t_p), intent(inout) :: receive_buffer(10,*), send_buffer(10,*)
         type(channel_t), intent(in) :: channel
         logical, intent(in) :: use_gpu
         logical, intent(in), optional :: send_phase
      end subroutine polar_post_recdir_fields

      module subroutine polar_finish_recdir_fields(ticket, values, buffer, channel, use_gpu)
         type(polar_comm_t), intent(inout) :: ticket
         real(t_p), intent(inout) :: values(10,*)
         real(t_p), intent(in) :: buffer(10,*)
         type(channel_t), intent(in) :: channel
         logical, intent(in) :: use_gpu
      end subroutine polar_finish_recdir_fields

      module subroutine commfield(nrhs, ef)
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: ef(3,nrhs,*)
      end subroutine commfield

      module subroutine commfieldgpu(nrhs, ef)
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: ef(3,nrhs,*)
      end subroutine commfieldgpu

      module subroutine commfieldshort(nrhs, ef, channel)
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: ef(3,nrhs,*)
         type(channel_t), intent(in) :: channel
      end subroutine commfieldshort

      module subroutine commchgflx(pdelta)
         real(t_p), intent(inout) :: pdelta(*)
      end subroutine commchgflx

      module subroutine commpot(pot, rule)
         real*8, intent(inout) :: pot(*)
         integer, intent(in) :: rule
      end subroutine commpot

      module subroutine commpotsum(pot, rule)
         real*8, intent(inout) :: pot(*)
         integer, intent(in) :: rule
      end subroutine commpotsum

      module subroutine commpotrec(pot, potrec)
         real*8, intent(inout) :: pot(*), potrec(*)
      end subroutine commpotrec

      module subroutine commpos()
      end subroutine commpos

      module subroutine commposrec()
      end subroutine commposrec

      module subroutine commposrespa(fast)
         logical, intent(in) :: fast
      end subroutine commposrespa

      module subroutine reassign_respa(ialt, nalt)
         integer, intent(in) :: ialt, nalt
      end subroutine reassign_respa

      module subroutine parallel_ctf(xr, yr, zr)
         real(t_p), intent(inout) :: xr, yr, zr
      end subroutine parallel_ctf

      module subroutine parallel_image(xr, yr, zr)
         real(t_p), intent(inout) :: xr, yr, zr
      end subroutine parallel_image

      module subroutine gather_positions()
      end subroutine gather_positions

      module subroutine assign_direct()
      end subroutine assign_direct

      module subroutine initialize_direct_communication()
      end subroutine initialize_direct_communication

      module subroutine assign_reciprocal()
      end subroutine assign_reciprocal

      module subroutine orderbuffer(init)
         logical, intent(in) :: init
      end subroutine orderbuffer

      module subroutine orderbuffer_gpu(init)
         logical, intent(in) :: init
      end subroutine orderbuffer_gpu

      module subroutine orderbufferrec()
      end subroutine orderbufferrec

      module subroutine orderbufrec_gpu()
      end subroutine orderbufrec_gpu

      module subroutine reset_polar_ticket(ticket)
         type(polar_comm_t), intent(inout) :: ticket
      end subroutine reset_polar_ticket

      module function active_polar_tickets() result(count)
         integer :: count
      end function active_polar_tickets

      module subroutine replicate_legacy_rp(global_ids, owned_values, global_values, components)
         integer, intent(in) :: global_ids(:)
         real(r_p), intent(in) :: owned_values(:)
         real(r_p), intent(out) :: global_values(:)
         integer, intent(in) :: components
      end subroutine replicate_legacy_rp

      module subroutine replicate_full_rp(global_ids, global_values, components)
         integer, intent(in) :: global_ids(:)
         real(r_p), intent(inout) :: global_values(:)
         integer, intent(in) :: components
      end subroutine replicate_full_rp

#if TINKER_MIXED_PREC
      module subroutine reduce_sum_legacy_tp_vector(values, visibility, root)
         real(t_p), intent(inout) :: values(:)
         integer, intent(in) :: visibility
         integer, intent(in), optional :: root
      end subroutine reduce_sum_legacy_tp_vector
#endif

      module subroutine reduce_sum_legacy_i_vector(values, visibility, root)
         integer, intent(inout) :: values(:)
         integer, intent(in) :: visibility
         integer, intent(in), optional :: root
      end subroutine reduce_sum_legacy_i_vector

      module subroutine reduce_sum_legacy_rp_vector(values, visibility, root)
         real(r_p), intent(inout) :: values(:)
         integer, intent(in) :: visibility
         integer, intent(in), optional :: root
      end subroutine reduce_sum_legacy_rp_vector

      module subroutine reduce_sum_legacy_rp_matrix(values, visibility, root)
         real(r_p), intent(inout) :: values(:,:)
         integer, intent(in) :: visibility
         integer, intent(in), optional :: root
      end subroutine reduce_sum_legacy_rp_matrix

#if defined(USE_DETERMINISTIC_REDUCTION)
      module subroutine reduce_sum_legacy_md_vector(values, visibility, root)
         mdyn_rtyp, intent(inout) :: values(:)
         integer, intent(in) :: visibility
         integer, intent(in), optional :: root
      end subroutine reduce_sum_legacy_md_vector
#endif

   end interface

   interface assignment(=)
      module procedure reject_parallel_copy, reject_domain_copy, reject_polar_ticket_copy
   end interface

   interface
      module subroutine extend_legacy_tp(values, components, halo, channel)
         real(t_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
         type(channel_t), intent(in) :: channel
      end subroutine extend_legacy_tp

      module subroutine accumulate_legacy_tp(values, components, halo, channel)
         real(t_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
         type(channel_t), intent(in) :: channel
      end subroutine accumulate_legacy_tp

      module subroutine accumulate_rec_to_dir_legacy_tp(direct_values, reciprocal_values, components, channel)
         real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
         integer, intent(in) :: components
         type(channel_t), intent(in) :: channel
      end subroutine accumulate_rec_to_dir_legacy_tp

#if TINKER_MIXED_PREC
      module subroutine extend_legacy_rp(values, components, halo, channel)
         real(r_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
         type(channel_t), intent(in) :: channel
      end subroutine extend_legacy_rp

      module subroutine accumulate_legacy_rp(values, components, halo, channel)
         real(r_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
         type(channel_t), intent(in) :: channel
      end subroutine accumulate_legacy_rp

#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
      module subroutine extend_legacy_md(values, components, halo, channel)
         mdyn_rtyp, intent(inout) :: values(*)
         integer, intent(in) :: components, halo
         type(channel_t), intent(in) :: channel
      end subroutine extend_legacy_md

      module subroutine accumulate_legacy_md(values, components, halo, channel)
         mdyn_rtyp, intent(inout) :: values(*)
         integer, intent(in) :: components, halo
         type(channel_t), intent(in) :: channel
      end subroutine accumulate_legacy_md
#endif
   end interface
contains

   integer function domain_handle(platform)
      type(parallel_t), intent(in) :: platform
      domain_handle = platform%domain%communicator
   end function domain_handle

   integer function domain_size(platform)
      type(parallel_t), intent(in) :: platform
      domain_size = platform%domain%size
   end function domain_size

   integer function host_handle(platform)
      type(parallel_t), intent(in) :: platform
      host_handle = platform%host%communicator
   end function host_handle

   integer function replica_handle(platform)
      type(parallel_t), intent(in) :: platform
      replica_handle = platform%replica%communicator
   end function replica_handle

   integer function replica_size(platform)
      type(parallel_t), intent(in) :: platform
      replica_size = platform%replica%size
   end function replica_size

   integer function replica_rank(platform)
      type(parallel_t), intent(in) :: platform
      replica_rank = platform%replica%rank
   end function replica_rank

   integer function root_exchange_handle(platform)
      type(parallel_t), intent(in) :: platform
      root_exchange_handle = platform%root_exchange%communicator
   end function root_exchange_handle

   integer function get_domain_comm()
      if (associated(active_platform)) then
         get_domain_comm = active_platform%domain%communicator
      else
         get_domain_comm = MPI_COMM_NULL
      end if
   end function get_domain_comm

   integer function get_host_comm()
      if (associated(active_platform)) then
         get_host_comm = active_platform%host%communicator
      else
         get_host_comm = MPI_COMM_NULL
      end if
   end function get_host_comm

   integer function get_replica_comm()
      if (associated(active_platform)) then
         get_replica_comm = active_platform%replica%communicator
      else
         get_replica_comm = MPI_COMM_NULL
      end if
   end function get_replica_comm

   integer function get_root_comm()
      if (associated(active_platform)) then
         get_root_comm = active_platform%root_exchange%communicator
      else
         get_root_comm = MPI_COMM_NULL
      end if
   end function get_root_comm

   integer function get_direct_comm()
      if (associated(active_platform)) then
         get_direct_comm = active_platform%domain%direct_comm
      else
         get_direct_comm = MPI_COMM_NULL
      end if
   end function get_direct_comm

   integer function get_recip_comm()
      if (associated(active_platform)) then
         get_recip_comm = active_platform%domain%recip_comm
      else
         get_recip_comm = MPI_COMM_NULL
      end if
   end function get_recip_comm

   subroutine extend_tp(values, components, halo, channel)
      real(t_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo
      type(channel_t), intent(in) :: channel

      call extend_legacy_tp(values, components, halo, channel)
   end subroutine extend_tp

#if TINKER_MIXED_PREC
   subroutine extend_rp(values, components, halo, channel)
      real(r_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo
      type(channel_t), intent(in) :: channel

      call extend_legacy_rp(values, components, halo, channel)
   end subroutine extend_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
   subroutine extend_md(values, components, halo, channel)
      mdyn_rtyp, intent(inout) :: values(*)
      integer, intent(in) :: components, halo
      type(channel_t), intent(in) :: channel

      call extend_legacy_md(values, components, halo, channel)
   end subroutine extend_md
#endif

   subroutine accumulate_tp(values, components, halo, channel)
      real(t_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo
      type(channel_t), intent(in) :: channel

      call accumulate_legacy_tp(values, components, halo, channel)
   end subroutine accumulate_tp

#if TINKER_MIXED_PREC
   subroutine accumulate_rp(values, components, halo, channel)
      real(r_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo
      type(channel_t), intent(in) :: channel

      call accumulate_legacy_rp(values, components, halo, channel)
   end subroutine accumulate_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
   subroutine accumulate_md(values, components, halo, channel)
      mdyn_rtyp, intent(inout) :: values(*)
      integer, intent(in) :: components, halo
      type(channel_t), intent(in) :: channel

      call accumulate_legacy_md(values, components, halo, channel)
   end subroutine accumulate_md
#endif

   subroutine accumulate_rec_to_dir_tp(direct_values, reciprocal_values, components, channel)
      real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
      integer, intent(in) :: components
      type(channel_t), intent(in) :: channel

      call accumulate_rec_to_dir_legacy_tp(direct_values, reciprocal_values, components, channel)
   end subroutine accumulate_rec_to_dir_tp

   subroutine migrate_atoms()

      call migrate_atoms_legacy()
   end subroutine migrate_atoms

   subroutine polar_post_comm_default(ticket, nrhs, values, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*)
      type(channel_t), intent(in) :: channel
      call polar_post_comm_legacy(ticket, nrhs, values, .false., values, .false., channel)
   end subroutine polar_post_comm_default

   subroutine polar_post_comm_short(ticket, nrhs, values, use_short_range, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*)
      logical, intent(in) :: use_short_range
      type(channel_t), intent(in) :: channel
      call polar_post_comm_legacy(ticket, nrhs, values, use_short_range, values, .false., channel)
   end subroutine polar_post_comm_short

   subroutine polar_post_comm_send(ticket, nrhs, values, use_short_range, send_values, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*), send_values(3,nrhs,*)
      logical, intent(in) :: use_short_range
      type(channel_t), intent(in) :: channel
      call polar_post_comm_legacy(ticket, nrhs, values, use_short_range, send_values, .false., channel)
   end subroutine polar_post_comm_send

   subroutine polar_post_group_comm(ticket, nrhs, values, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*)
      type(channel_t), intent(in) :: channel
      call polar_post_comm_legacy(ticket, nrhs, values, .false., values, .true., channel)
   end subroutine polar_post_group_comm

   subroutine polar_finish_direct(ticket, nrhs, values, status_value, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*)
      type(status_t), intent(out), optional :: status_value
      type(channel_t), intent(in) :: channel

      if (ticket%use_group .or. ticket%use_short_range .or. ticket%use_gpu) then
         call polar_finish_legacy(ticket, nrhs, values, status_value, channel)
      else
         call polar_finish_direct_legacy(ticket, nrhs, values, status_value, channel)
      end if
   end subroutine polar_finish_direct

   subroutine polar_post_receives_direct(ticket, nrhs, values, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*)
      type(channel_t), intent(in) :: channel

      call polar_post_receives_direct_legacy(ticket, nrhs, values, channel)
   end subroutine polar_post_receives_direct

   subroutine polar_post_sends_direct(ticket, nrhs, values, channel)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(3,nrhs,*)
      type(channel_t), intent(in) :: channel

      call polar_post_sends_direct_legacy(ticket, nrhs, values, channel)
   end subroutine polar_post_sends_direct

   subroutine polar_post_receives_dip(ticket, nrhs, diprec, dip, buffer_send, buffer_recv)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)

      call polar_post_receives_dip_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv)
   end subroutine polar_post_receives_dip

   subroutine polar_post_sends_dip(ticket, nrhs, diprec, dip, buffer_send, buffer_recv)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)

      call polar_post_sends_dip_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv)
   end subroutine polar_post_sends_dip

   subroutine polar_finish_dip(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, status_value)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      type(status_t), intent(out), optional :: status_value

      call polar_finish_dip_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, status_value)
   end subroutine polar_finish_dip

   subroutine polar_set_finish_failure(enabled)
      logical, intent(in) :: enabled

      call polar_set_finish_failure_legacy(enabled)
   end subroutine polar_set_finish_failure

   subroutine polar_post_dip_comm_default(ticket, nrhs, diprec, dip, buffer_send, buffer_recv)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      call polar_post_dip_comm_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, .true.)
   end subroutine polar_post_dip_comm_default

   subroutine polar_post_dip_comm_gpu(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, use_gpu)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      logical, intent(in) :: use_gpu
      call polar_post_dip_comm_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, use_gpu)
   end subroutine polar_post_dip_comm_gpu

   subroutine polar_finish_dip_comm_default(ticket, nrhs, diprec, dip, buffer_send, buffer_recv)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      call polar_finish_dip_comm_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, .true.)
   end subroutine polar_finish_dip_comm_default

   subroutine polar_finish_dip_comm_gpu(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, use_gpu)
      type(polar_comm_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: diprec(3,nrhs,*), dip(3,nrhs,*)
      real(t_p), intent(inout) :: buffer_send(3,nrhs,*), buffer_recv(3,nrhs,*)
      logical, intent(in) :: use_gpu
      call polar_finish_dip_comm_legacy(ticket, nrhs, diprec, dip, buffer_send, buffer_recv, use_gpu)
   end subroutine polar_finish_dip_comm_gpu

   subroutine replicate_rp(global_ids, owned_values, global_values, components)
      integer, intent(in) :: global_ids(:)
      real(r_p), intent(in) :: owned_values(:)
      real(r_p), intent(out) :: global_values(:)
      integer, intent(in) :: components

      call replicate_legacy_rp(global_ids, owned_values, global_values, components)
   end subroutine replicate_rp

#if TINKER_MIXED_PREC
   subroutine reduce_sum_tp_vector(values, visibility, root)
      real(t_p), intent(inout) :: values(:)
      integer, intent(in) :: visibility
      integer, intent(in), optional :: root

      call reduce_sum_legacy_tp_vector(values, visibility, root)
   end subroutine reduce_sum_tp_vector
#endif

   subroutine reduce_sum_i_vector(values, visibility, root)
      integer, intent(inout) :: values(:)
      integer, intent(in) :: visibility
      integer, intent(in), optional :: root

      call reduce_sum_legacy_i_vector(values, visibility, root)
   end subroutine reduce_sum_i_vector

   subroutine reduce_sum_rp_vector(values, visibility, root)
      real(r_p), intent(inout) :: values(:)
      integer, intent(in) :: visibility
      integer, intent(in), optional :: root

      call reduce_sum_legacy_rp_vector(values, visibility, root)
   end subroutine reduce_sum_rp_vector

   subroutine reduce_sum_rp_matrix(values, visibility, root)
      real(r_p), intent(inout) :: values(:,:)
      integer, intent(in) :: visibility
      integer, intent(in), optional :: root

      call reduce_sum_legacy_rp_matrix(values, visibility, root)
   end subroutine reduce_sum_rp_matrix

#if defined(USE_DETERMINISTIC_REDUCTION)
   subroutine reduce_sum_md_vector(values, visibility, root)
      mdyn_rtyp, intent(inout) :: values(:)
      integer, intent(in) :: visibility
      integer, intent(in), optional :: root

      call reduce_sum_legacy_md_vector(values, visibility, root)
   end subroutine reduce_sum_md_vector
#endif

end module parallel
