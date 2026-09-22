!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
module parallel
   use config, only: config_t
   use mpi, only: MPI_COMM_NULL
   use status, only: status_t
   use structure, only: structure_t
   implicit none
   private

   type, public :: domain_t
      private
      ! Noncopyable owner; reset_topology releases and restores its state explicitly.
      integer :: communicator = MPI_COMM_NULL
      integer :: size = 0
      integer :: rank = 0
      integer :: nproc = 0, nloc = 0
      integer :: nbig_send = 0, nbig_recep = 0
      integer :: nneig_send = 0, nneig_recep = 0
      integer, allocatable :: domlen(:), bufbeg(:)
      integer, allocatable :: pbig_send(:), pbig_recep(:)
      integer, allocatable :: pneig_send(:), pneig_recep(:)
      integer, allocatable :: reqs(:), reqr(:)
      real(t_p), allocatable :: work_tp(:)
      real(r_p), allocatable :: work_rp(:)
#if defined(USE_DETERMINISTIC_REDUCTION)
      mdyn_rtyp, allocatable :: work_md(:)
#endif
      logical :: layout_ready = .false.
   end type domain_t

   type, public :: polar_ticket_t
      private
      ! Noncopyable owner of in-flight MPI requests; polar_finish releases it.
      integer, allocatable :: receive_requests(:), send_requests(:)
      logical :: active = .false.
   end type polar_ticket_t

   type, public :: host_t
      private
      integer :: communicator = MPI_COMM_NULL
      integer :: rank = 0
   end type host_t

   type, public :: replica_t
      private
      integer :: communicator = MPI_COMM_NULL
      integer :: size = 0
      integer :: rank = 0
   end type replica_t

   type, public :: root_exch_t
      private
      integer :: communicator = MPI_COMM_NULL
   end type root_exch_t

   type, public :: parallel_t
      private
      ! Noncopyable MPI topology owner. Private components are loaned only through
      ! the handle queries and legacy adapter; reset_topology releases every loan.
      type(domain_t) :: domain
      type(host_t) :: host
      type(replica_t) :: replica
      type(root_exch_t) :: root_exchange
   end type parallel_t

   public :: construct_topology, configure_layout, reset_topology, release
   public :: domain_handle, domain_size, host_handle, replica_handle, replica_size, replica_rank, &
      &root_exchange_handle
   integer, parameter, public :: topology_classical = 1, topology_polymer = 2
   integer, parameter, public :: halo_large = 1, halo_neighbour = 2
   integer, parameter, public :: reduction_root = 1, reduction_all = 2
   public :: extend, accumulate, migrate_atoms, polar_post_comm, polar_finish, &
      &replicate, reduce_sum
   public :: assignment(=)

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
      module subroutine construct_topology(platform, configuration, topology_kind, status_value)
         type(parallel_t), target, intent(inout) :: platform
         type(config_t), intent(in) :: configuration
         integer, intent(in) :: topology_kind
         type(status_t), intent(out) :: status_value
      end subroutine construct_topology

      module subroutine configure_layout(platform, configuration, structure_value, status_value)
         type(parallel_t), target, intent(inout) :: platform
         type(config_t), intent(in) :: configuration
         type(structure_t), intent(in) :: structure_value
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
         type(polar_ticket_t), intent(out) :: lhs
         type(polar_ticket_t), intent(in) :: rhs
      end subroutine reject_polar_ticket_copy

      module subroutine capture_halo_layout(domain)
         type(domain_t), intent(inout) :: domain
      end subroutine capture_halo_layout

      module subroutine refresh_domain_layout(domain)
         type(domain_t), intent(inout) :: domain
      end subroutine refresh_domain_layout

      module subroutine release_halo_layout(domain)
         type(domain_t), intent(inout) :: domain
      end subroutine release_halo_layout

      module subroutine activate_domain(platform)
         type(parallel_t), target, intent(inout) :: platform
      end subroutine activate_domain

      module subroutine deactivate_domain(platform)
         type(parallel_t), target, intent(inout) :: platform
      end subroutine deactivate_domain

      module subroutine migrate_atoms_legacy()
      end subroutine migrate_atoms_legacy

      module subroutine polar_post_comm_legacy(ticket, nrhs, values)
         type(polar_ticket_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(:,:,:)
      end subroutine polar_post_comm_legacy

      module subroutine polar_finish_legacy(ticket, nrhs, values)
         type(polar_ticket_t), intent(inout) :: ticket
         integer, intent(in) :: nrhs
         real(t_p), intent(inout) :: values(:,:,:)
      end subroutine polar_finish_legacy

      module subroutine reset_polar_ticket(ticket)
         type(polar_ticket_t), intent(inout) :: ticket
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
      module subroutine extend_legacy_tp(values, components, halo)
         real(t_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
      end subroutine extend_legacy_tp

      module subroutine accumulate_legacy_tp(values, components, halo)
         real(t_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
      end subroutine accumulate_legacy_tp

      module subroutine accumulate_rec_to_dir_legacy_tp(direct_values, reciprocal_values, components)
         real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
         integer, intent(in) :: components
      end subroutine accumulate_rec_to_dir_legacy_tp

#if TINKER_MIXED_PREC
      module subroutine extend_legacy_rp(values, components, halo)
         real(r_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
      end subroutine extend_legacy_rp

      module subroutine accumulate_legacy_rp(values, components, halo)
         real(r_p), intent(inout) :: values(*)
         integer, intent(in) :: components, halo
      end subroutine accumulate_legacy_rp

#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
      module subroutine extend_legacy_md(values, components, halo)
         mdyn_rtyp, intent(inout) :: values(*)
         integer, intent(in) :: components, halo
      end subroutine extend_legacy_md

      module subroutine accumulate_legacy_md(values, components, halo)
         mdyn_rtyp, intent(inout) :: values(*)
         integer, intent(in) :: components, halo
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

   subroutine extend_tp(values, components, halo)
      real(t_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo

      call extend_legacy_tp(values, components, halo)
   end subroutine extend_tp

#if TINKER_MIXED_PREC
   subroutine extend_rp(values, components, halo)
      real(r_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo

      call extend_legacy_rp(values, components, halo)
   end subroutine extend_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
   subroutine extend_md(values, components, halo)
      mdyn_rtyp, intent(inout) :: values(*)
      integer, intent(in) :: components, halo

      call extend_legacy_md(values, components, halo)
   end subroutine extend_md
#endif

   subroutine accumulate_tp(values, components, halo)
      real(t_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo

      call accumulate_legacy_tp(values, components, halo)
   end subroutine accumulate_tp

#if TINKER_MIXED_PREC
   subroutine accumulate_rp(values, components, halo)
      real(r_p), intent(inout) :: values(*)
      integer, intent(in) :: components, halo

      call accumulate_legacy_rp(values, components, halo)
   end subroutine accumulate_rp
#endif

#if defined(USE_DETERMINISTIC_REDUCTION)
   subroutine accumulate_md(values, components, halo)
      mdyn_rtyp, intent(inout) :: values(*)
      integer, intent(in) :: components, halo

      call accumulate_legacy_md(values, components, halo)
   end subroutine accumulate_md
#endif

   subroutine accumulate_rec_to_dir_tp(direct_values, reciprocal_values, components)
      real(t_p), intent(inout) :: direct_values(*), reciprocal_values(*)
      integer, intent(in) :: components

      call accumulate_rec_to_dir_legacy_tp(direct_values, reciprocal_values, components)
   end subroutine accumulate_rec_to_dir_tp

   subroutine migrate_atoms()

      call migrate_atoms_legacy()
   end subroutine migrate_atoms

   subroutine polar_post_comm(ticket, nrhs, values)
      type(polar_ticket_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(:,:,:)

      call polar_post_comm_legacy(ticket, nrhs, values)
   end subroutine polar_post_comm

   subroutine polar_finish(ticket, nrhs, values)
      type(polar_ticket_t), intent(inout) :: ticket
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) :: values(:,:,:)

      call polar_finish_legacy(ticket, nrhs, values)
   end subroutine polar_finish

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
