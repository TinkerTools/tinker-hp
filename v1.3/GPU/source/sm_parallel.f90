!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (parallel) parallel_lifecycle
   use beads     , only: contract, nbeads, nbeads_ctr
   use domdec
   use inform    , only: initDebugEnv
   use mpi
   use replicas  , only: COMM_ROOT2ROOT, nreps, rank_reploc, use_reps
   use status    , only: clear_status, failed, set_status, status_invalid_argument
   use tinmemory
   implicit none

contains

   module subroutine construct_topology(platform, configuration, topology_kind, status_value)
      implicit none
      type(parallel_t), target, intent(inout) :: platform
      type(config_t), intent(in) :: configuration
      integer, intent(in) :: topology_kind
      type(status_t), intent(out) :: status_value
      integer :: ierr, color, nbeads_para

      call clear_status(status_value)
      call reset_topology(platform)
      if (topology_kind /= topology_classical .and. topology_kind /= topology_polymer) then
         call set_status(status_value, status_invalid_argument, 'Unknown MPI topology kind')
         return
      end if
      call MPI_Comm_size(MPI_COMM_WORLD, nproctot, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, ranktot, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'Unable to inspect MPI world')
         return
      end if

      if (topology_kind == topology_polymer .and. configuration%replicas_enabled) then
         call set_status(status_value, status_invalid_argument, 'PIMD and classical replicas cannot be combined')
         return
      end if
      if (configuration%replicas_enabled) then
         if (configuration%replica_count <= 0 .or. nproctot < configuration%replica_count .or. &
               mod(nproctot, configuration%replica_count) /= 0) then
            call set_status(status_value, status_invalid_argument, 'Invalid classical replica MPI topology')
            return
         end if
         nproc = nproctot / configuration%replica_count
      else if (topology_kind == topology_polymer) then
         nbeads_para = merge(nbeads_ctr, nbeads, contract .and. nbeads_ctr > 1)
         if (nbeads_para <= 0 .or. (nproctot >= nbeads_para .and. mod(nproctot, nbeads_para) /= 0)) then
            call set_status(status_value, status_invalid_argument, 'Invalid PIMD MPI topology')
            return
         end if
         nproc = merge(1, nproctot / nbeads_para, nproctot < nbeads_para)
      else
         nproc = nproctot
      end if

      rank_reploc = ranktot / nproc
      call MPI_Comm_split(MPI_COMM_WORLD, rank_reploc, ranktot, platform%domain%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         call reset_topology(platform)
         call set_status(status_value, status_invalid_argument, 'Unable to create MPI domain topology')
         return
      end if
      call MPI_Comm_size(platform%domain%communicator, platform%domain%size, ierr)
      if (ierr /= MPI_SUCCESS) then
         call reset_topology(platform)
         call set_status(status_value, status_invalid_argument, 'Unable to inspect MPI domain topology')
         return
      end if
      call MPI_Comm_rank(platform%domain%communicator, platform%domain%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         call reset_topology(platform)
         call set_status(status_value, status_invalid_argument, 'Unable to inspect MPI domain topology')
         return
      end if
      call MPI_Comm_split_type(platform%domain%communicator, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &
         &platform%host%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         call reset_topology(platform)
         call set_status(status_value, status_invalid_argument, 'Unable to create MPI host topology')
         return
      end if
      call MPI_Comm_rank(platform%host%communicator, platform%host%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         call reset_topology(platform)
         call set_status(status_value, status_invalid_argument, 'Unable to inspect MPI host topology')
         return
      end if
      if (topology_kind == topology_polymer) then
         call MPI_Comm_split(MPI_COMM_WORLD, platform%domain%rank, ranktot, platform%replica%communicator, ierr)
         if (ierr /= MPI_SUCCESS) then
            call reset_topology(platform)
            call set_status(status_value, status_invalid_argument, 'Unable to create MPI polymer topology')
            return
         end if
         call MPI_Comm_size(platform%replica%communicator, platform%replica%size, ierr)
         call MPI_Comm_rank(platform%replica%communicator, platform%replica%rank, ierr)
      else if (configuration%replicas_enabled) then
         color = merge(0, MPI_UNDEFINED, platform%domain%rank == 0)
         call MPI_Comm_split(MPI_COMM_WORLD, color, ranktot, platform%root_exchange%communicator, ierr)
         if (ierr /= MPI_SUCCESS) then
            call reset_topology(platform)
            call set_status(status_value, status_invalid_argument, 'Unable to create MPI root exchange topology')
            return
         end if
      end if

      call publish_legacy_topology(platform, configuration)
      call initDebugEnv
      call activate_domain(platform)
      extra_alloc = platform%domain%size > 1
   end subroutine construct_topology

   module subroutine configure_layout(platform, configuration, structure_value, status_value)
      implicit none
      type(parallel_t), target, intent(inout) :: platform
      type(config_t), intent(in) :: configuration
      type(structure_t), intent(in) :: structure_value
      type(status_t), intent(out) :: status_value

      call clear_status(status_value)
      if (.not.configuration%key_discovered .or. structure_value%atom_count <= 0) then
         call set_status(status_value, status_invalid_argument, 'MPI layout requires complete configuration and geometry')
         return
      end if
      call drivermpi
      call reinitnl(0)
      call refresh_domain_layout(platform%domain)
      call activate_domain(platform)
      platform%domain%size = nproc
      platform%domain%rank = rank
   end subroutine configure_layout

   module subroutine reset_topology(platform)
      implicit none
      type(parallel_t), target, intent(inout) :: platform
      call release_topology_unchecked(platform)
   end subroutine reset_topology

   module subroutine release(platform, status_value)
      implicit none
      type(parallel_t), target, intent(inout) :: platform
      type(status_t), intent(inout) :: status_value
      integer :: ierr, local_failure, global_failure
      integer :: local_domain, minimum_domain, maximum_domain
      integer :: local_host, minimum_host, maximum_host
      integer :: local_replica, minimum_replica, maximum_replica
      logical :: mpi_is_initialized, mpi_is_finalized

      call MPI_Initialized(mpi_is_initialized, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'MPI platform release preflight failed')
         return
      end if
      if (.not. mpi_is_initialized) then
         if (has_live_topology(platform)) then
            call set_status(status_value, status_invalid_argument, &
               &'MPI platform release requires initialized MPI')
         else
            call clear_status(status_value)
         end if
         return
      end if
      call MPI_Finalized(mpi_is_finalized, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'MPI platform release preflight failed')
         return
      end if
      if (mpi_is_finalized) then
         if (has_live_topology(platform)) then
            call set_status(status_value, status_invalid_argument, &
               &'MPI platform release occurred after MPI_Finalize')
         else
            call clear_status(status_value)
         end if
         return
      end if

      local_failure = merge(1, 0, failed(status_value))
      if (active_polar_tickets() > 0) local_failure = 1
      if (has_inconsistent_topology(platform)) local_failure = 1
      call MPI_Allreduce(local_failure, global_failure, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'MPI platform release preflight failed')
         return
      end if
      if (global_failure /= 0) then
         if (.not.failed(status_value)) then
            call set_status(status_value, status_invalid_argument, &
               &'MPI platform release rejected a rank-local lifecycle failure')
         end if
         return
      end if

      local_domain = merge(1, 0, platform%domain%communicator /= MPI_COMM_NULL)
      call MPI_Allreduce(local_domain, minimum_domain, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent domain state')
         return
      end if
      call MPI_Allreduce(local_domain, maximum_domain, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent domain state')
         return
      end if
      if (minimum_domain /= maximum_domain) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent domain state')
         return
      end if

      local_host = merge(1, 0, platform%host%communicator /= MPI_COMM_NULL)
      call MPI_Allreduce(local_host, minimum_host, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent host state')
         return
      end if
      call MPI_Allreduce(local_host, maximum_host, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS .or. minimum_host /= maximum_host) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent host state')
         return
      end if

      local_replica = merge(1, 0, platform%replica%communicator /= MPI_COMM_NULL)
      call MPI_Allreduce(local_replica, minimum_replica, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent replica state')
         return
      end if
      call MPI_Allreduce(local_replica, maximum_replica, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS .or. minimum_replica /= maximum_replica) then
         call set_status(status_value, status_invalid_argument, &
            &'MPI platform release found inconsistent replica state')
         return
      end if

      call release_topology_unchecked(platform, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'MPI platform communicator release failed')
      else
         call clear_status(status_value)
      end if
   end subroutine release

   subroutine release_topology_unchecked(platform, cleanup_status)
      implicit none
      type(parallel_t), target, intent(inout) :: platform
      integer, intent(out), optional :: cleanup_status
      integer :: first_error

      first_error = MPI_SUCCESS
      call deactivate_domain(platform)
      call release_halo_layout(platform%domain)
      call free_communicator(platform%root_exchange%communicator, first_error)
      call free_communicator(platform%replica%communicator, first_error)
      call free_communicator(platform%host%communicator, first_error)
      call free_communicator(comm_dir, first_error)
      call free_communicator(comm_rec, first_error)
      call free_communicator(platform%domain%communicator, first_error)
      if (present(cleanup_status)) cleanup_status = first_error
      if (first_error /= MPI_SUCCESS) then
         COMM_TINKER = platform%domain%communicator
         COMM_POLYMER = platform%replica%communicator
         COMM_ROOT2ROOT = platform%root_exchange%communicator
         hostcomm = platform%host%communicator
         return
      end if
      call reset_domain_state(platform%domain)
      platform%host = host_t()
      platform%replica = replica_t()
      platform%root_exchange = root_exch_t()
      COMM_TINKER = MPI_COMM_NULL
      COMM_POLYMER = MPI_COMM_NULL
      COMM_ROOT2ROOT = MPI_COMM_NULL
      nproc = 0
      rank = 0
      nproc_polymer = 0
      rank_polymer = 0
      hostcomm = MPI_COMM_NULL
      hostrank = 0
      rank_reploc = 0
      use_reps = .false.
      nreps = 0
      extra_alloc = .false.
   end subroutine release_topology_unchecked

   subroutine free_communicator(communicator, first_error)
      implicit none
      integer, intent(inout) :: communicator
      integer, intent(inout) :: first_error
      integer :: ierr

      if (communicator /= MPI_COMM_NULL .and. communicator /= MPI_COMM_WORLD) then
         call MPI_Comm_free(communicator, ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
   end subroutine free_communicator

   logical function has_live_topology(platform)
      implicit none
      type(parallel_t), intent(in) :: platform

      has_live_topology = platform%domain%communicator /= MPI_COMM_NULL .or. &
         &platform%host%communicator /= MPI_COMM_NULL .or. &
         &platform%replica%communicator /= MPI_COMM_NULL .or. &
         &platform%root_exchange%communicator /= MPI_COMM_NULL .or. &
         &comm_dir /= MPI_COMM_NULL .or. comm_rec /= MPI_COMM_NULL .or. &
         &COMM_TINKER /= MPI_COMM_NULL .or. COMM_POLYMER /= MPI_COMM_NULL .or. &
         &COMM_ROOT2ROOT /= MPI_COMM_NULL .or. hostcomm /= MPI_COMM_NULL
   end function has_live_topology

   logical function has_inconsistent_topology(platform)
      implicit none
      type(parallel_t), intent(in) :: platform

      has_inconsistent_topology = .false.
      if (platform%domain%communicator == MPI_COMM_NULL .and. &
            &platform%host%communicator /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%domain%communicator /= MPI_COMM_NULL .and. &
            &platform%host%communicator == MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%domain%communicator == MPI_COMM_NULL .and. &
            &COMM_TINKER /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%domain%communicator /= MPI_COMM_NULL .and. &
            &COMM_TINKER == MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%host%communicator == MPI_COMM_NULL .and. &
            &hostcomm /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%host%communicator /= MPI_COMM_NULL .and. &
            &hostcomm == MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%replica%communicator == MPI_COMM_NULL .and. &
            &COMM_POLYMER /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%replica%communicator /= MPI_COMM_NULL .and. &
            &COMM_POLYMER == MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%root_exchange%communicator == MPI_COMM_NULL .and. &
            &COMM_ROOT2ROOT /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%root_exchange%communicator /= MPI_COMM_NULL .and. &
            &COMM_ROOT2ROOT == MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (platform%domain%communicator == MPI_COMM_NULL .and. &
            &(comm_dir /= MPI_COMM_NULL .or. comm_rec /= MPI_COMM_NULL)) has_inconsistent_topology = .true.
      if (comm_dir /= MPI_COMM_NULL .and. comm_rec /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (comm_dir /= MPI_COMM_NULL .or. comm_rec /= MPI_COMM_NULL) then
         if (rank < ndir .and. comm_dir == MPI_COMM_NULL) has_inconsistent_topology = .true.
         if (rank >= ndir .and. comm_dir /= MPI_COMM_NULL) has_inconsistent_topology = .true.
         if (rank < ndir .and. comm_rec /= MPI_COMM_NULL) has_inconsistent_topology = .true.
         if (rank >= ndir .and. comm_rec == MPI_COMM_NULL) has_inconsistent_topology = .true.
      end if
      if (use_reps .and. platform%domain%rank == 0 .and. &
            &platform%root_exchange%communicator == MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (use_reps .and. platform%domain%rank /= 0 .and. &
            &platform%root_exchange%communicator /= MPI_COMM_NULL) has_inconsistent_topology = .true.
      if (.not.use_reps .and. platform%root_exchange%communicator /= MPI_COMM_NULL) &
         &has_inconsistent_topology = .true.
   end function has_inconsistent_topology

   subroutine reset_domain_state(domain)
      implicit none
      type(domain_t), intent(inout) :: domain

      domain%communicator = MPI_COMM_NULL
      domain%size = 0
      domain%rank = 0
      domain%nproc = 0
      domain%nloc = 0
      domain%nbig_send = 0
      domain%nbig_recep = 0
      domain%nneig_send = 0
      domain%nneig_recep = 0
      domain%layout_ready = .false.
   end subroutine reset_domain_state

   module subroutine reject_parallel_copy(lhs, rhs)
      implicit none
      type(parallel_t), intent(out) :: lhs
      type(parallel_t), intent(in) :: rhs

      error stop 'parallel_t is noncopyable'
   end subroutine reject_parallel_copy

   module subroutine reject_domain_copy(lhs, rhs)
      implicit none
      type(domain_t), intent(out) :: lhs
      type(domain_t), intent(in) :: rhs

      error stop 'domain_t is noncopyable'
   end subroutine reject_domain_copy

   module subroutine reject_polar_ticket_copy(lhs, rhs)
      implicit none
      type(polar_ticket_t), intent(out) :: lhs
      type(polar_ticket_t), intent(in) :: rhs

      error stop 'polar_ticket_t is noncopyable'
   end subroutine reject_polar_ticket_copy

   subroutine publish_legacy_topology(platform, configuration)
      implicit none
      type(parallel_t), intent(in) :: platform
      type(config_t), intent(in) :: configuration

      COMM_TINKER = platform%domain%communicator
      nproc = platform%domain%size
      rank = platform%domain%rank
      hostcomm = platform%host%communicator
      hostrank = platform%host%rank
      COMM_POLYMER = platform%replica%communicator
      nproc_polymer = platform%replica%size
      rank_polymer = platform%replica%rank
      COMM_ROOT2ROOT = platform%root_exchange%communicator
      use_reps = configuration%replicas_enabled
      nreps = configuration%replica_count
      comm_dir = MPI_COMM_NULL
      comm_rec = MPI_COMM_NULL
   end subroutine publish_legacy_topology
end submodule parallel_lifecycle
