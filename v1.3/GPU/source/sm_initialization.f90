!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  submodule initialization_lifecycle  --  runtime setup      ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
submodule (initialization) initialization_lifecycle
   use atoms , only: n
   use beads , only: read_bead_config
   use config, only: config_t, discover_key, discover_prm, apply_early_configuration, &
                     apply_defaults_prm_key
   use cflux , only: release_cflux_workspaces
   use cutoff, only: use_ewald
   use energi, only: delete_energi_on_device
   use fft   , only: free_FFTgrid_p, fft_final, is_fftInit
   use inform, only: app_id, analyze_a, analyze_beads_a, bar_a, dynamic_a, minimize_a, &
                     pibar_a, pimd_a, radial_a, testgrad_a
   use interfaces, only: mechanic
   use mdevice   , only: channel_t, direct_work, initialize_device, &
                     recip_work => reciprocal_work, default_work, sync_work, finalize_device
   use moldyn, only: a, aalt, aalt2, v
   use mpi, only: MPI_COMM_NULL, MPI_Finalized, MPI_Initialized, MPI_SUCCESS
   use parallel, only: active_polar_tickets, construct_topology, configure_layout, release, host_handle, &
      &domain_handle, topology_classical, topology_polymer
   use pme, only: pme_device_ready, release_pme_device
   use polar, only: release_polar_workspace
   use potent, only: use_charge, use_mpole, use_vdw
   use reduction, only: delete_reduction_buffers
   use status, only: clear_status, enrich_status, failed, set_status, status_invalid_argument, status_ok
   use structure, only: load_initial
   use tmatxb, only: release_tmatxb_workspace
   use tinMemory, only: prmem_requestm, sd_prmem
   use utils, only: set_to_zero2m
   use virial, only: delete_vir_device
#ifdef _OPENACC
   use thrust, only: thrust_clear_cache_memory
#endif
#ifdef USE_NVSHMEM
   use nvshmem, only: nvshmem_finalize
#endif
   implicit none
contains

   module subroutine initialize_defaults(raw_cli, status_value)
      use atoms
      use beads
      use bath
      use bound
      use config, only: publish_legacy_cli
      use deriv
      use domdec
      use energi
      use files
      use inform
      use iounit
      use keys
      use linmin
      use minima
      use moldyn
      use mpi
      use neigh
      use output
      use params
      use precis, only: p_tiny => tiny, small, p_huge => huge
      use timestat
      use tinMemory
      use tinheader
      use uprior
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      type(status_t), intent(out) :: status_value
      integer :: ierr
      logical :: mpi_ready, mpi_is_finalized
      real(t_p) :: precise

      call clear_status(status_value)
      call MPI_Initialized(mpi_ready, ierr)
      if (ierr /= MPI_SUCCESS .or. .not.mpi_ready) then
         call set_status(status_value, status_invalid_argument, &
            'Default initialization requires MPI_Init')
         return
      end if
      call MPI_Finalized(mpi_is_finalized, ierr)
      if (ierr /= MPI_SUCCESS .or. mpi_is_finalized) then
         call set_status(status_value, status_invalid_argument, &
            'Default initialization requires active MPI')
         return
      end if

      nbeads = 1
      nbeads_ctr = 0
      nproc = 1
      nthread = 1
      nrec = 0
      input = 5
      iout = 6
      call publish_legacy_cli(raw_cli)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproctot,ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'Unable to inspect MPI world size')
         return
      end if
      call MPI_COMM_RANK(MPI_COMM_WORLD,ranktot,ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'Unable to inspect MPI world rank')
         return
      end if
      if (ranktot.eq.0) call promo
      p_tiny = precise(1)
      small = precise(2)
      p_huge = precise(3)
      nkey = 0
      nprm = 0
      n = 0
      nprior = 0
      silent = .false.
      verbose = .false.
      debug = .false.
      abort = .false.
      use_bounds = .false.
      use_polymer = .false.
      etot_ave = 0.0
      etot_std = 0.0
      isothermal = .false.
      isobaric = .false.
      auto_lst = .true.
      dovlst = .true.
      domlst = .true.
      coordtype = 'NONE'
      call initatom
      polpred = 'ASPC'
      use_pred = (app_id == dynamic_a)
      call initres
      call init_memhandle
      nalt = i_init
      nalt2 = i_init
      stepfast = i_init
      stepint = i_init
      fctmin = 0.0_ti_p
      maxiter = 0
      nextiter = 0
      iprint = -1
      iwrite = -1
      n_fwriten = 0
      stpmax = 0.0_ti_p
      call initiate_timers()
      call timer_enter(timer_prog)
      dotstgrad = .false.
      calc_e = .true.
      call set_status(status_value, status_ok, 'Initialization defaults established')
   end subroutine initialize_defaults

   module subroutine initialize(application, raw_cli, configuration, molecular_structure, &
      &platform, reciprocal_work, status_value, probe, mode)
      implicit none
      type(application_t), intent(in) :: application
      type(raw_cli_t), intent(in) :: raw_cli
      type(config_t), intent(out) :: configuration
      type(structure_t), intent(out) :: molecular_structure
      type(parallel_t), intent(inout) :: platform
      type(channel_t), intent(out) :: reciprocal_work
      type(status_t), intent(out) :: status_value
      type(init_probe_t), intent(inout), optional :: probe
      integer, intent(in), optional :: mode
      integer :: ierr
      integer :: topology_kind
      logical :: mpi_ready, mpi_finalized, platform_ready, device_ready
      logical :: deferred
      external :: initDevice

      call clear_status(status_value)
      reciprocal_work = channel_t()
      platform_ready = .false.
      device_ready = .false.
      if (present(probe)) then
         probe%count = 0
         probe%phase = 0
      end if
      call MPI_Initialized(mpi_ready, ierr)
      if (ierr /= MPI_SUCCESS .or. .not.mpi_ready) then
         call set_status(status_value, status_invalid_argument, &
            'Common initialization requires MPI_Init')
         return
      end if
      call MPI_Finalized(mpi_finalized, ierr)
      if (ierr /= MPI_SUCCESS .or. mpi_finalized) then
         call set_status(status_value, status_invalid_argument, &
            'Common initialization requires active MPI')
         return
      end if
      select case (application%identity)
      case (dynamic_application%identity)
         topology_kind = topology_classical
         app_id = dynamic_a
      case (minimize_application%identity)
         topology_kind = topology_classical
         app_id = minimize_a
      case (radial_application%identity)
         topology_kind = topology_classical
         app_id = radial_a
      case (analyze_application%identity)
         topology_kind = topology_classical
         app_id = analyze_a
      case (pimd_application%identity)
         topology_kind = topology_polymer
         app_id = pimd_a
      case (analyze_beads_application%identity)
         topology_kind = topology_polymer
         app_id = analyze_beads_a
      case (bar_application%identity)
         topology_kind = topology_classical
         app_id = bar_a
      case (pibar_application%identity)
         topology_kind = topology_polymer
         app_id = pibar_a
      case (testgrad_application%identity)
         topology_kind = topology_classical
         app_id = testgrad_a
      case default
         call set_status(status_value, status_invalid_argument, 'Unknown application identity')
         return
      end select

      deferred = application%identity == bar_application%identity .or. &
         &application%identity == pibar_application%identity
      if (deferred) then
         if (.not.present(mode) .or. (mode /= energy_generation .and. mode /= bar_file_only .and. &
            &mode /= platform_bootstrap)) then
            call set_status(status_value, status_invalid_argument, 'BAR/PIBAR require a valid mode')
            return
         end if
      else if (present(mode)) then
         call set_status(status_value, status_invalid_argument, &
            'Only BAR and PIBAR accept an initialization mode')
         return
      end if

      if (deferred .and. present(mode)) then
         if (mode == energy_generation .and. domain_handle(platform) /= MPI_COMM_NULL) then
            call initialize_energy_state(application, raw_cli, configuration, molecular_structure, platform, &
               &status_value, 1)
            return
         end if
      end if

      call run_phase(phase_identity, probe, status_value)
      if (failed(status_value)) return
      call run_phase(phase_legacy_initial, probe, status_value)
      if (failed(status_value)) return
      call initialize_defaults(raw_cli, status_value)
      if (failed(status_value)) return
      if (deferred) then
         configuration = config_t()
         if (application%identity == pibar_application%identity) then
            call discover_key(raw_cli, configuration, status_value)
            if (failed(status_value)) return
            call apply_early_configuration(configuration)
            call read_bead_config()
         end if
         call construct_topology(platform, configuration, topology_kind, status_value)
         if (failed(status_value)) return
         platform_ready = .true.
         call initialize_device(host_handle(platform), status_value)
         if (failed(status_value)) then
            call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
            return
         end if
         device_ready = .true.
         reciprocal_work = recip_work
         call initDevice(.true.)
         if (mode == bar_file_only .or. mode == platform_bootstrap) return
         call initialize_energy_state(application, raw_cli, configuration, molecular_structure, platform, &
            &status_value, 1)
         if (failed(status_value)) then
            call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         end if
         return
      end if
      call run_phase(phase_key_discovery, probe, status_value)
      if (failed(status_value)) return
      call discover_key(raw_cli, configuration, status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic Configuration')
         else
            call enrich_status(status_value, 'Configuration')
         end if
         return
      end if
      call run_phase(phase_prm_discovery, probe, status_value)
      if (failed(status_value)) return
      call discover_prm(raw_cli, configuration, status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic Configuration')
         else
            call enrich_status(status_value, 'Configuration')
         end if
         return
      end if
      call run_phase(phase_topology_configuration, probe, status_value)
      if (failed(status_value)) return
      call apply_early_configuration(configuration)
      if (topology_kind == topology_polymer) call read_bead_config()
      call run_phase(phase_mpi_platform, probe, status_value)
      if (failed(status_value)) return
      call construct_topology(platform, configuration, topology_kind, status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic MPI Platform')
         else
            call enrich_status(status_value, 'MPI Platform')
         end if
         return
      end if
      platform_ready = .true.
      call run_phase(phase_gpu_platform, probe, status_value)
      if (failed(status_value)) then
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call initialize_device(host_handle(platform), status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic GPU Platform')
         else
            call enrich_status(status_value, 'GPU Platform')
         end if
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      device_ready = .true.
      reciprocal_work = recip_work
      call initDevice(.true.)
      call run_phase(phase_structure, probe, status_value)
      if (failed(status_value)) then
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call load_initial(configuration, molecular_structure, status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic Structure')
         else
            call enrich_status(status_value, 'Structure')
         end if
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call run_phase(phase_complete_configuration, probe, status_value)
      if (failed(status_value)) then
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call apply_defaults_prm_key(configuration, status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic Configuration')
         else
            call enrich_status(status_value, 'Configuration')
         end if
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call run_phase(phase_geometry, probe, status_value)
      if (failed(status_value)) then
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      if (application%identity == radial_application%identity) then
         use_vdw = .true.
         use_charge = .false.
         use_mpole = .false.
         use_ewald = .false.
      end if
      call unitcell
      call cutoffs
      call lattice
      call run_phase(phase_layout_neighbours, probe, status_value)
      if (failed(status_value)) then
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call configure_layout(platform, configuration, molecular_structure, status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic MPI Platform')
         else
            call enrich_status(status_value, 'MPI Platform')
         end if
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      if (application%identity == dynamic_application%identity .or. &
         &application%identity == pimd_application%identity) then
         call prmem_requestm(v,3,n)
         call prmem_requestm(a,3,n)
         call prmem_requestm(aalt,3,n)
         call prmem_requestm(aalt2,3,n)
         call set_to_zero2m(a,v,3*n,reciprocal_work%queue)
         call set_to_zero2m(aalt,aalt2,3*n,reciprocal_work%queue)
      end if
      call run_phase(phase_mechanics, probe, status_value)
      if (failed(status_value)) then
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
         return
      end if
      call mechanic(configuration,status_value)
      if (failed(status_value)) then
         if (application%identity == dynamic_application%identity) then
            call enrich_status(status_value, 'Dynamic Configuration')
         else
            call enrich_status(status_value, 'Configuration')
         end if
         call cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
      end if
   end subroutine initialize

   module subroutine initialize_state(application, raw_cli, configuration, molecular_structure, &
      &platform, status_value, state)
      implicit none
      type(application_t), intent(in) :: application
      type(raw_cli_t), intent(in) :: raw_cli
      type(config_t), intent(inout) :: configuration
      type(structure_t), intent(out) :: molecular_structure
      type(parallel_t), intent(inout) :: platform
      type(status_t), intent(out) :: status_value
      integer, intent(in) :: state

      call clear_status(status_value)
      if (application%identity /= bar_application%identity .and. &
         &application%identity /= pibar_application%identity) then
         call set_status(status_value, status_invalid_argument, &
            'Only BAR and PIBAR expose repeated energy-state initialization')
         return
      end if
      if (state < 2 .or. domain_handle(platform) == MPI_COMM_NULL .or. &
         &.not.configuration%complete) then
         call set_status(status_value, status_invalid_argument, &
            'Repeated energy state requires a completed BAR/PIBAR bootstrap')
         return
      end if
      call initialize_energy_state(application, raw_cli, configuration, molecular_structure, platform, &
         &status_value, state)
   end subroutine initialize_state

   subroutine initialize_energy_state(application, raw_cli, configuration, molecular_structure, platform, &
      &status_value, state)
      implicit none
      type(application_t), intent(in) :: application
      type(raw_cli_t), intent(in) :: raw_cli
      type(config_t), intent(inout) :: configuration
      type(structure_t), intent(out) :: molecular_structure
      type(parallel_t), intent(inout) :: platform
      type(status_t), intent(out) :: status_value
      integer, intent(in) :: state

      call clear_status(status_value)
      if (state < 1) then
         call set_status(status_value, status_invalid_argument, 'Energy-state index must be positive')
         return
      end if
      configuration = config_t()
      call discover_key(raw_cli, configuration, status_value)
      if (failed(status_value)) return
      call discover_prm(raw_cli, configuration, status_value)
      if (failed(status_value)) return
      call apply_early_configuration(configuration)
      if (application%identity == pibar_application%identity) call read_bead_config()
      call load_initial(configuration, molecular_structure, status_value)
      if (failed(status_value)) return
      call apply_defaults_prm_key(configuration, status_value)
      if (failed(status_value)) return
      call unitcell
      call cutoffs
      call lattice
      call configure_layout(platform, configuration, molecular_structure, status_value)
      if (failed(status_value)) return
      call mechanic(configuration, status_value)
      if (failed(status_value)) return
   end subroutine initialize_energy_state

   subroutine run_phase(phase, probe, status_value)
      implicit none
      integer, intent(in) :: phase
      type(init_probe_t), intent(inout), optional :: probe
      type(status_t), intent(out) :: status_value

      call clear_status(status_value)
      if (.not.present(probe)) return
      probe%count = probe%count + 1
      probe%phase(probe%count) = phase
      if (probe%fail_at == phase) then
         call set_status(status_value, status_invalid_argument, 'Injected initialization phase failure')
      end if
   end subroutine run_phase

   subroutine cleanup(application, platform, reciprocal_work, platform_ready, device_ready, status_value)
      implicit none
      type(application_t), intent(in) :: application
      type(parallel_t), intent(inout) :: platform
      type(channel_t), intent(inout) :: reciprocal_work
      logical, intent(in) :: platform_ready, device_ready
      type(status_t), intent(inout) :: status_value
      type(status_t) :: cleanup_status
      type(status_t) :: initial_status
      type(status_t) :: dev_clean_stat
      integer :: pme_cleanup_code

      call clear_status(cleanup_status)
      if (device_ready) then
         if (is_fftInit) then
            call fft_final
            call free_FFTgrid_p
            is_fftInit = .false.
         end if
         call delete_reduction_buffers(direct_work%queue)
         call release_cflux_workspaces(direct_work%queue)
         call release_polar_workspace(direct_work%queue)
         call release_tmatxb_workspace(direct_work%queue)
!$acc wait
#ifdef _OPENACC
         call thrust_clear_cache_memory(sd_prmem)
#endif
         pme_cleanup_code = 0
         call release_pme_device(reciprocal_work%queue, pme_cleanup_code)
         if (pme_cleanup_code /= MPI_SUCCESS) then
            call set_status(cleanup_status, status_invalid_argument, &
               &'Initialization cleanup could not drain PME MPI requests')
         else if (pme_device_ready) then
            call set_status(cleanup_status, status_invalid_argument, &
               &'Initialization cleanup retained PME device state')
         end if
         call delete_energi_on_device()
         call delete_vir_device()
#ifdef USE_NVSHMEM
         call nvshmem_finalize()
#endif
         call finalize_device(dev_clean_stat)
         if (failed(dev_clean_stat)) cleanup_status = dev_clean_stat
         if (direct_work%queue >= 0 .or. recip_work%queue >= 0 .or. &
               &default_work%queue >= 0 .or. sync_work%queue >= 0 .or. &
               &direct_work%stream /= 0 .or. recip_work%stream /= 0 .or. &
               &default_work%stream /= 0 .or. sync_work%stream /= 0) then
            call set_status(cleanup_status, status_invalid_argument, &
               &'Initialization cleanup retained a GPU execution channel')
         end if
         if (active_polar_tickets() /= 0) then
            call set_status(cleanup_status, status_invalid_argument, &
               &'Initialization cleanup retained an opaque MPI ticket')
         end if
         if (failed(cleanup_status)) then
            if (application%identity == dynamic_application%identity) then
               call enrich_status(status_value, 'Dynamic cleanup')
            else
               call enrich_status(status_value, 'Initialization cleanup')
            end if
         end if
      end if
      if (failed(cleanup_status) .and. .not.failed(status_value)) status_value = cleanup_status
      if (platform_ready) then
         initial_status = status_value
         call release(platform, status_value)
         if (failed(status_value)) then
            call clear_status(status_value)
            call release(platform, status_value)
            if (.not.failed(status_value)) status_value = initial_status
         end if
      end if
      reciprocal_work = channel_t()
   end subroutine cleanup
end submodule initialization_lifecycle
