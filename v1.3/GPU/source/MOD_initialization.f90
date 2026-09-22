!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module initialization  --  common runtime lifecycle       ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module initialization
   use config, only: raw_cli_t, config_t
   use mdevice, only: channel_t
   use parallel, only: parallel_t
   use status, only: status_t
   use structure, only: structure_t
   implicit none

   type, public :: application_t
      integer :: identity = 0
   end type application_t

   type(application_t), parameter, public :: dynamic_application = application_t(1)
   type(application_t), parameter, public :: minimize_application = application_t(2)
   type(application_t), parameter, public :: radial_application = application_t(3)
   type(application_t), parameter, public :: analyze_application = application_t(4)
   type(application_t), parameter, public :: pimd_application = application_t(5)
   type(application_t), parameter, public :: analyze_beads_application = application_t(6)
   type(application_t), parameter, public :: bar_application = application_t(7)
   type(application_t), parameter, public :: pibar_application = application_t(8)
   type(application_t), parameter, public :: testgrad_application = application_t(9)

   integer, parameter, public :: energy_generation = 1
   integer, parameter, public :: bar_file_only = 2
   integer, parameter, public :: platform_bootstrap = 3

   integer, parameter, public :: initialization_phase_count = 12
   integer, parameter, public :: phase_identity = 1
   integer, parameter, public :: phase_legacy_initial = 2
   integer, parameter, public :: phase_key_discovery = 3
   integer, parameter, public :: phase_prm_discovery = 4
   integer, parameter, public :: phase_topology_configuration = 5
   integer, parameter, public :: phase_mpi_platform = 6
   integer, parameter, public :: phase_gpu_platform = 7
   integer, parameter, public :: phase_structure = 8
   integer, parameter, public :: phase_complete_configuration = 9
   integer, parameter, public :: phase_geometry = 10
   integer, parameter, public :: phase_layout_neighbours = 11
   integer, parameter, public :: phase_mechanics = 12

   ! Test-only observer; production callers do not provide it.
   type, public :: init_probe_t
      integer :: fail_at = 0
      integer :: count = 0
      integer :: phase(initialization_phase_count) = 0
   end type init_probe_t

   public :: initialize
   public :: initialize_defaults
   public :: initialize_state
   interface
      module subroutine initialize_defaults(raw_cli, status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         type(status_t), intent(out) :: status_value
      end subroutine initialize_defaults

      module subroutine initialize(application, raw_cli, configuration, molecular_structure, &
         &platform, reciprocal_work, status_value, probe, mode)
         type(application_t), intent(in) :: application
         type(raw_cli_t), intent(in) :: raw_cli
         type(config_t), intent(out) :: configuration
         type(structure_t), intent(out) :: molecular_structure
         type(parallel_t), intent(inout) :: platform
         type(channel_t), intent(out) :: reciprocal_work
         type(status_t), intent(out) :: status_value
         type(init_probe_t), intent(inout), optional :: probe
         integer, intent(in), optional :: mode
      end subroutine initialize

      module subroutine initialize_state(application, raw_cli, configuration, molecular_structure, &
         &platform, status_value, state)
         type(application_t), intent(in) :: application
         type(raw_cli_t), intent(in) :: raw_cli
         type(config_t), intent(inout) :: configuration
         type(structure_t), intent(out) :: molecular_structure
         type(parallel_t), intent(inout) :: platform
         type(status_t), intent(out) :: status_value
         integer, intent(in) :: state
      end subroutine initialize_state
   end interface
end module initialization
