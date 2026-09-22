!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module config  --  immutable command line parsing surface   ##
!     ##                                                              ##
!     ##################################################################
#include "tinker_macro.h"

module config
   use sizes, only: maxkey, maxprm
   use status, only: status_t
   implicit none

   integer, parameter :: max_raw_cli_args = 20
   integer, parameter :: raw_cli_length = 240
   integer, parameter :: key_from_none = 0
   integer, parameter :: key_from_cli = 1
   integer, parameter :: key_from_xyz = 2
   integer, parameter :: key_from_working_dir = 3

   type :: raw_cli_t
      integer :: count = 0
      character(len=raw_cli_length) :: argument(0:max_raw_cli_args) = ''
   end type raw_cli_t

   type :: config_t
      logical :: key_discovered = .false.
      logical :: key_found = .false.
      logical :: complete = .false.
      logical :: minimize_xyz_enabled = .true.
      logical :: dcd_enabled = .false.
      logical :: replicas_enabled = .false.
      logical :: pme_enabled = .false.
      integer :: pme_processes = 0
      integer :: replica_count = 0
      integer :: key_count = 0
      integer :: prm_count = 0
      integer :: key_origin = key_from_none
      integer :: xyz_base_length = 0
      integer :: xyz_directory_length = 0
      character(len=raw_cli_length) :: xyz_path = ''
      character(len=raw_cli_length) :: key_path = ''
      character(len=raw_cli_length) :: key_line(maxkey) = ''
      character(len=raw_cli_length) :: prm_path = ''
      character(len=raw_cli_length) :: prm_line(maxprm) = ''
   end type config_t

   interface

      module subroutine parse_positional_cli(raw_cli, values, value_count, status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         character(len=*), intent(out) :: values(:)
         integer, intent(out) :: value_count
         type(status_t), intent(out) :: status_value
      end subroutine parse_positional_cli

      module subroutine select_positional_cli(raw_cli, coordinate_index, key_index, selected_cli, &
         &status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         integer, intent(in) :: coordinate_index, key_index
         type(raw_cli_t), intent(out) :: selected_cli
         type(status_t), intent(out) :: status_value
      end subroutine select_positional_cli

      module subroutine publish_legacy_cli(raw_cli)
         type(raw_cli_t), intent(in) :: raw_cli
      end subroutine publish_legacy_cli

      module subroutine parse_common_cli(raw_cli, xyz_path, values, value_count, status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         character(*), intent(out) :: xyz_path
         character(*), intent(out) :: values(:)
         integer, intent(out) :: value_count
         type(status_t), intent(out) :: status_value
      end subroutine parse_common_cli

      module subroutine discover_key(raw_cli, configuration, status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         type(config_t), intent(inout) :: configuration
         type(status_t), intent(out) :: status_value
      end subroutine discover_key

      module subroutine apply_early_configuration(configuration)
         type(config_t), intent(in) :: configuration
      end subroutine apply_early_configuration

      module subroutine discover_prm(raw_cli, configuration, status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         type(config_t), intent(inout) :: configuration
         type(status_t), intent(out) :: status_value
      end subroutine discover_prm

      module subroutine apply_defaults_prm_key(configuration, status_value)
         type(config_t), intent(inout) :: configuration
         type(status_t), intent(out) :: status_value
      end subroutine apply_defaults_prm_key

      module subroutine set_potential_defaults()
      end subroutine set_potential_defaults

      module subroutine parse_dynamic_cli(raw_cli, periodic, nstep, timestep_fs, &
            dump_interval_ps, ensemble, temperature, pressure, status_value)
         type(raw_cli_t), intent(in) :: raw_cli
         logical, intent(in) :: periodic
         integer, intent(out) :: nstep
         real(r_p), intent(out) :: timestep_fs
         real(r_p), intent(out) :: dump_interval_ps
         integer, intent(out) :: ensemble
         real(r_p), intent(out) :: temperature
         real(r_p), intent(out) :: pressure
         type(status_t), intent(out) :: status_value
      end subroutine parse_dynamic_cli
   end interface
contains

   subroutine capture_raw_cli(raw_cli)
      implicit none
      type(raw_cli_t), intent(out) :: raw_cli
      integer :: i

      raw_cli%count = min(command_argument_count(), max_raw_cli_args)
      raw_cli%argument = ''
      do i = 0, raw_cli%count
         call get_command_argument(i, raw_cli%argument(i))
      end do
   end subroutine capture_raw_cli
end module config
