!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"

submodule (config) config_cli
   use ascii,     only: backslash
   use iounit,    only: iout
   use tinheader ,only: re_p
   use status    ,only: status_t, status_missing_argument, status_invalid_argument &
                 ,clear_status, set_status, failed
   implicit none
contains

   module subroutine parse_positional_cli(raw_cli, values, value_count, status_value)
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      character(len=*), intent(out) :: values(:)
      integer, intent(out) :: value_count
      type(status_t), intent(out) :: status_value
      integer :: i

      call clear_status(status_value)
      values = ''
      value_count = 0
      i = 1
      do while (i <= raw_cli%count)
         if (letter_option(raw_cli%argument(i))) then
            i = i + 2
         else
            if (value_count == size(values)) then
               call set_status(status_value, status_invalid_argument, &
                  'Too Many Positional Command Line Arguments')
               return
            end if
            value_count = value_count + 1
            values(value_count) = raw_cli%argument(i)
            i = i + 1
         end if
      end do
   end subroutine parse_positional_cli

   module subroutine select_positional_cli(raw_cli, coordinate_index, key_index, selected_cli, &
      &status_value)
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      integer, intent(in) :: coordinate_index, key_index
      type(raw_cli_t), intent(out) :: selected_cli
      type(status_t), intent(out) :: status_value
      character(len=raw_cli_length) :: values(max_raw_cli_args)
      character(len=raw_cli_length) :: key_path
      integer :: value_count, i
      logical :: exists

      call parse_positional_cli(raw_cli, values, value_count, status_value)
      if (failed(status_value)) return
      if (coordinate_index < 1 .or. coordinate_index > value_count .or. &
         &key_index < 1 .or. key_index > value_count) then
         call set_status(status_value, status_missing_argument, &
            'BAR/PIBAR requires both coordinate state arguments')
         return
      end if

      selected_cli = raw_cli_t()
      selected_cli%count = 1
      selected_cli%argument(0) = raw_cli%argument(0)
      selected_cli%argument(1) = values(coordinate_index)

      i = 1
      do while (i <= raw_cli%count)
         if (key_option(raw_cli%argument(i))) then
            if (i == raw_cli%count) then
               call set_status(status_value, status_missing_argument, &
                  'You Need To Enter the Keyfile Name')
               return
            end if
            selected_cli%count = selected_cli%count + 1
            selected_cli%argument(selected_cli%count) = raw_cli%argument(i)
            selected_cli%count = selected_cli%count + 1
            selected_cli%argument(selected_cli%count) = raw_cli%argument(i+1)
            return
         end if
         i = i + 1
      end do

      call state_key_path(values(key_index), key_path)
      inquire(file=trim(key_path), exist=exists)
      if (exists) then
         selected_cli%count = selected_cli%count + 1
         selected_cli%argument(selected_cli%count) = '-K'
         selected_cli%count = selected_cli%count + 1
         selected_cli%argument(selected_cli%count) = key_path
      end if
   end subroutine select_positional_cli

   subroutine state_key_path(xyz_path, key_path)
      implicit none
      character(*), intent(in) :: xyz_path
      character(*), intent(out) :: key_path
      integer :: base_length, directory_length
      logical :: exists

      call xyz_base(xyz_path, base_length, directory_length)
      key_path = xyz_path(1:base_length)//'.key'
      inquire(file=trim(key_path), exist=exists)
      if (.not.exists) then
         if (index(trim(xyz_path), '/') == 0) then
            key_path = 'tinker.key'
         else
            key_path = xyz_path(1:directory_length)//'tinker.key'
         end if
      end if
   end subroutine state_key_path

   module subroutine publish_legacy_cli(raw_cli)
      use argue, only: maxarg, narg, listarg, arg
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      integer :: i

      narg = min(raw_cli%count, maxarg)
      arg = ''
      do i = 0, narg
         arg(i) = raw_cli%argument(i)
      end do

      listarg = .false.
      do i = 1, narg
         listarg(i) = .true.
      end do
      do i = 1, narg
         if (letter_option(arg(i))) then
            listarg(i) = .false.
            if (i < narg) listarg(i+1) = .false.
         end if
      end do
   end subroutine publish_legacy_cli

   module subroutine parse_common_cli(raw_cli, xyz_path, values, value_count, status_value)
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      character(*), intent(out) :: xyz_path
      character(*), intent(out) :: values(:)
      integer, intent(out) :: value_count
      type(status_t), intent(out) :: status_value
      integer :: i

      call clear_status(status_value)
      xyz_path = ''
      values = ''
      value_count = 0
      i = 1
      do while (i <= raw_cli%count)
         if (letter_option(raw_cli%argument(i))) then
            i = i + 2
         else if (len_trim(xyz_path) == 0) then
            xyz_path = raw_cli%argument(i)
            i = i + 1
         else
            if (value_count == size(values)) then
               call set_status(status_value, status_invalid_argument, &
                  'Too Many Command Line Arguments')
               return
            end if
            value_count = value_count + 1
            values(value_count) = raw_cli%argument(i)
            i = i + 1
         end if
      end do

      if (len_trim(xyz_path) == 0) then
         call set_status(status_value, status_missing_argument, &
            'You Need To Enter the Coordinate File Name')
      end if
   end subroutine parse_common_cli

   module subroutine discover_key(raw_cli, configuration, status_value)
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      type(config_t), intent(inout) :: configuration
      type(status_t), intent(out) :: status_value
      character(len=raw_cli_length) :: values(max_raw_cli_args)
      character(len=raw_cli_length) :: explicit_key
      integer :: value_count, i, candidate_origin
      logical :: explicit_key_requested, exists

      call clear_status(status_value)
      if (configuration%key_discovered) return

      call parse_common_cli(raw_cli, configuration%xyz_path, values, value_count, status_value)
      if (failed(status_value)) return

      call xyz_base(configuration%xyz_path, configuration%xyz_base_length, &
         configuration%xyz_directory_length)
      explicit_key_requested = .false.
      explicit_key = ''
      do i = 1, raw_cli%count
         if (.not.key_option(raw_cli%argument(i))) cycle
         explicit_key_requested = .true.
         if (i == raw_cli%count) then
            call set_status(status_value, status_missing_argument, &
               'You Need To Enter the Keyfile Name')
            return
         end if
         explicit_key = raw_cli%argument(i+1)
      end do

      configuration%key_count = 0
      configuration%key_line = ''
      configuration%key_found = .false.
      configuration%key_origin = key_from_none
      configuration%key_path = ''
      if (explicit_key_requested) then
         configuration%key_path = explicit_key
         call suffix(configuration%key_path, 'key', 'old')
         inquire(file=configuration%key_path, exist=exists)
         if (.not.exists) then
            call set_status(status_value, status_invalid_argument, &
               'Keyfile Specified on Command Line was not Found')
            return
         end if
         candidate_origin = key_from_cli
      else
         configuration%key_path = configuration%xyz_path(1:configuration%xyz_base_length)//'.key'
         call version(configuration%key_path, 'old')
         inquire(file=configuration%key_path, exist=exists)
         candidate_origin = key_from_xyz
         if (.not.exists) then
            if (configuration%xyz_directory_length == 0) then
               configuration%key_path = 'tinker.key'
            else
               configuration%key_path = configuration%xyz_path(1: &
                  configuration%xyz_directory_length)//'tinker.key'
            end if
            call version(configuration%key_path, 'old')
            inquire(file=configuration%key_path, exist=exists)
            candidate_origin = key_from_working_dir
         end if
      end if

      configuration%key_found = exists
      if (exists) configuration%key_origin = candidate_origin
      if (exists) then
         call read_key_file(configuration, status_value)
         if (failed(status_value)) return
         call echo_key_comments(configuration)
         call discover_early_key_values(configuration, status_value)
         if (failed(status_value)) return
      end if
      configuration%key_discovered = .true.
   end subroutine discover_key

   module subroutine apply_early_configuration(configuration)
      use dcdmod, only: dcdio
      use domdec, only: nrec
      use files, only: filename, leng, ldir
      use keys, only: nkey, keyline
      use potent, only: use_pmecore
      use replicas, only: nreps, use_reps
      implicit none
      type(config_t), intent(in) :: configuration

      filename = configuration%xyz_path
      leng = configuration%xyz_base_length
      ldir = configuration%xyz_directory_length
      nkey = configuration%key_count
      keyline = ''
      if (nkey > 0) keyline(1:nkey) = configuration%key_line(1:nkey)
      dcdio = configuration%dcd_enabled
      use_reps = configuration%replicas_enabled
      nreps = configuration%replica_count
      use_pmecore = configuration%pme_enabled
      nrec = configuration%pme_processes
      call control
   end subroutine apply_early_configuration

   module subroutine discover_prm(raw_cli, configuration, status_value)
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      type(config_t), intent(inout) :: configuration
      type(status_t), intent(out) :: status_value
      character(len=raw_cli_length) :: values(max_raw_cli_args), prm_path, paramdir
      character(len=20) :: keyword
      character(len=raw_cli_length) :: record, string
      integer :: value_count, i, next, length
      logical :: exists

      call clear_status(status_value)
      configuration%prm_count = 0
      configuration%prm_path = ''
      configuration%prm_line = ''
      if (len_trim(configuration%xyz_path) == 0) then
         call set_status(status_value, status_missing_argument, &
            'You Need To Enter the Coordinate File Name')
         return
      end if
      prm_path = configuration%xyz_path(1:configuration%xyz_base_length)//'.prm'
      do i = 1, configuration%key_count
         record = configuration%key_line(i)
         next = 1
         call gettext(record, keyword, next)
         call upcase(keyword)
         if (trim(keyword) /= 'PARAMETERS') cycle
         string = record(next:raw_cli_length)
         next = 1
         call getstring(string, prm_path, next)
         if (next == 1) call gettext(string, prm_path, next)
      end do
      call get_environment_variable('TINKER_PARAMS', paramdir, length)
      if (length /= 0) prm_path = trim(paramdir)//'/'//trim(prm_path)
      call suffix(prm_path, 'prm', 'old')
      inquire(file=prm_path, exist=exists)
      if (.not.exists) then
         call parse_common_cli(raw_cli, record, values, value_count, status_value)
         if (failed(status_value)) return
         if (value_count > 0) then
            prm_path = values(1)
            call suffix(prm_path, 'prm', 'old')
            inquire(file=prm_path, exist=exists)
         end if
      end if
      if (.not.exists) then
         call set_status(status_value, status_invalid_argument, &
            'You need to Specify the parameter file in the keyfile or in the command line !')
         return
      end if
      configuration%prm_path = prm_path
      configuration%prm_count = 0
      configuration%prm_line = ''
      call read_prm_file(configuration, status_value)
   end subroutine discover_prm

   module subroutine apply_defaults_prm_key(configuration, status_value)
      use params, only: nprm, prmline
      use potent, only: use_mpole, use_polar
      implicit none
      type(config_t), intent(inout) :: configuration
      type(status_t), intent(out) :: status_value
      integer :: i

      call clear_status(status_value)
      if (len_trim(configuration%prm_path) == 0) then
         call set_status(status_value, status_invalid_argument, 'Parameter File was not Discovered')
         return
      end if
      call validate_ndir_keywords(configuration, status_value)
      if (failed(status_value)) return
      call set_potential_defaults()
      call initprm
      nprm = configuration%prm_count
      prmline = ''
      if (nprm > 0) prmline(1:nprm) = configuration%prm_line(1:nprm)
      call readprm
      do i = 1, configuration%key_count
         call prmkey(configuration%key_line(i))
      end do
!$acc update device(use_mpole,use_polar)
      configuration%complete = .true.
   end subroutine apply_defaults_prm_key

   module subroutine set_potential_defaults()
      use potent
      implicit none

      use_bond       = .true.
      use_angle      = .true.
      use_strbnd     = .true.
      use_urey       = .true.
      use_angang     = .true.
      use_opbend     = .true.
      use_opdist     = .true.
      use_improp     = .true.
      use_imptor     = .true.
      use_tors       = .true.
      use_pitors     = .true.
      use_angtor     = .true.
      use_strtor     = .true.
      use_tortor     = .true.
      use_vdw        = .true.
      use_vdwshort   = .false.
      use_vdwlong    = .false.
      use_charge     = .true.
      use_cshortreal = .false.
      use_clong      = .false.
      use_creal      = .true.
      use_crec       = .true.
      use_cself      = .true.
      use_mpole      = .true.
      use_mpoleshortreal = .false.
      use_mpolelong      = .false.
      use_mreal      = .true.
      use_mrec       = .true.
      use_mself      = .true.
      use_polar      = .true.
      use_preal      = .true.
      use_prec       = .true.
      use_pself      = .true.
      use_polarshortreal = .false.
      use_solv       = .false.
      use_geom       = .true.
      use_extra      = .false.
      use_emtp       = .false.
      use_repulsshort    = .false.
      use_repulslong     = .false.
      use_repuls         = .true.
      use_disp           = .true.
      use_dispshortreal  = .false.
      use_displong   = .false.
      use_dispreal   = .true.
      use_disprec    = .true.
      use_dispself   = .true.
      use_chgtrn      = .true.
      use_chgflx      = .true.
      use_lambdadyn   = .false.
      use_OSRW        = .false.
      fuse_chglj      = .false.
      fuse_bonded     = .false.
      use_mlpot       = .false.

      use_mlpot_only    = .false.
      use_ml_embedding  = .false.
      use_embd_potoff   = .false.
      use_embd_bond     = .true.
      use_embd_angle    = .true.
      use_embd_strbnd   = .true.
      use_embd_urey     = .true.
      use_embd_angang   = .true.
      use_embd_opbend   = .true.
      use_embd_opdist   = .true.
      use_embd_improp   = .true.
      use_embd_imptor   = .true.
      use_embd_tors     = .true.
      use_embd_pitors   = .true.
      use_embd_strtor   = .true.
      use_embd_tortor   = .true.

      bonded_l             = .true.
      shortnonbonded_l     = .false.
      nonbonded_l          = .true.
      PotentialAll         = .true.
      PotentialAmoeba      = .false.
      PotentialAmoeba18    = .false.
      PotentialWaterAmoeba = .false.
      PotentialWaterCharmm = .false.
      PotentialCharmm      = .false.
   end subroutine set_potential_defaults

   module subroutine parse_dynamic_cli(raw_cli, periodic, nstep, timestep_fs, &
         dump_interval_ps, ensemble, temperature, pressure, status_value)
      implicit none
      type(raw_cli_t), intent(in) :: raw_cli
      logical, intent(in) :: periodic
      integer, intent(out) :: nstep
      real(r_p), intent(out) :: timestep_fs
      real(r_p), intent(out) :: dump_interval_ps
      integer, intent(out) :: ensemble
      real(r_p), intent(out) :: temperature
      real(r_p), intent(out) :: pressure
      type(status_t), intent(out) :: status_value
      character(len=raw_cli_length) :: xyz_path
      character(len=raw_cli_length) :: values(max_raw_cli_args)
      integer :: value_count

      nstep = 0
      timestep_fs = 0.0_re_p
      dump_interval_ps = 0.0_re_p
      ensemble = 1
      temperature = 0.0_re_p
      pressure = 0.0_re_p

      call parse_common_cli(raw_cli, xyz_path, values, value_count, status_value)
      if (failed(status_value)) return

      call parse_dynamic_values(values, value_count, periodic, nstep, timestep_fs, &
         dump_interval_ps, ensemble, temperature, pressure, status_value)
   end subroutine parse_dynamic_cli

   subroutine parse_dynamic_values(values, value_count, periodic, nstep, timestep_fs, &
         dump_interval_ps, ensemble, temperature, pressure, status_value)
      implicit none
      character(*), intent(in) :: values(:)
      integer, intent(in) :: value_count
      logical, intent(in) :: periodic
      integer, intent(out) :: nstep
      real(r_p), intent(out) :: timestep_fs
      real(r_p), intent(out) :: dump_interval_ps
      integer, intent(out) :: ensemble
      real(r_p), intent(out) :: temperature
      real(r_p), intent(out) :: pressure
      type(status_t), intent(out) :: status_value
      integer :: position

      call clear_status(status_value)
      nstep = 0
      timestep_fs = 0.0_re_p
      dump_interval_ps = 0.0_re_p
      ensemble = 1
      temperature = 0.0_re_p
      pressure = 0.0_re_p
      position = 1

      call read_required_integer(values, value_count, position, nstep, &
         'You Need To Enter the Number of Steps', 'Invalid Number of Steps', status_value)
      if (failed(status_value)) return
      call read_required_real(values, value_count, position, timestep_fs, &
         'You Need To Enter the Time Step Length in Femtoseconds', &
         'Invalid Time Step Length in Femtoseconds', status_value)
      if (failed(status_value)) return
      call read_required_real(values, value_count, position, dump_interval_ps, &
         'You Need To Enter the Time Between Dumps in Picoseconds', &
         'Invalid Time Between Dumps in Picoseconds', status_value)
      if (failed(status_value)) return

      if (.not.periodic) return

      call read_required_integer(values, value_count, position, ensemble, &
         'You Need To Enter the Statistical Ensemble', 'Invalid Statistical Ensemble', status_value)
      if (failed(status_value)) return
      if (ensemble <= 0) ensemble = 1
      if (ensemble == 2 .or. ensemble == 4) then
         call read_required_real(values, value_count, position, temperature, &
            'You Need To Enter the Temperature in Kelvin', 'Invalid Temperature in Kelvin', status_value)
         if (failed(status_value)) return
         if (temperature <= 0.0_re_p) temperature = 298.0_re_p
      end if
      if (ensemble == 4) then
         call read_required_real(values, value_count, position, pressure, &
            'You Need To Enter the Pressure in Atm', 'Invalid Pressure in Atm', status_value)
         if (failed(status_value)) return
         if (pressure <= 0.0_re_p) pressure = 1.0_re_p
      end if
   end subroutine parse_dynamic_values

   subroutine read_required_integer(values, value_count, position, value, missing, invalid, status_value)
      implicit none
      character(*), intent(in) :: values(:)
      integer, intent(in) :: value_count
      integer, intent(inout) :: position
      integer, intent(out) :: value
      character(*), intent(in) :: missing, invalid
      type(status_t), intent(out) :: status_value
      integer :: io_status

      call clear_status(status_value)
      if (position > value_count) then
         call set_status(status_value, status_missing_argument, missing)
         return
      end if
      read(values(position), *, iostat=io_status) value
      if (io_status /= 0) then
         call set_status(status_value, status_invalid_argument, invalid)
         return
      end if
      position = position + 1
   end subroutine read_required_integer

   subroutine read_required_real(values, value_count, position, value, missing, invalid, status_value)
      implicit none
      character(*), intent(in) :: values(:)
      integer, intent(in) :: value_count
      integer, intent(inout) :: position
      real(r_p), intent(out) :: value
      character(*), intent(in) :: missing, invalid
      type(status_t), intent(out) :: status_value
      integer :: io_status

      call clear_status(status_value)
      if (position > value_count) then
         call set_status(status_value, status_missing_argument, missing)
         return
      end if
      read(values(position), *, iostat=io_status) value
      if (io_status /= 0) then
         call set_status(status_value, status_invalid_argument, invalid)
         return
      end if
      position = position + 1
   end subroutine read_required_real

   logical function letter_option(argument)
      implicit none
      character(*), intent(in) :: argument
      character :: letter

      letter_option = .false.
      if (len_trim(argument) < 2 .or. argument(1:1) /= '-') return
      letter = argument(2:2)
      if (letter >= 'a' .and. letter <= 'z') letter = achar(iachar(letter) - 32)
      letter_option = letter >= 'A' .and. letter <= 'Z'
   end function letter_option

   logical function key_option(argument)
      implicit none
      character(*), intent(in) :: argument

      key_option = len_trim(argument) >= 2 .and. argument(1:1) == '-'
      if (.not.key_option) return
      key_option = argument(2:2) == 'K' .or. argument(2:2) == 'k'
   end function key_option

   subroutine xyz_base(xyz_path, base_length, directory_length)
      implicit none
      character(*), intent(in) :: xyz_path
      integer, intent(out) :: base_length, directory_length
      integer :: i, path_length

      path_length = len_trim(xyz_path)
      base_length = path_length
      do i = 1, path_length
         select case (xyz_path(i:i))
         case ('/', ']', ':', '~', achar(backslash))
            base_length = path_length
         case ('.')
            base_length = i - 1
         end select
      end do
      directory_length = 0
      do i = base_length, 1, -1
         select case (xyz_path(i:i))
         case ('/', ']', ':', '~', achar(backslash))
            directory_length = i
            return
         end select
      end do
   end subroutine xyz_base

   subroutine read_key_file(configuration, status_value)
      implicit none
      type(config_t), intent(inout) :: configuration
      type(status_t), intent(out) :: status_value
      integer :: key_unit, io_status

      call clear_status(status_value)
      open(newunit=key_unit, file=configuration%key_path, status='old', &
         action='read', iostat=io_status)
      if (io_status /= 0) then
         call set_status(status_value, status_invalid_argument, 'Unable to Read Keyfile')
         return
      end if
      do
         read(key_unit, '(a)', iostat=io_status) configuration%key_line(configuration%key_count+1)
         if (io_status < 0) exit
         if (io_status /= 0) then
            close(key_unit)
            call set_status(status_value, status_invalid_argument, 'Unable to Read Keyfile')
            return
         end if
         configuration%key_count = configuration%key_count + 1
         if (configuration%key_count >= maxkey) then
            close(key_unit)
            call set_status(status_value, status_invalid_argument, 'Keyfile Too Large; Increase MAXKEY')
            return
         end if
      end do
      close(key_unit)
   end subroutine read_key_file

   subroutine echo_key_comments(configuration)
      implicit none
      type(config_t), intent(in) :: configuration
      character(len=240) :: record, keyword, comment
      integer :: i, next, length, trimtext
      logical :: header

      header = .true.
      do i = 1, configuration%key_count
         next = 1
         record = configuration%key_line(i)
         call gettext(record, keyword, next)
         call upcase(keyword)
         if (keyword(1:5) /= 'ECHO ') cycle
         comment = record(next:240)
         length = trimtext(comment)
         if (header) then
            header = .false.
            write (iout, 10)
10          format ()
         end if
         if (length == 0) then
            write (iout, 20)
20          format ()
         else
            write (iout, '(a)') comment(1:length)
         end if
      end do
   end subroutine echo_key_comments

   subroutine read_prm_file(configuration, status_value)
      implicit none
      type(config_t), intent(inout) :: configuration
      type(status_t), intent(out) :: status_value
      integer :: prm_unit, io_status

      call clear_status(status_value)
      open(newunit=prm_unit, file=configuration%prm_path, status='old', &
         action='read', iostat=io_status)
      if (io_status /= 0) then
         call set_status(status_value, status_invalid_argument, 'Unable to Read Parameter File')
         return
      end if
      do
         if (configuration%prm_count >= maxprm) then
            close(prm_unit)
            call set_status(status_value, status_invalid_argument, &
               'Parameter File Too Large; Increase MAXPRM')
            return
         end if
         read(prm_unit, '(a)', iostat=io_status) &
            configuration%prm_line(configuration%prm_count+1)
         if (io_status < 0) exit
         if (io_status /= 0) then
            close(prm_unit)
            call set_status(status_value, status_invalid_argument, 'Unable to Read Parameter File')
            return
         end if
         configuration%prm_count = configuration%prm_count + 1
      end do
      close(prm_unit)
   end subroutine read_prm_file

   subroutine validate_ndir_keywords(configuration, status_value)
      implicit none
      type(config_t), intent(in) :: configuration
      type(status_t), intent(out) :: status_value
      character(len=20) :: keyword
      character(len=raw_cli_length) :: record, value
      integer :: i, next, io_status, ndir_value

      call clear_status(status_value)
      do i = 1, configuration%key_count
         record = configuration%key_line(i)
         next = 1
         call gettext(record, keyword, next)
         call upcase(keyword)
         if (trim(keyword) /= 'NDIR-BLOCK' .and. trim(keyword) /= 'NDIR-ASYNC-BLOCK') cycle
         value = record(next:raw_cli_length)
         read(value, *, iostat=io_status) ndir_value
         if (io_status /= 0) then
            call set_status(status_value, status_invalid_argument, 'Invalid NDIR Value')
            return
         end if
      end do
   end subroutine validate_ndir_keywords

   subroutine discover_early_key_values(configuration, status_value)
      implicit none
      type(config_t), intent(inout) :: configuration
      type(status_t), intent(out) :: status_value
      character(len=20) :: keyword
      character(len=raw_cli_length) :: record, value
      integer :: i, next, io_status

      call clear_status(status_value)
      configuration%minimize_xyz_enabled = .true.
      configuration%dcd_enabled = .false.
      configuration%replicas_enabled = .false.
      configuration%pme_enabled = .false.
      configuration%pme_processes = 0
      configuration%replica_count = 0
      do i = 1, configuration%key_count
         record = configuration%key_line(i)
         next = 1
         call gettext(record, keyword, next)
         call upcase(keyword)
         value = record(next:raw_cli_length)
         select case (trim(keyword))
         case ('NO-MINIMIZE-XYZ')
            configuration%minimize_xyz_enabled = .false.
         case ('DCDIO')
            configuration%dcd_enabled = .true.
         case ('PME-PROCS')
            read(value, *, iostat=io_status) configuration%pme_processes
            if (io_status /= 0 .or. configuration%pme_processes == 0) then
               call set_status(status_value, status_invalid_argument, 'Invalid PME-PROCS Value')
               return
            end if
            configuration%pme_enabled = .true.
         case ('REPLICAS')
            read(value, *, iostat=io_status) configuration%replica_count
            if (io_status /= 0) then
               call set_status(status_value, status_invalid_argument, 'Invalid REPLICAS Value')
               return
            end if
            configuration%replicas_enabled = .true.
         end select
      end do
   end subroutine discover_early_key_values
end submodule config_cli
