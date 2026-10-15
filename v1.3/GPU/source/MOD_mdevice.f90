!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module mdevice  --  GPU platform lifecycle                ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module mdevice
   use mdevice_channel, only: channel_t
   use status, only: status_t
   implicit none

   public :: channel_t

   integer, parameter :: phase_new = 0
   integer, parameter :: phase_bound = 1
   integer, parameter :: phase_active = 2
   integer, parameter :: phase_finalized = 3
   integer, protected :: phase = phase_new
   logical, protected :: capabilities_selected = .false.
   integer, protected :: gpu_sm_count = 1
   integer, protected :: cores_per_sm = 1
   type(channel_t), protected, target :: direct_work, reciprocal_work, sync_work
   type(status_t), protected :: platform_error

   public :: bind_device, initialize_device
   public :: after, synchronize, finalize_device
   public :: begin_direct_overlap, end_direct_overlap
   public :: disable_overlap
   public :: grid_dim, warning_text
   public :: publish_warning, set_gpu_config
   public :: device_number, device_name, precomp_allowed
   public :: check_launch_kernel, record_platform_error, query_platform_status

   interface
      module subroutine bind_device(status_value)
         type(status_t), intent(out) :: status_value
      end subroutine bind_device

      module subroutine initialize_device(host_communicator, status_value)
         integer, intent(in) :: host_communicator
         type(status_t), intent(out) :: status_value
      end subroutine initialize_device

      module subroutine after(consumer, producer, status_value)
         type(channel_t), intent(in) :: consumer, producer
         type(status_t), intent(out) :: status_value
      end subroutine after

      module subroutine synchronize(channel, status_value)
         type(channel_t), intent(in) :: channel
         type(status_t), intent(out) :: status_value
      end subroutine synchronize

      module subroutine finalize_device(status_value)
         type(status_t), intent(out) :: status_value
      end subroutine finalize_device

      module subroutine begin_direct_overlap(status_value)
         type(status_t), intent(out), optional :: status_value
      end subroutine begin_direct_overlap

      module subroutine end_direct_overlap(status_value)
         type(status_t), intent(out), optional :: status_value
      end subroutine end_direct_overlap

      module subroutine disable_overlap
      end subroutine disable_overlap

      module subroutine record_platform_error(platform_code, message)
         integer, intent(in) :: platform_code
         character(*), intent(in) :: message
      end subroutine record_platform_error

      module subroutine check_launch_kernel(message, wait_for_completion)
         character(*), optional, intent(in) :: message
         logical, optional, intent(in) :: wait_for_completion
      end subroutine check_launch_kernel

      module subroutine query_platform_status(status_value)
         type(status_t), intent(out) :: status_value
      end subroutine query_platform_status

      module function grid_dim(n, block_dim, opt) result(grid_size)
         integer, intent(in) :: n, block_dim
         integer, intent(in), optional :: opt
         integer :: grid_size
      end function grid_dim

      module function warning_text() result(text)
         character(len=128) :: text
      end function warning_text

      module subroutine publish_warning
      end subroutine publish_warning

      module subroutine set_gpu_config(gangs, workers, vector, modulo)
         integer, intent(in), optional :: gangs, workers, vector
         logical, intent(in), optional :: modulo
      end subroutine set_gpu_config

      module function device_number() result(value)
         integer :: value
      end function device_number

      module subroutine device_name(name)
         character(*), intent(out) :: name
      end subroutine device_name

      module function precomp_allowed(n) result(authority)
         logical :: authority
         integer, intent(in) :: n
      end function precomp_allowed
   end interface
end module mdevice
