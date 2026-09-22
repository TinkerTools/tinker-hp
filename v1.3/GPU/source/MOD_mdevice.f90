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
#ifdef _OPENACC
   use cudafor, only: cuda_stream_kind
#else
   integer, parameter :: cuda_stream_kind = 8
#endif
   use status, only: status_t
   implicit none
   private

   type, public :: channel_t
      character(len=12) :: name = 'unbound'
      integer :: queue = -1
      integer(cuda_stream_kind) :: stream = 0
   end type channel_t

   integer, parameter, private :: phase_new = 0
   integer, parameter, private :: phase_bound = 1
   integer, parameter, private :: phase_active = 2
   integer, parameter, private :: phase_finalized = 3
   integer, private :: phase = phase_new
   logical, private :: capabilities_selected = .false.
   type(channel_t), private :: direct_work, reciprocal_work, default_work

   public :: bind_device, initialize_device, direct_channel, reciprocal_channel, default_channel
   public :: after, synchronize, finalize_device

   interface
      module subroutine bind_device(status_value)
         type(status_t), intent(out) :: status_value
      end subroutine bind_device

      module subroutine initialize_device(host_communicator, status_value)
         integer, intent(in) :: host_communicator
         type(status_t), intent(out) :: status_value
      end subroutine initialize_device

      module function direct_channel() result(channel)
         type(channel_t) :: channel
      end function direct_channel

      module function reciprocal_channel() result(channel)
         type(channel_t) :: channel
      end function reciprocal_channel

      module function default_channel() result(channel)
         type(channel_t) :: channel
      end function default_channel

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
   end interface
end module mdevice
