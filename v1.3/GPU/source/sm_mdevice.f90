!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  submodule mdevice_lifecycle  --  GPU platform lifecycle    ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
submodule (mdevice) mdevice_lifecycle
   use mpi, only: MPI_COMM_NULL, MPI_SUCCESS, MPI_Comm_rank
   use status, only: clear_status, set_status, status_invalid_argument
   use utilgpu, only: bind_gpu, selectDevice, initDevice_stream, &
      &rec_queue, dir_queue, def_queue, rec_stream, dir_stream, def_stream, &
      &stream_wait_async, syncWithStream, ngpus, nSMP, nSPcores, cores_SMP, &
      &dir_event, rec_event, def_event, devicenum, getDeviceProp
#ifdef _OPENACC
   use cudafor, only: cudaEventDestroy, cudaStreamDestroy
#endif
   implicit none

contains

   module subroutine bind_device(status_value)
      implicit none
      type(status_t), intent(out) :: status_value

      if (phase == phase_finalized) then
         call set_status(status_value, status_invalid_argument, 'A finalized device cannot be rebound')
         return
      end if
      if (phase == phase_new) then
#ifdef _OPENACC
         call bind_gpu
         if (devicenum >= 0) then
            call getDeviceProp
            capabilities_selected = .true.
         end if
#endif
         phase = phase_bound
      end if
      call clear_status(status_value)
   end subroutine bind_device

   module subroutine initialize_device(host_communicator, status_value)
      implicit none
      integer, intent(in) :: host_communicator
      type(status_t), intent(out) :: status_value
      integer :: ierr, host_rank, device_status

      if (phase /= phase_bound) then
         call set_status(status_value, status_invalid_argument, 'Device binding must precede post-MPI initialization')
         return
      end if
      if (host_communicator == MPI_COMM_NULL) then
         call set_status(status_value, status_invalid_argument, 'Device initialization requires a host communicator')
         return
      end if
      call MPI_Comm_rank(host_communicator, host_rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         call set_status(status_value, status_invalid_argument, 'Unable to inspect the host communicator')
         return
      end if

#ifdef _OPENACC
      if (.not.capabilities_selected) then
         ngpus = 0
         nSMP = 1
         nSPcores = 1
         cores_SMP = 1
         call selectDevice
      end if
#else
      ngpus = 0
      nSMP = 1
      nSPcores = 1
      cores_SMP = 1
#endif
#ifdef _OPENACC
      device_status = 0
      call initDevice_stream(device_status)
      if (device_status /= 0) then
         call release_platform_resources
         call set_status(status_value, status_invalid_argument, 'Unable to initialize GPU execution channels')
         return
      end if
#endif
      reciprocal_work = channel_t('reciprocal', rec_queue, rec_stream)
      direct_work = channel_t('direct', dir_queue, dir_stream)
      default_work = channel_t('default', rec_queue, rec_stream)
      call publish_legacy_platform_context
      phase = phase_active
      call clear_status(status_value)
   end subroutine initialize_device

   module function direct_channel() result(channel)
      implicit none
      type(channel_t) :: channel

      if (phase < phase_active) then
         channel = channel_t('direct',dir_queue,dir_stream)
      else
         channel = direct_work
      end if
   end function direct_channel

   module function reciprocal_channel() result(channel)
      implicit none
      type(channel_t) :: channel

      if (phase < phase_active) then
         channel = channel_t('reciprocal',rec_queue,rec_stream)
      else
         channel = reciprocal_work
      end if
   end function reciprocal_channel

   module function default_channel() result(channel)
      implicit none
      type(channel_t) :: channel

      if (phase < phase_active) then
         channel = channel_t('default',def_queue,def_stream)
      else
         channel = default_work
      end if
   end function default_channel

   module subroutine after(consumer, producer, status_value)
      implicit none
      type(channel_t), intent(in) :: consumer, producer
      type(status_t), intent(out) :: status_value

      if (phase /= phase_active) then
         call set_status(status_value, status_invalid_argument, 'A device dependency requires an active device')
         return
      end if
      if (consumer%stream /= producer%stream) then
#ifdef _OPENACC
         call stream_wait_async(producer%stream, consumer%stream)
#endif
      end if
      call clear_status(status_value)
   end subroutine after

   module subroutine synchronize(channel, status_value)
      implicit none
      type(channel_t), intent(in) :: channel
      type(status_t), intent(out) :: status_value

      if (phase /= phase_active) then
         call set_status(status_value, status_invalid_argument, 'Device synchronization requires an active device')
         return
      end if
#ifdef _OPENACC
      call syncWithStream(channel%stream)
#endif
      call clear_status(status_value)
   end subroutine synchronize

   module subroutine finalize_device(status_value)
      implicit none
      type(status_t), intent(out) :: status_value
      integer :: cleanup_status

      if (phase == phase_finalized) then
         call clear_status(status_value)
         return
      end if
      if (phase /= phase_active) then
         call set_status(status_value, status_invalid_argument, 'Only an active device can be finalized')
         return
      end if
      call release_platform_resources(cleanup_status)
      phase = phase_finalized
      if (cleanup_status /= 0) then
         call set_status(status_value, status_invalid_argument, 'Unable to finalize GPU execution channels')
      else
         call clear_status(status_value)
      end if
   end subroutine finalize_device

   subroutine release_platform_resources(status_code)
      implicit none
      integer, intent(out), optional :: status_code
      integer :: cuda_status, cleanup_status

      cleanup_status = 0

#ifdef _OPENACC
      if (dir_queue /= rec_queue .and. rec_stream /= 0) then
         call syncWithStream(rec_stream)
         call syncWithStream(dir_stream)
         cuda_status = cudaEventDestroy(dir_event)
         if (cuda_status /= 0 .and. cleanup_status == 0) cleanup_status = cuda_status
         cuda_status = cudaEventDestroy(rec_event)
         if (cuda_status /= 0 .and. cleanup_status == 0) cleanup_status = cuda_status
         cuda_status = cudaEventDestroy(def_event)
         if (cuda_status /= 0 .and. cleanup_status == 0) cleanup_status = cuda_status
         cuda_status = cudaStreamDestroy(dir_stream)
         if (cuda_status /= 0 .and. cleanup_status == 0) cleanup_status = cuda_status
         cuda_status = cudaStreamDestroy(rec_stream)
         if (cuda_status /= 0 .and. cleanup_status == 0) cleanup_status = cuda_status
      end if
#endif
      if (present(status_code)) status_code = cleanup_status
      rec_queue = -1
      dir_queue = -1
      def_queue = -1
      rec_stream = 0
      dir_stream = 0
      def_stream = 0
      direct_work = channel_t()
      reciprocal_work = channel_t()
      default_work = channel_t()
   end subroutine release_platform_resources

   subroutine publish_legacy_platform_context
      implicit none

      ! The compatibility adapter republishes only mdevice's queue and stream values.
      rec_queue = reciprocal_work%queue
      dir_queue = direct_work%queue
      def_queue = default_work%queue
      def_stream = rec_stream
   end subroutine publish_legacy_platform_context
end submodule mdevice_lifecycle
