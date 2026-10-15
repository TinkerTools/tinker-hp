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
   use mpi    , only: MPI_COMM_NULL, MPI_SUCCESS, MPI_Comm_rank
   use status , only: clear_status, failed, set_status, StatusDeviceError, status_invalid_argument
   use inform , only: Tinker_shellEnv, warning, set_warning
   use domdec  , only: rank, nproc, hostrank
   use lifecycle_observer, only: LifecycleDevice, record_lifecycle_event
   use tinheader, only: ti_p
   use tinMemory, only: mipk
#ifdef _OPENACC
   use openacc
   use cudafor
#endif
#ifndef _OPENACC
   integer, parameter :: cuda_stream_kind_local = 8
#endif
   implicit none

   integer ngpus
   integer :: devicenum = -1
   integer gpu_gangs, gpu_workers, gpu_vector
   integer rec_queue, dir_queue, def_queue, syn_queue
   logical mod_gangs
   integer nSMP, nSPcores, cores_SMP
#ifdef _OPENACC
   type(cudaDeviceProp) devProp
   integer(cuda_stream_kind) dir_stream, rec_stream, syn_stream, def_stream
   type(cudaEvent) dir_event, rec_event, def_event
   integer, parameter :: zero_flag = 0
#else
   integer(cuda_stream_kind_local) dir_stream, rec_stream, syn_stream, def_stream
   integer(cuda_stream_kind_local) dir_event, rec_event, def_event
#endif

!$acc declare create(gpu_gangs,gpu_workers,gpu_vector,devicenum,ngpus)

contains

#ifdef _OPENACC
   subroutine bind_gpu
      character(len=6) :: local_rank_env
      integer :: local_rank_env_status, local_rank
      integer device_type
      logical local_rank_found

      call get_environment_variable(name="PGI_ACC_DEVICE_NUM", &
         &value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         call capture_bound_device
         return
      end if

      call get_environment_variable(name="ACC_DEVICE_NUM", &
         &value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         call capture_bound_device
         return
      end if

      call get_environment_variable(name="ACC_DEVICE_HOST_RANK", &
         &value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         call capture_bound_device
         return
      end if

      call get_environment_variable(name="CUDA_VISIBLE_DEVICE", &
         &value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         read(local_rank_env, *) local_rank
         local_rank_found = .true.
      else
         local_rank_found = .false.
      end if

      if (.not.local_rank_found) then
         call get_environment_variable(name="SLURM_LOCALID", &
            &value=local_rank_env, status=local_rank_env_status)
         if (local_rank_env_status == 0) then
            read(local_rank_env, *) local_rank
            local_rank_found = .true.
         end if
      end if

      if (.not.local_rank_found) then
         call get_environment_variable(name="OMPI_COMM_WORLD_LOCAL_RANK", &
            &value=local_rank_env, status=local_rank_env_status)
         if (local_rank_env_status == 0) then
            read(local_rank_env, *) local_rank
            local_rank_found = .true.
         end if
      end if

      if (.not.local_rank_found) then
         call get_environment_variable(name="MPI_LOCALRANKID", &
            &value=local_rank_env, status=local_rank_env_status)
         if (local_rank_env_status == 0) then
            read(local_rank_env, *) local_rank
            local_rank_found = .true.
         end if
      end if

      if (.not.local_rank_found) then
         write(0,'(a,/,a)') &
            &"Warning : Could not bind GPU device before MPI_Init", &
            &"          Switched to default binding procedure"
         return
      end if

      device_type = acc_get_device_type()
      ngpus = acc_get_num_devices(device_type)
      if (local_rank.ge.ngpus) local_rank = mod(local_rank,ngpus)
      call acc_set_device_num(local_rank,device_type)
      devicenum = local_rank
   end subroutine bind_gpu

   subroutine capture_bound_device
      integer :: device_type

      device_type = acc_get_device_type()
      ngpus = acc_get_num_devices(device_type)
      devicenum = acc_get_device_num(device_type)
   end subroutine capture_bound_device

   subroutine selectDevice(status_code)
      use mpi
      implicit none
      integer, intent(out), optional :: status_code
      character*64 value,value1,value2
      character*22 name
      integer length,length1,length2,status,status1,status2
      integer device_type, device_start, host_comm, host_rank, ierr
      integer device_status, getpid, hostnm

      if (present(status_code)) status_code = 0
      host_rank = hostrank
      if (devicenum == -1) then
         device_type = acc_get_device_type()
         device_start = 0
         call get_environment_variable('ACC_DEVICE_NUM',value, &
            &length,status)
         call get_environment_variable('PGI_ACC_DEVICE_NUM',value1, &
            &length1,status1)
         call get_environment_variable('ACC_DEVICE_HOST_RANK',value2, &
            &length2,status2)
         if (status1.eq.0) then
            ngpus = 1
            read(value1,'(i)') devicenum
         else if (status.eq.0) then
            ngpus = 1
            read(value,'(i)') devicenum
         else
            if (status2.eq.0) read(value2,*) device_start
            ngpus = acc_get_num_devices(device_type)
            call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
               &0, MPI_INFO_NULL, host_comm, ierr)
            call MPI_Comm_rank(host_comm,host_rank,ierr)
            call MPI_Comm_free(host_comm,ierr)
            devicenum = mod(device_start+host_rank,ngpus)
         end if
         call acc_set_device_num(devicenum,device_type)
      end if

      device_status = 0
      call getDeviceProp(device_status)
      if (present(status_code)) status_code = device_status
      if (device_status /= 0) return
      status = hostnm(name)
      if (nproc.gt.16.or.nproc.eq.1) return
      write(0,'(a,i2,a,i2,a,i8,a,a)') "rank ", host_rank, &
         &" attach to device ", devicenum, " pid : ", getpid(), &
         &" on ", name
   end subroutine selectDevice

   subroutine initDevice_stream(status_code)
      implicit none
      integer, intent(out), optional :: status_code
      integer cuda_success, int_val, length, status
      integer high_priority, low_priority
      character*64 value

      if (present(status_code)) status_code = 0
      rec_queue = acc_async_noval
      dir_queue = acc_async_noval
      def_queue = acc_async_noval
      syn_queue = acc_async_sync
      rec_stream = acc_get_cuda_stream(acc_async_noval)
      dir_stream = rec_stream
      def_stream = rec_stream
      syn_stream = acc_get_cuda_stream(acc_async_sync)

      call Tinker_shellEnv("ASYNC_COVER",int_val,0)
      if (int_val.eq.0.or.(nproc.eq.1.and.int_val.eq.1)) return
      if (rec_stream.eq.0) then
         if (rank.eq.0) write(*,'(/,a,/,a)') &
            &'| WARNING !! cannot recover CUDA stream from OpenACC async', &
            &'             Skiping multi queues feature activation'
         return
      end if

      if (rank.eq.0) write(*,'(a,I5,a)') &
         &' ***** Asynchronous Computation Overlapping enable ', &
         &acc_async_noval, ' *****'
      cuda_success = cudaDeviceGetStreamPriorityRange(low_priority, &
         &high_priority)
      if (cuda_success.ne.0) then
         if (present(status_code)) status_code = cuda_success
         return
      end if
      dir_queue = acc_async_noval + 1
      cuda_success = cuda_success + &
         &cudastreamcreate(rec_stream,cudaStreamNonBlocking,high_priority)
      call acc_set_cuda_stream(rec_queue,rec_stream)
      if (dir_queue.ne.rec_queue) then
         cuda_success = cuda_success + &
            &cudastreamcreate(dir_stream,cudaStreamNonBlocking,low_priority)
         call acc_set_cuda_stream(dir_queue,dir_stream)
      else
         dir_stream = rec_stream
         def_stream = rec_stream
      end if
      if (cuda_success.ne.0) then
         if (present(status_code)) status_code = cuda_success
         return
      end if
      cuda_success = cuda_success + cudaEventCreateWithFlags( &
         &dir_event,cudaEventDisableTiming)
      cuda_success = cuda_success + cudaEventCreateWithFlags( &
         &rec_event,cudaEventDisableTiming)
      cuda_success = cuda_success + cudaEventCreateWithFlags( &
         &def_event,cudaEventDisableTiming)
      if (present(status_code)) status_code = cuda_success
   end subroutine initDevice_stream

   subroutine syncWithStream(stream,status_code)
      integer(cuda_stream_kind), optional :: stream
      integer, intent(out), optional :: status_code
      integer(cuda_stream_kind) stream_

      if (present(status_code)) status_code = 0
      if (present(stream)) then
         stream_ = stream
      else
         stream_ = rec_stream
      end if
      if (present(status_code)) status_code = cudaStreamSynchronize(stream_)
   end subroutine syncWithStream

   subroutine Finalize_async_recover
      implicit none

      if (dir_queue.ne.rec_queue) then
         if (rank.eq.0) write(*,'(/,1x,a,/)') &
            &'***** Asynchronous Computation Overlapping disabled *****'
         dir_stream = rec_stream
         def_stream = rec_stream
         dir_queue = rec_queue
         def_queue = rec_queue
      end if
   end subroutine Finalize_async_recover

   subroutine stream_wait_async(wait_s,async_s,event_,status_code)
      implicit none
      integer(cuda_stream_kind), intent(in) :: wait_s, async_s
      type(cudaEvent), optional :: event_
      integer, intent(out), optional :: status_code
      type(cudaEvent) :: event
      integer cuda_success

      if (present(status_code)) status_code = 0
      event = def_event
      if (present(event_)) event = event_
      cuda_success = cudaEventRecord(event,wait_s)
      cuda_success = cuda_success + cudaStreamWaitEvent(async_s,event,zero_flag)
      if (present(status_code)) status_code = cuda_success
   end subroutine stream_wait_async

   subroutine start_dir_stream_cover(status_code)
      implicit none
      integer, intent(out), optional :: status_code
      integer cuda_success

      cuda_success = cudaEventRecord(dir_event,rec_stream)
      cuda_success = cuda_success + cudaStreamWaitEvent( &
         &dir_stream,dir_event,zero_flag)
      if (present(status_code)) status_code = cuda_success
   end subroutine start_dir_stream_cover

   subroutine end_dir_stream_cover(status_code)
      implicit none
      integer, intent(out), optional :: status_code
      integer cuda_success

      cuda_success = cudaEventRecord(rec_event,dir_stream)
      cuda_success = cuda_success + cudaStreamWaitEvent( &
         &rec_stream,rec_event,zero_flag)
      if (present(status_code)) status_code = cuda_success
   end subroutine end_dir_stream_cover

   subroutine getDeviceProp(status_code)
      implicit none
      integer, intent(out), optional :: status_code
      integer ierr

      ierr = CUDAGETDEVICEPROPERTIES(devProp,devicenum)
      if (ierr.ne.0) then
         if (present(status_code)) status_code = ierr
         return
      end if
      nSMP = devProp%multiProcessorCount
      select case (devProp%major)
      case (2)
         if (devProp%minor.eq.1) then
            cores_SMP = 48
         else
            cores_SMP = 32
         end if
      case (3)
         cores_SMP = 192
      case (5)
         cores_SMP = 128
      case (6)
         if (devProp%minor.eq.1) then
            cores_SMP = 128
         else if (devProp%minor.eq.0) then
            cores_SMP = 64
         else
            cores_SMP = 64
            ierr = 1
         end if
      case (7,8)
         cores_SMP = 64
      case (9)
         cores_SMP = 128
      case default
         cores_SMP = 128
         ierr = 1
      end select
      nSPcores = cores_SMP*nSMP
      if (present(status_code)) status_code = ierr
   end subroutine getDeviceProp
#endif

   logical function allow_precomp(n) result(auth)
      implicit none
      integer, intent(in) :: n
#ifdef _OPENACC
      integer(mipk) free,total,esti
      integer ierr
      esti = 200*n*7*ti_p
      ierr = cudaMemGetInfo(free,total)
      auth = (esti.lt.(free*9)/10)
#else
      auth = .true.
#endif
   end function allow_precomp

   integer function get_GridDim(n,bDim,opt) result(GDim)
      implicit none
      integer, intent(in) :: bDim,n
      integer, optional :: opt
      integer, parameter :: maxB=2**14, warp_width=32

      if (present(opt)) then
         GDim = min(((n-1)/ishft(bDim,1))+1,maxB)
      else
         GDim = min(((n-1)*warp_width/bDim)+1,maxB)
      end if
   end function get_GridDim

   module subroutine bind_device(status_value)
      implicit none
      type(status_t), intent(out) :: status_value
      integer :: device_status

      if (phase == phase_finalized) then
         call set_status(status_value, status_invalid_argument, 'A finalized device cannot be rebound')
         return
      end if
      if (phase == phase_new) then
#ifdef _OPENACC
         call bind_gpu
         if (devicenum >= 0) then
            device_status = 0
            call getDeviceProp(device_status)
            if (device_status /= 0) then
               call set_status(status_value, StatusDeviceError, 'Unable to query GPU device properties')
               return
            end if
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
         device_status = 0
         call selectDevice(device_status)
         if (device_status /= 0) then
            call set_status(status_value, StatusDeviceError, 'Unable to query GPU device properties')
            return
         end if
      end if
#else
      ngpus = 0
      nSMP = 1
      nSPcores = 1
      cores_SMP = 1
#endif
      gpu_sm_count = nSMP
      cores_per_sm = cores_SMP
#ifdef _OPENACC
      device_status = 0
      call initDevice_stream(device_status)
      if (device_status /= 0) then
         call release_platform_resources
         call set_status(status_value, StatusDeviceError, 'Unable to initialize GPU execution channels')
         return
      end if
#endif
      reciprocal_work = channel_t('reciprocal', rec_queue, rec_stream)
      direct_work     = channel_t('direct',     dir_queue, dir_stream)
      default_work    = channel_t('default',    rec_queue, rec_stream)
      sync_work       = channel_t('synchronous',syn_queue, syn_stream)
      phase = phase_active
      call clear_status(status_value)
   end subroutine initialize_device

   module subroutine after(consumer, producer, status_value)
      implicit none
      type(channel_t), intent(in) :: consumer, producer
      type(status_t), intent(out) :: status_value
      integer :: dep_status

      if (phase /= phase_active) then
         call set_status(status_value, status_invalid_argument, 'A device dependency requires an active device')
         return
      end if
      if (consumer%stream /= producer%stream) then
#ifdef _OPENACC
         dep_status = 0
         select case (consumer%name)
         case ('direct')
            call stream_wait_async(producer%stream, consumer%stream, rec_event, dep_status)
         case ('reciprocal')
            call stream_wait_async(producer%stream, consumer%stream, dir_event, dep_status)
         case default
            call stream_wait_async(producer%stream, consumer%stream, def_event, dep_status)
         end select
         if (dep_status /= 0) then
            call set_status(status_value, StatusDeviceError, 'Unable to establish GPU stream dependency')
            return
         end if
#endif
      end if
      call clear_status(status_value)
   end subroutine after

   module subroutine synchronize(channel, status_value)
      implicit none
      type(channel_t), intent(in) :: channel
      type(status_t), intent(out) :: status_value
      integer :: sync_status

      if (phase /= phase_active) then
         call set_status(status_value, status_invalid_argument, 'Device synchronization requires an active device')
         return
      end if
#ifdef _OPENACC
      sync_status = 0
      call syncWithStream(channel%stream, sync_status)
      if (sync_status /= 0) then
         call set_status(status_value, StatusDeviceError, 'Unable to synchronize GPU execution channel')
         return
      end if
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
      if (failed(platform_error)) then
         status_value = platform_error
      else if (cleanup_status /= 0) then
         call set_status(status_value, status_invalid_argument, 'Unable to finalize GPU execution channels')
      else
         call clear_status(status_value)
      end if
      call clear_status(platform_error)
   end subroutine finalize_device

   module subroutine record_platform_error(platform_code, message)
      implicit none
      integer, intent(in) :: platform_code
      character(*), intent(in) :: message
      character(len=240) :: detail

      if (failed(platform_error)) return
      write(detail,'(A," (platform code ",I0,")")') trim(message), platform_code
      call set_status(platform_error, StatusDeviceError, trim(detail))
   end subroutine record_platform_error

   module subroutine query_platform_status(status_value)
      implicit none
      type(status_t), intent(out) :: status_value

      status_value = platform_error
   end subroutine query_platform_status

   module subroutine begin_direct_overlap(status_value)
      implicit none
      type(status_t), intent(out), optional :: status_value
      integer :: overlap_status

      overlap_status = 0
      call start_dir_stream_cover(overlap_status)
      if (overlap_status /= 0) call record_platform_error(overlap_status, &
         &'Unable to begin direct GPU stream overlap')
      if (present(status_value)) then
         if (overlap_status /= 0) then
            call set_status(status_value, StatusDeviceError, 'Unable to begin direct GPU stream overlap')
         else
            call clear_status(status_value)
         end if
      end if
   end subroutine begin_direct_overlap

   module subroutine end_direct_overlap(status_value)
      implicit none
      type(status_t), intent(out), optional :: status_value
      integer :: overlap_status

      overlap_status = 0
      call end_dir_stream_cover(overlap_status)
      if (overlap_status /= 0) call record_platform_error(overlap_status, &
         &'Unable to end direct GPU stream overlap')
      if (present(status_value)) then
         if (overlap_status /= 0) then
            call set_status(status_value, StatusDeviceError, 'Unable to end direct GPU stream overlap')
         else
            call clear_status(status_value)
         end if
      end if
   end subroutine end_direct_overlap

   module subroutine disable_overlap
      implicit none
#ifdef _OPENACC
      call Finalize_async_recover
      reciprocal_work = channel_t('reciprocal', rec_queue, rec_stream)
      direct_work     = channel_t('direct',     dir_queue, dir_stream)
      default_work    = channel_t('default',    def_queue, def_stream)
      sync_work       = channel_t('synchronous',syn_queue, syn_stream)
#endif
   end subroutine disable_overlap

   module function grid_dim(n, block_dim, opt) result(grid_size)
      implicit none
      integer, intent(in) :: n, block_dim
      integer, intent(in), optional :: opt
      integer :: grid_size

      if (present(opt)) then
         grid_size = get_GridDim(n, block_dim, opt)
      else
         grid_size = get_GridDim(n, block_dim)
      end if
   end function grid_dim

   module function warning_text() result(text)
      implicit none
      character(len=128) :: text
      text = warning
   end function warning_text

   module subroutine init_legacy(platform_ready)
      implicit none
      logical, intent(in), optional :: platform_ready
      logical :: platform_rdy
      integer :: device_status

      platform_rdy = .false.
      if (present(platform_ready)) platform_rdy = platform_ready
      if (.not.platform_rdy) then
         ngpus = 0
         nSMP = 1
         nSPcores = 1
         cores_SMP = 1
#ifdef _OPENACC
         device_status = 0
         call selectDevice(device_status)
         if (device_status /= 0) then
            call record_platform_error(device_status, 'Unable to query GPU device properties')
            return
         end if
         device_status = 0
         call initDevice_stream(device_status)
         if (device_status /= 0) then
            call record_platform_error(device_status, 'Unable to initialize GPU execution channels')
            return
         end if
#endif
      end if
   end subroutine init_legacy

   module subroutine publish_warning
      implicit none
      call set_warning
   end subroutine publish_warning

   module subroutine set_gpu_config(gangs, workers, vector, modulo)
      implicit none
      integer, intent(in), optional :: gangs, workers, vector
      logical, intent(in), optional :: modulo

      if (present(gangs)) gpu_gangs = gangs
      if (present(workers)) gpu_workers = workers
      if (present(vector)) gpu_vector = vector
      if (present(modulo)) mod_gangs = modulo
!$acc update device(gpu_gangs,gpu_workers,gpu_vector)
   end subroutine set_gpu_config

   module function device_number() result(value)
      implicit none
      integer :: value
      value = devicenum
   end function device_number

   module subroutine device_name(name)
      implicit none
      character(*), intent(out) :: name
      name = ''
#ifdef _OPENACC
      name = devProp%name
#endif
   end subroutine device_name

   module function precomp_allowed(n) result(authority)
      implicit none
      logical :: authority
      integer, intent(in) :: n
      authority = allow_precomp(n)
   end function precomp_allowed

   subroutine release_platform_resources(status_code)
      implicit none
      integer, intent(out), optional :: status_code
      integer :: cuda_status, cleanup_status

      call record_lifecycle_event(LifecycleDevice)
      cleanup_status = 0

#ifdef _OPENACC
      if (dir_queue /= rec_queue .and. rec_stream /= 0) then
         call syncWithStream(rec_stream, cuda_status)
         if (cuda_status /= 0) cleanup_status = cuda_status
         call syncWithStream(dir_stream, cuda_status)
         if (cuda_status /= 0 .and. cleanup_status == 0) cleanup_status = cuda_status
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
      gpu_sm_count = 1
      cores_per_sm = 1
      direct_work = channel_t()
      reciprocal_work = channel_t()
      default_work = channel_t()
      sync_work    = channel_t()
   end subroutine release_platform_resources

end submodule mdevice_lifecycle
