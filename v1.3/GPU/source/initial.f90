!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"

!     Initialize GPU
!
subroutine initDevice(platform_ready)
   use deriv  ,only: delambdae,delambdav,delambdaesave,&
      &delambdavsave
   use domdec ,only: nproc,rank,ranktot
   use parallel, only: get_host_comm
   use energi
   use inter  ,only: einter
   use inform
   use mamd   ,only: aMDwattype
   use mdevice,only: reciprocal_work,direct_work,init_legacy,publish_warning,&
      &set_gpu_config,device_number
   use nvshmem
   use pme, only: initialize_pme_device
   use reduction, only: create_reduction_buffers
#ifdef _OPENACC
   use interfaces,only:C_init_env
#endif
   use sizes  ,only: tinkerdebug
   use virial
   !use tinMemory
   !use utils
   implicit none
   integer,volatile:: ierr,i
   logical, intent(in), optional :: platform_ready
   logical :: platform_is_ready
   platform_is_ready = .false.
   if (present(platform_ready)) platform_is_ready = platform_ready
!
!     Inititialize gpu stuff
!
   if (.not.platform_is_ready) then
      call init_legacy
   end if
   call init_nvshmem(get_host_comm())
   !call mHalte
!
!     For Mpi debug purpose
!
   if (btest(tinkerdebug,tindGdb)) then
      if (ranktot.eq.0) then
         print*
         write(*,*) "*** Enabled Debug with [cuda-]Gdb ***"
         i = 0
         do while (i.eq.0)
            call sleep(2)
         end do
      end if
      call MPI_BARRIER(get_host_comm(),ierr)
   end if

   call publish_warning
!
!     For blocking computation of direct space
!
!$acc enter data create(aMDwattype)
!
!     create virial component on device
!
   call create_vir_device
!
!     create energy data on device
!
   call create_energi_on_device()
!$acc enter data create(einter,delambdae,delambdav,delambdaesave, &
!$acc delambdavsave)
!
!     set to zero reductions buffers for (action,energy and virial)
!
   call create_reduction_buffers(direct_work%queue)
   call initialize_pme_device(reciprocal_work%queue)
!
!     mapping parallelism for gpu execution
!
#ifdef _OPENACC
   call C_init_env(device_number(),nproc,rank,tinkerdebug)
#endif

   call set_gpu_config(gangs=80,workers=1,vector=32,modulo=.true.)

end subroutine
