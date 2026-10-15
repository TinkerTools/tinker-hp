!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  CUDA submodule mdevice_diagnostics -- launch diagnostics  ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
submodule (mdevice) mdevice_diagnostics
   use cudafor, only: cudaDeviceSynchronize, cudaGetLastError, cudaSuccess
   use sizes, only: tinkerdebug
   implicit none

contains

   module subroutine check_launch_kernel(message, wait_for_completion)
      implicit none
      character(*), optional, intent(in) :: message
      logical, optional, intent(in) :: wait_for_completion
      integer :: error_code
      logical :: synchronize

      synchronize = tinkerdebug > 0
      if (present(wait_for_completion)) synchronize = synchronize .or. wait_for_completion
      if (synchronize) then
         error_code = cudaDeviceSynchronize()
      else
         error_code = cudaGetLastError()
      end if
      if (error_code /= cudaSuccess) then
         if (present(message)) then
            call record_platform_error(error_code, 'CUDA launch: '//trim(message))
         else
            call record_platform_error(error_code, 'CUDA launch')
         end if
      end if
   end subroutine check_launch_kernel

end submodule mdevice_diagnostics
