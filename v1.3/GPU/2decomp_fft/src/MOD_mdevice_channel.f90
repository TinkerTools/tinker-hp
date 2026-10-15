!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
module mdevice_channel
#ifdef _OPENACC
   use cudafor, only: cuda_stream_kind
#else
   integer, parameter :: cuda_stream_kind = 8
#endif
   implicit none

   type, public :: channel_t
      character(len=12) :: name = 'unbound'
      integer :: queue = -1
      integer(cuda_stream_kind) :: stream = 0
   end type channel_t
end module mdevice_channel
