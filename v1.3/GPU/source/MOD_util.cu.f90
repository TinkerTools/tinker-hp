!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module utilcu  -- CUDA utility functions and params           ##
!     ##                                                                ##
!     ####################################################################
!
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
module utilcu
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
   use libm ,only: f_sign=>copysignf, f_abs=>fabsf&
      &, f_sqrt=>sqrtf, f_floor=>floorf, f_erfc=>erfcf&
      &, f_exp=>expf, f_min=>fminf, f_max=>fmaxf
   use cudadevice ,only: f_inv=>__frcp_rn
#  ifdef USE_DETERMINISTIC_REDUCTION
   use cudadevice ,only: __tp2ll_rz=>__float2ll_rz&
      &, __rp2ll_rz=>__double2ll_rz
#  endif
#else
   use libm ,only: f_sign=>copysign, f_abs=>fabs&
      &, f_floor=>floor, f_erfc=> erfc&
      &, f_min=>fmin, f_max=>fmax
   use cudadevice ,only: f_inv=>__drcp_rn
#endif
   implicit none
   integer  ,parameter :: all_lanes=Z'ffffffff'

end module
#else
 ! For portability purpose
subroutine void_utilcu
end
#endif
