!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     "empole_cpencu" : driver for calculation of the multipole and dipole polarization
!     energy and derivatives with respect to Cartesian coordinates on device
!
!
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"
module empole_cpencu
   use utilcu   ,only: nproc,ngrp,BLOCK_DIM
   use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE,WARP_SIZE
   use tinTypes ,only: real3,real6,mdyn3_r,rpole_elt
   use tinheader,only: zeror,oner
   use sizes    ,only: maxgrp
   implicit none
   private
#include "atomicOp.h.f90"

   public:: emreal_cpen1_kcu, emreal_cpen1s_kcu
contains

#include "convert.inc.f90"
#include "image.inc.f90"
#include "midpointimage.inc.f90"
#include "groups.inc.f90"
#include "atomicOp.inc.f90"
#include "pair_mpole1.inc.f90"

   M_subroutine mdyn3r_zr(r3)
   type(mdyn3_r),intent(out)::r3
   r3%x=0;r3%y=0;r3%z=0;
end subroutine
M_subroutine real3_zr(r3)
type(real3),intent(out)::r3
r3%x=0;r3%y=0;r3%z=0;
end subroutine

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_chgpen__+__use_chgflx__)
#define __sufx__ 1_kcu
#include "empole_cpencu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_chgpen__+__use_shortRange__+__use_chgflx__)
#define __sufx__ 1s_kcu
#include "empole_cpencu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_chgpen__+__use_longRange__+__use_chgflx__)
#define __sufx__ 1l_kcu
#include "empole_cpencu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_chgpen__)
#define __sufx__ 3_kcu
#include "empole_cpencu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

end module
#endif
