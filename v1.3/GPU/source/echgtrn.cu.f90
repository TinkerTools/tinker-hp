!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     "echgtrncu" calculates the charge transfer energy and first
!     derivatives with respect to Cartesian coordinates on device
!
!
#ifndef TINKER_CUF
#define TINKER_CUF
#endif
#include "tinker_macro.h"
#include "tinker_cudart.h"
module echgtrncu
   use sizes    ,only: maxgrp
   use tinheader,only: zeror,oner
   use tinTypes ,only: real3,mdyn3_r
   use libm, only: f_abs=>TINKERHP_LIBM_ABS
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
   use libm, only: f_sqrt=>TINKERHP_LIBM_SQRT,f_exp=>TINKERHP_LIBM_EXP
#endif
   use reduction,only: RED_BUFF_SIZE
   private
   integer,parameter :: BLOCK_DIM=128,WARP_SIZE=32,WARP_SHFT=5
   public :: echgtrn1_kcu
   public :: BLOCK_DIM,WARP_SIZE,WARP_SHFT
contains
#include "image.inc.f90"
#include "midpointimage.inc.f90"
#include "groups.inc.f90"
#include "atomicOp.inc.f90"
#include "convert.inc.f90"
#include "pair_chgtrn.inc.f90"

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_lambdadyn__)
#define __sufx__ 1_kcu
#include "echgtrncu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_lambdadyn__+__use_shortRange__)
#define __sufx__ 1s_kcu
#include "echgtrncu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_lambdadyn__+__use_longRange__)
#define __sufx__ 1l_kcu
#include "echgtrncu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_groups__)
#define __sufx__ 3_kcu
#include "echgtrncu.tpl.f90"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

!#define __tver__ (__use_ene__)
!#define __tfea__ (__use_mpi__+__use_groups__)
!#define __sufx__ _kcu
!#include "echgtrncu.tpl.f90"
!#undef  __tver__
!#undef  __tfea__
!#undef  __sufx__

end module
