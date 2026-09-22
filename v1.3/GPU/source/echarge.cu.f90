!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     "echargecu" : driver for calculation of the point charge
!     energy and derivatives with respect to Cartesian coordinates on device
!
!
#ifndef TINKER_CUF
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"
module echargecu
   use sizes   ,only: maxgrp
   use utilcu  ,only: ALL_LANES
   use tinTypes,only: real3,real6,mdyn3_r,rpole_elt
   use reduction,only: RED_BUFF_SIZE

   implicit none
   private
   integer,parameter :: BLOCK_DIM=128,BLOCK_SIZE=32,WARP_SIZE=32,WARP_SHFT=5
#include "atomicOp.h.f90"

   public :: ecreal1_kcu,ecreal3_kcu
   public :: BLOCK_DIM,BLOCK_SIZE,WARP_SIZE,WARP_SHFT
contains

#include "convert.inc.f90"
#include "image.inc.f90"
#include "midpointimage.inc.f90"
#include "atomicOp.inc.f90"
#include "groups.inc.f90"
#include "pair_charge.inc.f90"

#define __tver__ __use_grd__+__use_ene__+__use_vir__
#define __tfea__ __use_mpi__+__use_groups__+__use_lambdadyn__
#define __sufx__ 1_kcu
#include "echargecu.tpl.f90"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ __use_grd__+__use_ene__+__use_vir__
#define __tfea__ __use_mpi__+__use_shortRange__+__use_groups__+__use_lambdadyn__
#define __sufx__ 1_sht_kcu
#include "echargecu.tpl.f90"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ __use_ene__+__use_act__
#define __tfea__ __use_mpi__+__use_groups__
#define __sufx__ 3_kcu
#include "echargecu.tpl.f90"
#undef __tver__
#undef __tfea__
#undef __sufx__

end module
#endif
