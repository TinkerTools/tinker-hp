c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "echargecu" : driver for calculation of the point charge
c     energy and derivatives with respect to Cartesian coordinates on device
c
c
#ifndef TINKER_CUF
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"
      module echargecu
      use sizes   ,only: maxgrp
      use utilcu  ,only: nproc,ndir,ngrp,BLOCK_DIM,ALL_LANES,use_virial
      use utilgpu ,only: real3,real6,mdyn3_r,rpole_elt
     &            ,BLOCK_SIZE,RED_BUFF_SIZE,WARP_SIZE

      implicit none
      private
#include "atomicOp.h.f"

      public :: ecreal1_kcu,ecreal3_kcu
      contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "atomicOp.inc.f"
#include "groups.inc.f"
#include "pair_charge.f.inc"

#define __tver__ __use_grd__+__use_ene__+__use_vir__
#define __tfea__ __use_mpi__+__use_groups__+__use_lambdadyn__
#define __sufx__ 1_kcu
#include "echargecu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ __use_grd__+__use_ene__+__use_vir__
#define __tfea__ __use_mpi__+__use_shortRange__+__use_groups__+__use_lambdadyn__
#define __sufx__ 1_sht_kcu
#include "echargecu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ __use_ene__+__use_act__
#define __tfea__ __use_mpi__+__use_groups__
#define __sufx__ 3_kcu
#include "echargecu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

      end module
#endif
