c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "echgtrncu" calculates the charge transfer energy and first
c     derivatives with respect to Cartesian coordinates on device
c
c
#ifndef TINKER_CUF
#define TINKER_CUF
#endif
#include "tinker_macro.h"
#include "tinker_cudart.h"
      module echgtrncu
      use sizes    ,only: maxgrp
      use tinheader,only: zeror,oner
      use tinTypes ,only: real3,mdyn3_r
      use utilcu   ,only: BLOCK_DIM,use_virial,nproc,ngrp,f_abs
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &             ,f_sqrt,f_exp
#endif
      use utilgpu  ,only: RED_BUFF_SIZE,WARP_SIZE
      private
      public :: echgtrn1_kcu
      contains
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"
#include "convert.f.inc"
#include "pair_chgtrn.inc.f"

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__)
#define __sufx__ 1_kcu
#include "echgtrncu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_shortRange__)
#define __sufx__ 1s_kcu
#include "echgtrncu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_longRange__)
#define __sufx__ 1l_kcu
#include "echgtrncu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_groups__)
#define __sufx__ 3_kcu
#include "echgtrncu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

c#define __tver__ (__use_ene__)
c#define __tfea__ (__use_mpi__+__use_groups__)
c#define __sufx__ _kcu
c#include "echgtrncu.tpl.f"
c#undef  __tver__
c#undef  __tfea__
c#undef  __sufx__

      end module
