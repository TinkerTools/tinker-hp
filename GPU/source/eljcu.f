c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj1cu  --  Lennard-Jones energy & derivatives ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the Lennard-Jones van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c
#define TINKER_CUF
#include  "tinker_macro.h"
      module eljcu
        use cudafor
        use tinheader,only: ti_p
        use tintypes ,only: real3
        use sizes    ,only: maxclass,maxvalue
        use utilcu   ,only: ndir,ngrp,VDW_BLOCK_DIM,ALL_LANES,use_virial
     &               ,skipvdw12,vcouple
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &               ,f_sqrt
#endif
        use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE
        use sizes    ,only: maxgrp
        use vdw      ,only: vdw_lcut2,vdweAbsurd
        use vdwpot   ,only: ARITHMETIC_RL,GEOMETRIC_RL

        implicit none
        private
        integer(1) one1,two1
#include "atomicOp.h.f"

        parameter( one1=1, two1=2 )

        public :: lj1_kcu,lj3_kcu 

        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "atomicOp.inc.f"
#include "groups.inc.f"
#include "pair_lj.inc.f"

        M_subroutine 
     &        get_rad(radi,radk,rv,rule)
        implicit none
        integer  ,intent(in)::rule
        real(t_p),intent(in)::radi,radk
        real(t_p) rv
        rv = merge(radi+radk,2.0_ti_p*f_sqrt(radi*radk)
     &            ,rule.eq.ARITHMETIC_RL)
        end subroutine
        M_subroutine 
     &        get_eps(epsi,epsk,ep,rule)
        implicit none
        integer  ,intent(in)::rule
        real(t_p),intent(in)::epsi,epsk
        real(t_p) ep
        ep = merge(f_sqrt(epsi*epsk),0.5_ti_p*(epsi+epsk)
     &            ,rule.eq.GEOMETRIC_RL)
        end subroutine

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__)
#define __sufx__ 1_kcu
#include "eljcu.tpl.f"
#undef __sufx__
#undef __tfea__
#undef __tver__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__+__use_shortRange__)
#define __sufx__ 1_shr_kcu
#include "eljcu.tpl.f"
#undef __sufx__
#undef __tfea__
#undef __tver__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__+__use_longRange__)
#define __sufx__ 1_lgr_kcu
#include "eljcu.tpl.f"
#undef __sufx__
#undef __tfea__
#undef __tver__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__)
#define __sufx__ 3_kcu
#include "eljcu.tpl.f"
#undef __sufx__
#undef __tfea__
#undef __tver__

      end module
