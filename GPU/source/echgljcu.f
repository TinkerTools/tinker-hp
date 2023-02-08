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
      module eChgLjcu
        use cudafor
        use tinheader,only: ti_p,oner
        use tintypes ,only: real3
        use sizes    ,only: maxclass,maxvalue
        use utilcu   ,only: ndir,VDW_BLOCK_DIM,ALL_LANES,use_virial
     &               ,skipvdw12
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &               ,f_sqrt
#endif
        use utilgpu  ,only: inf,BLOCK_SIZE,RED_BUFF_SIZE
        use vdw      ,only: vdw_lcut2,vdweAbsurd
        use vdwpot   ,only: ARITHMETIC_RL,GEOMETRIC_RL

        implicit none
        private
#include "atomicOp.h.f"

        !procedure(echg_lj1_kcu_v0),pointer,public :: echg_lj1_kcu 
        public:: echg_lj1_kcu_v0,echg_lj1_kcu_v1

        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "atomicOp.inc.f"
#include "pair_lj.inc.f"
#include "pair_charge.f.inc"

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

        ! Fill with a templates file for optimisations purposes

#define __tver__ __use_grd__
#define __tfea__ __radepsOpt__
#define __sufx__ _v0
#include "echgljcu.tpl.f"

#undef __tver__ 
#undef __tfea__
#undef __sufx__
#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__)
#define __sufx__ _v1
#include "echgljcu.tpl.f"

      end module
