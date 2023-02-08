c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar1  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar_cpencu" calculates the induced dipole polarization energy
c                     + Charge penetration
c     and derivatives with respect to Cartesian coordinates in a CUDA Fortran kernel
c
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"

      module epolar_cpencu
        use utilcu  ,only: nproc,ndir,BLOCK_DIM,ALL_LANES
        use utilgpu ,only: BLOCK_SIZE,RED_BUFF_SIZE
     &              ,WARP_SIZE
        use tinTypes,only: real3,real6,mdyn3_r,rpole_elt
        implicit none
        private

        public:: epreal_cpen1_kcu,epreal_cpen3_kcu
        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "atomicOp.inc.f"
#include "pair_polar_cpen.inc.f"

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_ewald__+__use_chgpen__+__use_chgflx__)
#define __sufx__ 1_kcu
#include "epolar_cpencu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_ewald__+__use_chgpen__)
#define __sufx__ 3_kcu
#include "epolar_cpencu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

        end module
#endif
