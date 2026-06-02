!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine epolar1  --  polarization energy & derivs  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "epolar_cpencu" calculates the induced dipole polarization energy
!                     + Charge penetration
!     and derivatives with respect to Cartesian coordinates in a CUDA Fortran kernel
!
!
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"

module epolar_cpencu
   use utilcu  ,only: nproc,ndir,BLOCK_DIM,ALL_LANES
   use utilgpu ,only: BLOCK_SIZE,RED_BUFF_SIZE&
      &,WARP_SIZE
   use tinTypes,only: real3,real6,mdyn3_r,rpole_elt
   implicit none
   private

   public:: epreal_cpen1_kcu,epreal_cpen3_kcu
contains

#include "convert.inc.f90"
#include "image.inc.f90"
#include "midpointimage.inc.f90"
#include "atomicOp.inc.f90"
#include "pair_polar_cpen.inc.f90"

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_ewald__+__use_chgpen__+__use_chgflx__)
#define __sufx__ 1_kcu
#include "epolar_cpencu.tpl.f90"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_ewald__+__use_chgpen__)
#define __sufx__ 3_kcu
#include "epolar_cpencu.tpl.f90"
#undef __tver__
#undef __tfea__
#undef __sufx__

end module
#endif
