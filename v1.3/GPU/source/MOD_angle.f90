!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module angle  --  bond angles within the current structure  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     ak       harmonic angle force constant (kcal/mole/rad**2)
!     winak    window object corresponding to ak
!     anat     ideal bond angle or phase shift angle (degrees)
!     winanat    window object corresponding to anat
!     afld     periodicity for Fourier bond angle term
!     winafld    window object corresponding to afld
!     nangle   total number of bond angles in the system
!     iang     numbers of the atoms in each bond angle
!     winiang    window object corresponding to iang
!     nangle_pe   total number of bond angles per process element in the system
!     nangleloc numbers of bond angles in the local domain
!     angleloc correspondance between global and local bond angles
!
!
#include "tinker_macro.h"
module angle
#ifdef USE_NVSHMEM_CUDA
   use tinTypes,only: i2dDPC=>Int2dDevPointerContainer&
      &,        rDPC=> RealDevPointerContainer
#endif
   implicit none
   integer nangle,nangle_pe,nangleloc
   !DIR$ ATTRIBUTES ALIGN:64 :: angleloc
   integer, allocatable :: angleloc(:)
   integer, pointer :: iang(:,:)
   real(t_p), pointer ::  ak(:),anat(:),afld(:)
   integer :: winiang,winak,winanat,winafld

#ifdef USE_NVSHMEM_CUDA
   type(rDPC)  ,device,pointer::d_ak(:),d_anat(:),d_afld(:)
   type(rDPC)  ,   allocatable::c_ak(:),c_anat(:),c_afld(:)
   type(i2dDPC),device,pointer::d_iang(:)
   type(i2dDPC),   allocatable::c_iang(:)
#endif
end
