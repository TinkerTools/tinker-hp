c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module angle  --  bond angles within the current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     ak       harmonic angle force constant (kcal/mole/rad**2)
c     winak    window object corresponding to ak
c     anat     ideal bond angle or phase shift angle (degrees)
c     winanat    window object corresponding to anat
c     afld     periodicity for Fourier bond angle term
c     winafld    window object corresponding to afld
c     nangle   total number of bond angles in the system
c     iang     numbers of the atoms in each bond angle
c     winiang    window object corresponding to iang
c     nangle_pe   total number of bond angles per process element in the system
c     nangleloc numbers of bond angles in the local domain
c     angleloc correspondance between global and local bond angles
c
c
#include "tinker_precision.h"
      module angle
#ifdef USE_NVSHMEM_CUDA
      use tinTypes,only: i2dDPC=>Int2dDevPointerContainer
     &            ,        rDPC=> RealDevPointerContainer
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
