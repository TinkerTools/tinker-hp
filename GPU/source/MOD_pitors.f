c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  module pitors  --  pi-orbital torsions in the current structure  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     npitors   total number of pi-orbital torsional interactions
c     npitorsloc   local number of pi-orbital torsional interactions
c     nbpitors  number of pi-orbital torsional interactions before each atom
c     winnbpitors    window object corresponding to nbpitors
c     kpit      2-fold pi-orbital torsional force constants
c     winkpit    window object corresponding to kpit
c     ipit      numbers of the atoms in each pi-orbital torsion
c     winipit    window object corresponding to ipit
c
c
#include "tinker_precision.h"
      module pitors
#ifdef USE_NVSHMEM_CUDA
      use tinTypes,only: i2dDPC=>Int2dDevPointerContainer 
     &            ,        iDPC=>  IntDevPointerContainer 
     &            ,        rDPC=> RealDevPointerContainer
#endif
      implicit none
      integer npitors,npitorsloc
      integer, pointer :: ipit(:,:), nbpitors(:)
      real(t_p), pointer :: kpit(:)
      integer :: winipit,winnbpitors,winkpit

#ifdef USE_NVSHMEM_CUDA
      ! nvshmem container of pitors data
      integer npitors_pe
      type(i2dDPC),device,pointer::d_ipit(:)
      type(i2dDPC),   allocatable::c_ipit(:)
      type(iDPC)  ,device,pointer::d_nbpitors(:)
      type(iDPC)  ,   allocatable::c_nbpitors(:)
      type(rDPC)  ,device,pointer::d_kpit(:)
      type(rDPC)  ,   allocatable::c_kpit(:)
#endif
      end
