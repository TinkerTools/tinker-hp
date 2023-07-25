c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module bond  --  covalent bonds in the current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     bk      bond stretch force constants (kcal/mole/Ang**2)
c     winbk    window object corresponding to bk
c     bl      ideal bond length values in Angstroms
c     winbl    window object corresponding to bk
c     nbond   total number of bond stretches in the system
c     nbond_pe   total number of bond stretches in the system per process element
c     nbondloc   local number of bond stretches in the system
c     ibnd    numbers of the atoms in each bond stretch
c     winibnd    window object corresponding to ibnd
c
#include "tinker_macro.h"
      module bond
#ifdef USE_NVSHMEM_CUDA
      use tinTypes,only: i2dDPC=>Int2dDevPointerContainer 
     &            ,        rDPC=> RealDevPointerContainer
#endif
      implicit none
      integer nbond,nbondloc
      integer  , pointer :: ibnd(:,:)
      integer :: winibnd
      real(t_p), pointer :: bk(:),bl(:),ba(:)
      integer :: winbk,winbl,winba

#ifdef USE_NVSHMEM_CUDA
      !d_*     device data type container for nvshmem feature
      !c_*     host data type container for nvshmem feature
      integer nbond_pe
      type(i2dDPC),device,pointer::d_ibnd(:)
      type(i2dDPC),   allocatable::c_ibnd(:)
      type(rDPC)  ,device,pointer::d_bk(:),d_bl(:),d_ba(:)
      type(rDPC)  ,   allocatable::c_bk(:),c_bl(:),c_ba(:)
#endif
      end
