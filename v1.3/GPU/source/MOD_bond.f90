!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module bond  --  covalent bonds in the current structure  ##
!     ##                                                            ##
!     ################################################################
!
!
!     bk      bond stretch force constants (kcal/mole/Ang**2)
!     winbk    window object corresponding to bk
!     bl      ideal bond length values in Angstroms
!     winbl    window object corresponding to bk
!     nbond   total number of bond stretches in the system
!     nbond_pe   total number of bond stretches in the system per process element
!     nbondloc   local number of bond stretches in the system
!     ibnd    numbers of the atoms in each bond stretch
!     winibnd    window object corresponding to ibnd
!
#include "tinker_macro.h"
module bond
#ifdef USE_NVSHMEM_CUDA
   use tinTypes,only: i2dDPC=>Int2dDevPointerContainer&
      &,        rDPC=> RealDevPointerContainer
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
