!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module strbnd  --  stretch-bends in the current structure  ##
!     ##                                                             ##
!     #################################################################
!
!
!     sbk       force constants for stretch-bend terms
!     winsbk    window object corresponding to sbk
!     nstrbnd   total number of stretch-bend interactions
!     nstrbndloc   local number of stretch-bend interactions
!     nbstrbnd   number of stretch-bend interactions before each atom
!     winnbstrbnd    window object corresponding to nbstrbnd
!     isb       angle and bond numbers used in stretch-bend
!     winisb    window object corresponding to isb
!
!
#include "tinker_macro.h"
module strbnd
#ifdef USE_NVSHMEM_CUDA
   use tinTypes,only: i2dDPC=> Int2dDevPointerContainer&
      &,        iDPC=>   IntDevPointerContainer&
      &,      r2dDPC=>Real2dDevPointerContainer
#endif
   implicit none
   integer nstrbnd,nstrbndloc
   integer, pointer :: isb(:,:),nbstrbnd(:)
   real(t_p), pointer ::  sbk(:,:)
   integer :: winisb,winnbstrbnd,winsbk

#ifdef USE_NVSHMEM_CUDA
   type(i2dDPC),device,pointer::d_isb(:)
   type(i2dDPC),   allocatable::c_isb(:)
   type(r2dDPC),device,pointer::d_sbk(:)
   type(r2dDPC),   allocatable::c_sbk(:)
   type(iDPC)  ,device,pointer::d_nbstrbnd(:)
   type(iDPC)  ,   allocatable::c_nbstrbnd(:)
#endif
end
