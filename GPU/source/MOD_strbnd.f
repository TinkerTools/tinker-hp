c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module strbnd  --  stretch-bends in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     sbk       force constants for stretch-bend terms
c     winsbk    window object corresponding to sbk
c     nstrbnd   total number of stretch-bend interactions
c     nstrbndloc   local number of stretch-bend interactions
c     nbstrbnd   number of stretch-bend interactions before each atom
c     winnbstrbnd    window object corresponding to nbstrbnd
c     isb       angle and bond numbers used in stretch-bend
c     winisb    window object corresponding to isb
c
c
#include "tinker_precision.h"
      module strbnd
#ifdef USE_NVSHMEM_CUDA
      use tinTypes,only: i2dDPC=> Int2dDevPointerContainer
     &            ,        iDPC=>   IntDevPointerContainer
     &            ,      r2dDPC=>Real2dDevPointerContainer
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
