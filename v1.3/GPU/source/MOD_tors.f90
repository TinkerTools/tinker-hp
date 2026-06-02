!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module tors  --  torsional angles within the current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     tors1   1-fold amplitude and phase for each torsional angle
!     wintors1    window object corresponding to tors1
!     tors2   2-fold amplitude and phase for each torsional angle
!     wintors2    window object corresponding to tors2
!     tors3   3-fold amplitude and phase for each torsional angle
!     wintors3    window object corresponding to tors3
!     tors4   4-fold amplitude and phase for each torsional angle
!     wintors4    window object corresponding to tors4
!     tors5   5-fold amplitude and phase for each torsional angle
!     wintors5    window object corresponding to tors5
!     tors6   6-fold amplitude and phase for each torsional angle
!     wintors6    window object corresponding to tors6
!     ntors   total number of torsional angles in the system
!     ntorsloc   local number of torsional angles in the system
!     nbtors   number of torsional angles before each atom
!     winnbtors    window object corresponding to nbtors
!     itors   numbers of the atoms in each torsional angle
!     winitors    window object corresponding to itors
!
!
#include "tinker_macro.h"
module tors
#ifdef USE_NVSHMEM_CUDA
   use tinTypes,only: i2dDPC=> Int2dDevPointerContainer&
      &,        iDPC=>   IntDevPointerContainer&
      &,      r2dDPC=>Real2dDevPointerContainer
#endif
   implicit none
   integer ntors,ntorsloc
   integer, pointer :: nbtors(:)
   integer, pointer :: itors(:,:)
   real(t_p), pointer :: tors1(:,:),tors2(:,:),tors3(:,:)
   real(t_p), pointer :: tors4(:,:),tors5(:,:),tors6(:,:)
   integer :: winnbtors,winitors,wintors1,wintors2
   integer :: wintors3,wintors4,wintors5,wintors6

#ifdef USE_NVSHMEM_CUDA
   integer ntors_pe
   type(iDPC  ),device,pointer::d_nbtors(:)
   type(iDPC  ),   allocatable::c_nbtors(:)
   type(i2dDPC),device,pointer::d_itors(:)
   type(i2dDPC),   allocatable::c_itors(:)
   type(r2dDPC),device,pointer::d_tors1(:),d_tors2(:),d_tors3(:)
   type(r2dDPC),   allocatable::c_tors1(:),c_tors2(:),c_tors3(:)
   type(r2dDPC),device,pointer::d_tors4(:),d_tors5(:),d_tors6(:)
   type(r2dDPC),   allocatable::c_tors4(:),c_tors5(:),c_tors6(:)
#endif
end
