c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module tors  --  torsional angles within the current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     tors1   1-fold amplitude and phase for each torsional angle
c     wintors1    window object corresponding to tors1
c     tors2   2-fold amplitude and phase for each torsional angle
c     wintors2    window object corresponding to tors2
c     tors3   3-fold amplitude and phase for each torsional angle
c     wintors3    window object corresponding to tors3
c     tors4   4-fold amplitude and phase for each torsional angle
c     wintors4    window object corresponding to tors4
c     tors5   5-fold amplitude and phase for each torsional angle
c     wintors5    window object corresponding to tors5
c     tors6   6-fold amplitude and phase for each torsional angle
c     wintors6    window object corresponding to tors6
c     ntors   total number of torsional angles in the system
c     ntorsloc   local number of torsional angles in the system
c     nbtors   number of torsional angles before each atom
c     winnbtors    window object corresponding to nbtors
c     itors   numbers of the atoms in each torsional angle
c     winitors    window object corresponding to itors
c
c
#include "tinker_precision.h"
      module tors
#ifdef USE_NVSHMEM_CUDA
      use tinTypes,only: i2dDPC=> Int2dDevPointerContainer 
     &            ,        iDPC=>   IntDevPointerContainer 
     &            ,      r2dDPC=>Real2dDevPointerContainer
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
