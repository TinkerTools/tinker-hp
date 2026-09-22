!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  submodule tmatxb_lifecycle  --  T-matrix workspace          ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
submodule(tmatxb) tmatxb_lifecycle
contains
   module subroutine transposez_r6(src,dst,loc,sloc,siz,queue)
      implicit none
      real(t_p),intent(inout) :: src(siz*6)
      real(t_p),intent(inout) :: dst(6*siz)
      integer,intent(in) :: loc(sloc),sloc,siz,queue
      integer :: i,j,idxd,idxs

!$acc host_data use_device(src,dst,loc)
!$acc parallel loop collapse(2) deviceptr(src,dst,loc) async(queue)
      do i = 1,sloc; do j = 0,5
         idxd = (loc(i)-1)*6+j+1
         idxs = j*siz+i
         dst(idxd) = dst(idxd) + src(idxs)
         src(idxs) = 0
      end do; end do
!$acc end host_data
   end subroutine

   module subroutine release_tmatxb_workspace(queue)
      use tinMemory, only: s_prmem,sd_prmem,szoTp
      implicit none
      integer,intent(in):: queue

      if (.not.allocated(tmatxb_workspace)) return
!$acc wait(queue)
!$acc exit data delete(tmatxb_workspace)
      s_prmem = s_prmem - size(tmatxb_workspace)*szoTp
      sd_prmem = sd_prmem - size(tmatxb_workspace)*szoTp
      deallocate(tmatxb_workspace)
   end subroutine
end submodule
