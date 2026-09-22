!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module tmatxb  --  T-matrix workspace owner               ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module tmatxb
   implicit none
   real(t_p),allocatable,target :: tmatxb_workspace(:)
   interface
      module subroutine transposez_r6(src,dst,loc,sloc,siz,queue)
         real(t_p),intent(inout) :: src(*)
         real(t_p),intent(inout) :: dst(*)
         integer,intent(in) :: loc(*),sloc,siz,queue
      end subroutine
      module subroutine release_tmatxb_workspace(queue)
         integer,intent(in):: queue
      end subroutine
   end interface
end module
