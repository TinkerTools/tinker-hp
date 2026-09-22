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
      module subroutine release_tmatxb_workspace(queue)
         integer,intent(in):: queue
      end subroutine
   end interface
end module
