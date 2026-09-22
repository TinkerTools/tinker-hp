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
