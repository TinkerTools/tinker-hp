!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  submodule polar_lifecycle  --  polarization workspace       ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
submodule(polar) polar_lifecycle
   use lifecycle_observer, only: LifecyclePolar, record_lifecycle_event
contains
   module subroutine release_polar_workspace(queue)
      use tinMemory, only: s_prmem,sd_prmem,szoTp
      implicit none
      integer,intent(in):: queue

      if (.not.allocated(polar_workspace)) return
      call record_lifecycle_event(LifecyclePolar)
!$acc wait(queue)
!$acc exit data delete(polar_workspace)
      s_prmem = s_prmem - size(polar_workspace)*szoTp
      sd_prmem = sd_prmem - size(polar_workspace)*szoTp
      deallocate(polar_workspace)
   end subroutine
end submodule
