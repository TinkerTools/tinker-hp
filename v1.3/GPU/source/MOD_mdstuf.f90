!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module mdstuf  --  control of molecular dynamics trajectory  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     nfree       total number of degrees of freedom for a system
!     irest       steps between removal of COM inertia (0=no removal)
!     bmnmix      mixing coefficient for use with Beeman integrator
!     dorest      logical flag to remove center of mass inertia
!     velsave     logical flag to save velocity vector components
!     frcsave     logical flag to save force vector components
!     uindsave    logical flag to save induced atomic dipoles
!     integrate   type of molecular dynamics integration algorithm
!
!
#include "tinker_macro.h"
module mdstuf
   implicit none
   integer nfree,irest
   integer bmnmix
   logical dorest
   logical velsave
   logical frcsave
   logical uindsave
   logical mts
   character*20 integrate
   save
end
!
!     ###############################################################
!     #                                                             #
!     #  -- integrator workSpace module --                          #
!     #  holds data structure to be used inside integrator routines #
!     #                                                             #
!     ###############################################################
!
!     derivs  stores forces computes by gradient routine
!     etot    holds the system total energy at current timestep
!     epot    holds the system potential energy at current timestep
!     eksum   holds the system kinetic energy at current timestep
!     temp    holds the system temperature at current timestep
!     pres    holds the system pressure at current timestep
!
module mdstuf1
   implicit none
   logical  ,private:: isDataAlloc=.false.
   real(r_p),allocatable::derivs(:,:)
   real(r_p) etot,epot,eksum,ealt,ealt2,eml,ealtml
   real(r_p) temp,pres
   real(r_p) ekin(3,3),stress(3,3),viralt(3,3),viralt2(3,3)

contains
   subroutine gpuAllocMdstuf1Data
      if (.not.isDataAlloc) then
!$acc enter data create(etot,epot,eksum,ealt,ealt2,eml,ealtml,ekin &
!$acc                 ,temp,pres,stress,viralt,viralt2)
         isDataAlloc=.true.
      end if
   end subroutine
   subroutine gpuFreeMdstuf1Data
      if (isDataAlloc) then
!$acc exit data delete(etot,epot,eksum,ealt,ealt2,eml,ealtml,ekin &
!$acc                ,temp,pres,stress,viralt,viralt2)
         isDataAlloc=.false.
      end if
   end subroutine
end module
