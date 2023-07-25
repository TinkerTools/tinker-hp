c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module mdstuf  --  control of molecular dynamics trajectory  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     nfree       total number of degrees of freedom for a system
c     irest       steps between removal of COM inertia (0=no removal)
c     bmnmix      mixing coefficient for use with Beeman integrator
c     dorest      logical flag to remove center of mass inertia
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     integrate   type of molecular dynamics integration algorithm
c
c
#include "tinker_macro.h"
      module mdstuf
      implicit none
      integer nfree,irest
      integer bmnmix
      logical dorest
      logical velsave
      logical frcsave
      logical uindsave
      character*20 integrate
      save
      end
c
c     ###############################################################
c     #                                                             #
c     #  -- integrator workSpace module --                          #
c     #  holds data structure to be used inside integrator routines #
c     #                                                             #
c     ###############################################################
c
c     derivs  stores forces computes by gradient routine
c     etot    holds the system total energy at current timestep
c     epot    holds the system potential energy at current timestep
c     eksum   holds the system kinetic energy at current timestep
c     temp    holds the system temperature at current timestep
c     pres    holds the system pressure at current timestep
c
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
!$acc enter data create(etot,epot,eksum,ealt,ealt2,eml,ealtml,ekin
!$acc&                 ,temp,pres,stress,viralt,viralt2)
         isDataAlloc=.true.
      end if
      end subroutine
      subroutine gpuFreeMdstuf1Data
      if (isDataAlloc) then
!$acc exit data delete(etot,epot,eksum,ealt,ealt2,eml,ealtml,ekin
!$acc&                ,temp,pres,stress,viralt,viralt2)
         isDataAlloc=.false.
      end if
      end subroutine
      end module
