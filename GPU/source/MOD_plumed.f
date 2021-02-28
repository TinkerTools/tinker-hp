c
c     ###################################################################
c     ##                                                               ##
c     ##  module plumed                                                ##
c     ##                                                               ##
c     ###################################################################
c
#include "tinker_precision.h"
      module plumed
      implicit none
      integer  ncount
      logical:: lplumed=.false.
      real(r_p) energyUnits
      real(r_p) lengthUnits
      real(r_p) timeUnits
      character(len=100):: pl_input
      character(len=100):: pl_output

      integer  ,allocatable:: pl_glob(:)
      real(r_p),allocatable:: pl_pos(:,:),pl_force(:,:),pl_mass(:)
      real(r_p) :: pl_virial(3,3), pl_epot

      interface
         module subroutine plumed_init(dt)
         real(r_p),intent(in):: dt
         end subroutine
         module subroutine eplumed(energy,derivs)
         implicit none
         real(r_p) ,intent(in)   :: energy
         real(r_p) ,intent(inout):: derivs(:,:)
         end subroutine
      end interface

      end
