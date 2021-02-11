c
c     ###################################################################
c     ##                                                               ##
c     ##  module plumed                                                ##
c     ##                                                               ##
c     ###################################################################
c
      module plumed
      implicit none
#ifdef PLUMED
      logical :: lplumed = .false.
      real*8 :: energyUnits
      real*8 :: lengthUnits
      real*8 :: timeUnits
      character(len=100) :: pl_input
      character(len=100) :: pl_output
      real*8, allocatable ::  pl_pos(:,:), pl_force(:,:), pl_mass(:)
      real*8 :: pl_virial(3,3), pl_epot
      integer :: ncount
      integer, allocatable, dimension(:) :: pl_glob
#endif
      save
      end
