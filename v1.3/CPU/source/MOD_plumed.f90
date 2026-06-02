!
!     ###################################################################
!     ##                                                               ##
!     ##  module plumed -- PLUMED/Tinker-HP interface                  ##
!     ##                                                               ##
!     ###################################################################
!
module plumed
   implicit none
#ifdef PLUMED
   integer :: ncount !<counter of calls to PLUMED
   integer, allocatable, dimension(:) :: pl_glob
   logical :: lplumed = .false. !<flag governing use of PLUMED
   real*8 :: energyUnits !<unit of energies
   real*8 :: lengthUnits !<unit of lengths
   real*8 :: timeUnits !<unit of time
   character(len=100) :: pl_input !<name of PLUMED input
   character(len=100) :: pl_output !<name of PLUMED output
   real*8, allocatable :: pl_pos(:,:) !<array of positions
   real*8, allocatable :: pl_force(:,:) !<array of forces
   real*8, allocatable :: pl_mass(:) !<array of masses
   real*8 :: pl_virial(3,3) !<virial contribution from PLUMED
   real*8 :: pl_epot !<energy contribution from PLUMED
#endif
   save
end
