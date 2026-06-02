!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module moldyn  --  velocity and acceleration on MD trajectory  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module moldyn
   implicit none
   real*8, allocatable :: v(:,:) !<current velocity of each atom along the x,y,z-axes
   real*8, allocatable :: a(:,:) !<current acceleration of each atom along x,y,z-axes
   real*8, allocatable :: aalt(:,:) !<alternate acceleration of each atom along x,y,z-axes
   real*8, allocatable :: aalt2(:,:) !<alternate 2 acceleration of each atom along x,y,z-axes
   real*8 :: dshort !<short range (bonded) forces time step for respa and respa1 integrator
   real*8 :: dinter !<intermediate range forces time step for respa1 integrator
   integer :: nalt !<number of inner short range time steps for respa integrator and intermediate for respa1

   integer :: nalt2 !<number of inner short range time steps for respa1 integrator
   save
end
