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
!
module mdstuf
   implicit none
   integer :: nfree !<total number of degrees of freedom for a system
   integer :: irest !<steps between removal of COM inertia (0=no removal)
   integer :: bmnmix !<mixing coefficient for use with Beeman integrator
   logical :: dorest !<logical flag to remove center of mass inertia
   logical :: velsave !<logical flag to save velocity vector components
   logical :: frcsave !<logical flag to save force vector components
   logical :: uindsave !<logical flag to save induced atomic dipoles
   logical :: mts !<logical flag regarding multi-timestep integration
   character*11 :: integrate !<type of molecular dynamics integration algorithm
   save
end
