!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module ring  --  number and location of small ring structures  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!
module ring
   use sizes
   implicit none
   integer :: nring3 !<total number of 3-membered rings in the system
   integer :: iring3(3,maxring) !<numbers of the atoms involved in each 3-ring
   integer :: nring4 !<total number of 4-membered rings in the system
   integer :: iring4(4,maxring) !<numbers of the atoms involved in each 4-ring
   integer :: nring5 !<total number of 5-membered rings in the system
   integer :: iring5(5,maxring) !<numbers of the atoms involved in each 5-ring
   integer :: nring6 !<total number of 6-membered rings in the system
   integer :: iring6(6,maxring) !<numbers of the atoms involved in each 6-ring
   save
end
