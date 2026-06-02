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
!     nring3   total number of 3-membered rings in the system
!     iring3   numbers of the atoms involved in each 3-ring
!     nring4   total number of 4-membered rings in the system
!     iring4   numbers of the atoms involved in each 4-ring
!     nring5   total number of 5-membered rings in the system
!     iring5   numbers of the atoms involved in each 5-ring
!     nring6   total number of 6-membered rings in the system
!     iring6   numbers of the atoms involved in each 6-ring
!
!
module ring
   use sizes
   implicit none
   integer nring3,iring3(3,maxring)
   integer nring4,iring4(4,maxring)
   integer nring5,iring5(5,maxring)
   integer nring6,iring6(6,maxring)
   save
end
