!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module ptable  --  atomic symbols for the chemical elements  ##
!     ##                                                               ##
!     ###################################################################
!
!
!
module ptable
   use sizes
   implicit none
   character*3 :: elemnt(maxele) !<atomic symbol for each chemical element
   save
end
