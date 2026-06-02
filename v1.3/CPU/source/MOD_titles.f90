!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module titles  --  title for the current molecular system  ##
!     ##                                                             ##
!     #################################################################
!
!
!
module titles
   implicit none
   integer :: ltitle !<length in characters of the nonblank title string
   character*240 :: title !<title used to describe the current structure
   save
end
