!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module precis  --  values of machine precision tolerances  ##
!     ##                                                             ##
!     #################################################################
!
!
!
module precis
   implicit none
   real*8 :: tiny !<the smallest positive floating point value
   real*8 :: small !<the smallest relative floating point spacing
   real*8 :: huge !<the largest relative floating point spacing
   save
end
