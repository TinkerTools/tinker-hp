!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kanang  --  forcefield parameters for angle-angle terms  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module kanang
   use sizes
   implicit none
   real*8 :: anan(3,maxclass) !<angle-angle cross term parameters for each atom class
   save
end
