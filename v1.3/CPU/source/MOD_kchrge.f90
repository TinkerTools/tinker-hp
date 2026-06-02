!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module kchrge  --  forcefield parameters for partial charges  ##
!     ##                                                                ##
!     ####################################################################
!
!
!
module kchrge
   use sizes
   implicit none
   real*8 :: chg(maxtyp) !<partial charge parameters for each atom type
   save
end
