!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module minima  --  general parameters for minimizations  ##
!     ##                                                           ##
!     ###############################################################
!
!
module minima
   implicit none
   integer :: maxiter !<maximum number of iterations during optimization
   integer :: nextiter !<iteration number to use for the first iteration
   real*8 :: fctmin !<value below which function is deemed optimized
   save
end
