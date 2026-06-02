!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module bound  --  control of periodic boundary conditions  ##
!     ##                                                             ##
!     #################################################################
!
!
module bound
   implicit none
   logical :: use_bounds !<flag to use periodic boundary conditions
   logical :: use_polymer !<flag to mark presence of infinite polymer
   save
end
