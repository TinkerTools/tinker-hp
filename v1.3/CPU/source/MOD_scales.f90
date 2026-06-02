!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module scales  --  parameter scale factors for optimization  ##
!     ##                                                               ##
!     ###################################################################
!
!
module scales
   implicit none
   logical :: set_scale !<logical flag to show if scale factors have been set
   real*8, pointer :: scale(:) !<multiplicative factor for each optimization parameter
   save
end
