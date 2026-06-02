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
!     scale      multiplicative factor for each optimization parameter
!     set_scale  logical flag to show if scale factors have been set
!
!
#include "tinker_macro.h"
module scales
   implicit none
   real(r_p), pointer :: scale(:)
   logical set_scale
   save
end
