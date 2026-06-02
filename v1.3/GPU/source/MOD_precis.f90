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
!     tiny    the smallest positive floating point value
!     small   the smallest relative floating point spacing
!     huge    the largest relative floating point spacing
!
!
#include "tinker_macro.h"
module precis
   implicit none
   real(t_p) tiny,small,huge
   save
end
