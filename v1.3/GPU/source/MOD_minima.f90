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
!     fctmin    value below which function is deemed optimized
!     maxiter   maximum number of iterations during optimization
!     nextiter  iteration number to use for the first iteration
!
!
#include "tinker_macro.h"
module minima
   implicit none
   integer maxiter,nextiter
   real(r_p) fctmin
   save
end
