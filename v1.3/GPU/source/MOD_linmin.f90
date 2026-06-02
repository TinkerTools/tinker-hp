!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module linmin  --  parameters for line search minimization  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     stpmin   minimum step length in current line search direction
!     stpmax   maximum step length in current line search direction
!     cappa    stringency of line search (0=tight < cappa < 1=loose)
!     slpmax   projected gradient above which stepsize is reduced
!     angmax   maximum angle between search direction and -gradient
!     intmax   maximum number of interpolations during line search
!
!
#include "tinker_macro.h"
module  linmin
   implicit none
   integer intmax
   real(r_p) stpmin,stpmax
   real(r_p) cappa,slpmax
   real(r_p) angmax
   save
end
