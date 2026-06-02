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
!     anan   angle-angle cross term parameters for each atom class
!
!
#include "tinker_macro.h"
module kanang
   use sizes ,only: maxclass
   implicit none
   real(t_p) anan(3,maxclass)
   save
end
