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
!     chg   partial charge parameters for each atom type
!
!
#include "tinker_macro.h"
module kchrge
   use sizes ,only: maxtyp
   implicit none
   real(t_p) chg(maxtyp)
   save
end
