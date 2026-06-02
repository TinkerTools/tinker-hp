!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################
!     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
!     ##                   All Rights Reserved                  ##
!     ############################################################
!
!     ##################################################################
!     ##                                                              ##
!     ##  module kcpen  --  charge penetration forcefield parameters  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     cpele     valence electron magnitude for each atom class
!     cpalp     alpha charge penetration parameter for each atom class
!
!
#include "tinker_macro.h"
module kcpen
   use sizes ,only: maxclass
   implicit none
   real(t_p) cpele(maxclass)
   real(t_p) cpalp(maxclass)
end
