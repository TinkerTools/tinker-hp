!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module kctrn  --  charge transfer forcefield parameters  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     ctchg     charge transfer magnitude for each atom class
!     ctdmp     alpha charge transfer parameter for each atom class
!
!
#include "tinker_macro.h"
module kctrn
   use sizes ,only: maxtyp
   implicit none
   real(t_p) ctchg(maxtyp)
   real(t_p) ctdmp(maxtyp)
end
