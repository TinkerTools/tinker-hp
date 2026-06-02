!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module kpolr  --  forcefield parameters for polarizability  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     polr   dipole polarizability parameters for each atom type
!     athl   Thole polarizability damping value for each atom type
!     pgrp   connected types in polarization group of each atom type
!     dthl   alternate Thole direct polarization damping values
!
!
#include "tinker_macro.h"
module kpolr
   use sizes ,only: maxtyp,maxvalue
   implicit none
   integer pgrp(maxvalue,maxtyp)
   real(t_p) polr(maxtyp),athl(maxtyp),dthl(maxtyp)
end
