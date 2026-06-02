!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module krepl  --  Pauli repulsion forcefield parameters  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     prsiz     Pauli repulsion size value for each atom class
!     prdmp     alpha Pauli repulsion parameter for each atom class
!     prele     number of valence electrons for each atom class
!
!
#include "tinker_macro.h"
module krepl
   use sizes ,only: maxtyp
   implicit none
   real(t_p) prsiz(maxtyp),prdmp(maxtyp),prele(maxtyp)
end
