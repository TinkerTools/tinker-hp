!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  module angtor  --  angle-torsions in current structure  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     nangtor   total number of angle-torsion interactions
!     nangtorloc   local number of angle-torsion interactions
!     nbangtor   number of angle-torsion interactions before each atom
!     winnbangtor    window object corresponding to nbangtor
!     iat       torsion and angle numbers used in angle-torsion
!     winiat    window object corresponding to iat
!     kant      1-, 2- and 3-fold angle-torsion force constants
!     winkant   window object corresponding to kant
!
!
#include "tinker_macro.h"
module angtor
   implicit none
   integer nangtor,nangtorloc
   integer winiat,winnbangtor,winkant
   integer  ,pointer :: iat(:,:),nbangtor(:)
   real(t_p),pointer :: kant(:,:)
end
