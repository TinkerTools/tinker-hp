!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module strtor  --  stretch-torsions in the current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     kst       1-, 2- and 3-fold stretch-torsion force constants
!     winkst    window object corresponding to kst
!     nstrtor   total number of stretch-torsion interactions
!     nstrtorloc   local number of stretch-torsion interactions
!     nbstrtor   number of stretch-torsion interactions before each atom
!     winnbstrtor    window object corresponding to nbstrtor
!     ist       torsion and bond numbers used in stretch-torsion
!     winist    window object corresponding to ist
!
!
#include "tinker_macro.h"
module strtor
   implicit none
   integer nstrtor,nstrtorloc
   integer, pointer :: ist(:,:),nbstrtor(:)
   real(t_p), pointer :: kst(:,:)
   integer :: winist,winnbstrtor,winkst
end
