!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module urey  --  Urey-Bradley interactions in the structure  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     uk      Urey-Bradley force constants (kcal/mole/Ang**2)
!     winuk    window object corresponding to uk
!     ul      ideal 1-3 distance values in Angstroms
!     winul    window object corresponding to ul
!     nurey   total number of Urey-Bradley terms in the system
!     nureyloc   local number of Urey-Bradley terms in the system
!     iury    numbers of the atoms in each Urey-Bradley interaction
!     winiury    window object corresponding to iury
!     nburey    numbers of Urey-Bradley interactions before each atom
!     winnburey    window object corresponding to nburey
!
!
#include "tinker_macro.h"
module urey
   implicit none
   integer nurey,nureyloc
   integer, pointer :: iury(:,:),nburey(:)
   real(t_p), pointer ::  uk(:),ul(:)
   integer :: winiury,winnburey,winuk,winul
end
