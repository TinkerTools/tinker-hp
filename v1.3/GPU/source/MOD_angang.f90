!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module angang  --  angle-angle terms in current structure  ##
!     ##                                                             ##
!     #################################################################
!
!
!     kaa       force constant for angle-angle cross terms
!     winkaa    window object corresponding to kaa
!     nangang   total number of angle-angle interactions
!     nangangloc   local number of angle-angle interactions
!     nbangang   total number of angle-angle interactions before each atom in the global index
!     winnbangang    window object corresponding to nbangang
!     iaa       angle numbers used in each angle-angle term
!     winiaa    window object corresponding to iaa
!
!
#include "tinker_macro.h"
module angang
   implicit none
   integer nangang,nangangloc
   integer, pointer :: nbangang(:)
   integer, pointer :: iaa(:,:)
   real(t_p), pointer ::  kaa(:)
   integer winnbangang,winiaa,winkaa
end
