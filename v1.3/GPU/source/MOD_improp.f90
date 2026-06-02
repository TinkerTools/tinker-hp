!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module improp  --  improper dihedrals in the current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     kprop    force constant values for improper dihedral angles
!     winkprop    window object corresponding to kprop
!     vprop    ideal improper dihedral angle value in degrees
!     winvprop    window object corresponding to vprop
!     niprop   total number of improper dihedral angles in the system
!     iiprop   numbers of the atoms in each improper dihedral angle
!     winiiprop    window object corresponding to iiprop
!     nbimprop number of improper diehedral before each atom
!     winnbimprop    window object corresponding to nbimprop
!
!
#include "tinker_macro.h"
module improp
   implicit none
   integer niprop,niproploc
   integer, pointer :: iiprop(:,:),nbimprop(:)
   real(t_p), pointer :: kprop(:),vprop(:)
   integer :: winiiprop,winnbimprop,winkprop,winvprop
   save
end
