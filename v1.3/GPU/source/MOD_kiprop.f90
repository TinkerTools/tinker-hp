!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kiprop  --  forcefield parameters for improper dihedral  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     maxndi   maximum number of improper dihedral parameter entries
!
!     dcon     force constant parameters for improper dihedrals
!     tdi      ideal dihedral angle values for improper dihedrals
!     kdi      string of atom classes for improper dihedral angles
!
!
#include "tinker_macro.h"
module kiprop
   implicit none
   integer maxndi
   parameter (maxndi=500)
   real(t_p) dcon(maxndi),tdi(maxndi)
   integer(8) kdi(maxndi)
   integer(8) kdi_sys(0:maxndi)
!$acc declare create(kdi_sys)
end
