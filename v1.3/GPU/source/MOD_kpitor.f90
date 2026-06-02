!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kpitor  --  forcefield parameters for pi-orbit torsions  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     maxnpt   maximum number of pi-orbital torsion parameter entries
!
!     ptcon    force constant parameters for pi-orbital torsions
!     kpt      integer siganture of atom classes for pi-orbital torsion terms
!     kpt_sys  system siganture of atom classes for pi-orbital torsion terms
!
!
#include "tinker_macro.h"
module kpitor
   implicit none
   integer maxnpt
   parameter (maxnpt=500)
   real(t_p) ptcon(maxnpt)
   integer(8) kpt(maxnpt)
   integer(8) kpt_sys(0:maxnpt)
   save
!$acc declare create(kpt,kpt_sys)
end
