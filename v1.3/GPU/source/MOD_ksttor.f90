!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module ksttor  --  forcefield parameters for stretch-torsions  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     maxnbt   maximum number of stretch-torsion parameter entries
!
!     btcon    force constant parameters for stretch-torsion
!     kbt      string of atom classes for stretch-torsion terms
!
!
#include "tinker_macro.h"
module ksttor
   implicit none
   integer maxnbt
   parameter (maxnbt=500)
   real(t_p) btcon(9,maxnbt)
   integer(8) kbt(maxnbt)
   integer(8) kbt_sys(0:maxnbt)
   save
!$acc declare create(kbt,kbt_sys)
end
