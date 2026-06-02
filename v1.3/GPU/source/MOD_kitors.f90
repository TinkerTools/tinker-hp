!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kitors  --  forcefield parameters for improper torsions  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     maxnti   maximum number of improper torsion parameter entries
!
!     ti1      torsional parameters for improper 1-fold rotation
!     ti2      torsional parameters for improper 2-fold rotation
!     ti3      torsional parameters for improper 3-fold rotation
!     kti      string of atom classes for improper torsional parameters
!
!
#include "tinker_macro.h"
module kitors
   implicit none
   integer maxnti
   parameter (maxnti=500)
   real(t_p) ti1(2,maxnti),ti2(2,maxnti),ti3(2,maxnti)
   integer(8) kti(maxnti)
   integer(8) kti_sys(0:maxnti)
!$acc declare create(kti_sys)
end
