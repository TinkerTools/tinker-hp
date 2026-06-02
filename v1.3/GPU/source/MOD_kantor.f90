!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  module kantor  --  angle-torsion forcefield parameters  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     maxnat   maximum number of angle-torsion parameter entries
!
!     atcon    torsional amplitude parameters for angle-torsion
!     kat      id of atom classes for angle-torsion terms of the model
!     kat_sys  id of atom classes for angle-torsion terms of the system
!
!
#include "tinker_macro.h"
module kantor
   implicit none
   integer maxnat
   parameter (maxnat=500)
   real(t_p) atcon(6,maxnat)
   integer(8) kat(maxnat)
   integer(8) kat_sys(0:maxnat)
end
