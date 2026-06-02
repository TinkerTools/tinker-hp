!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kvdwpr  --  forcefield parameters for special vdw terms  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     maxnvp   maximum number of special van der Waals pair entries
!     vdwpr_l  switch to detect the use of van der Waals special pairs
!
!     radpr    radius parameter for special van der Waals pairs
!     epspr    well depth parameter for special van der Waals pairs
!     kvpr     string of atom classes for special van der Waals pairs
!
!
#include "tinker_macro.h"
module kvdwpr
   implicit none
   integer maxnvp
   parameter (maxnvp=500)
   logical vdwpr_l
   real(t_p) radpr(maxnvp),epspr(maxnvp)
   character*8 kvpr(maxnvp)
   save
end
