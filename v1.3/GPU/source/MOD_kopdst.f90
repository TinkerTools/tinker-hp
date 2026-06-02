!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #######################################################################
!     ##                                                                   ##
!     ##  module kopdst  --  forcefield parameters for out-plane distance  ##
!     ##                                                                   ##
!     #######################################################################
!
!
!     maxnopd   maximum number of out-of-plane distance entries
!
!     opds      force constant parameters for out-of-plane distance
!     kopd      string of atom classes for out-of-plane distance
!
!
#include "tinker_macro.h"
module kopdst
   implicit none
   integer maxnopd
   parameter (maxnopd=500)
   real(t_p) opds(maxnopd)
   character*16 kopd(maxnopd)
   save
end
