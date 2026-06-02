!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module chgpen  --  charge penetration in current structure  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     ncp       total number of charge penetration sites in system
!     pcore     number of core electrons at each multipole site
!     pval      number of valence electrons at each multipole site
!     pval0     original number of valence electrons at each multipole site for chglx
!     palpha    charge penetration damping at each multipole site
!     winpcore  window corresponding to pcore array
!     winpval   window corresponding to pval array
!     winpalpha window corresponding to palpha array
!
!
#include "tinker_macro.h"
module chgpen
   implicit none
   integer ncp
   real(t_p), pointer :: pcore(:)
   real(t_p), pointer :: pval(:)
   real(t_p), pointer :: pval0(:)
   real(t_p), pointer :: palpha(:)
   integer :: winpcore,winpval,winpval0,winpalpha
end
