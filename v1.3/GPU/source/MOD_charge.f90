!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module charge  --  partial charges for the current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     nion      total number of partial charges in system
!     nionloc   local number of partial charges in system
!     nionlocloop  First multiple of 16 after nionloc (help vecto)
!     nionbloc   local+neighbors number of partial charges in system
!     nionlocnl  localnl number of partial charges in system
!     nionlocnlloop  First multiple of 16 after nionlocnl (help vecto)
!     nionrecloc   local reciprocal number of partial charges in system
!     iion      number of the atom site for each partial charge
!     winiion    window object corresponding to iion
!     jion      neighbor generation site for each partial charge
!     winjion    window object corresponding to jion
!     kion      cutoff switching site for each partial charge
!     winkion    window object corresponding to kion
!     chglist   partial charge site for each atom (0=no charge)
!     winchglist    window object corresponding to chglist
!     nbchg     number of charges before each index
!     winnbchg    window object corresponding to nbchg
!     chgloc    global-local charge correspondance
!     chglocnl  global-localnl charge correspondance
!     chgrecloc  global-local reciprocal charge correspondance
!     pchg      magnitude of the partial charges (e-)
!     winpchg    window object corresponding to pchg
!     pchg_orig   original magnitude of the partial charges (e-) (lambda dyn)
!     winpchg_orig    window object corresponding to pchg_orig
!     pchg0     original partial charge values for charge flux
!     winpchg0    window object corresponding to pchg0
!
!     nionlocnlb First multiple of BLOCK_SIZE after nionlocnl
!     nionlocnlb_pair  total number of charge pair blocks interaction
!     nionlocnlb2_pair  total number of charge pair blocks interaction from C2 nblist
!     nshortionlocnlb2_pair  total number of charge pair blocks interaction in short range interaction list
!
#include "tinker_macro.h"
module charge
   implicit none
   integer nion,nionloc,nionbloc,nionlocnl,nionrecloc
   integer nionlocloop,nionlocnlloop
   integer nionlocnlb,nionlocnlb_pair,nionlocnlb2_pair&
      &,nshortionlocnlb2_pair
   integer, allocatable :: chgloc(:),chglocnl(:)
   integer, allocatable :: chgrecloc(:)
   integer, pointer :: iion(:)
   integer, pointer :: jion(:),kion(:)
   integer, pointer :: chglist(:)
   integer, pointer :: nbchg(:)
   integer winiion,winjion,winkion,winchglist,winnbchg&
      &,winpchg,winpchg0,winpchg_orig
   real(t_p), pointer :: pchg(:),pchg_orig(:),pchg0(:)
end
