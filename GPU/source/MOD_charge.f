c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module charge  --  partial charges for the current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     nion      total number of partial charges in system
c     nionloc   local number of partial charges in system
c     nionlocloop  First multiple of 16 after nionloc (help vecto)
c     nionbloc   local+neighbors number of partial charges in system
c     nionlocnl  localnl number of partial charges in system
c     nionlocnlloop  First multiple of 16 after nionlocnl (help vecto)
c     nionrecloc   local reciprocal number of partial charges in system
c     iion      number of the atom site for each partial charge
c     winiion    window object corresponding to iion
c     jion      neighbor generation site for each partial charge
c     winjion    window object corresponding to jion
c     kion      cutoff switching site for each partial charge
c     winkion    window object corresponding to kion
c     chglist   partial charge site for each atom (0=no charge)
c     winchglist    window object corresponding to chglist
c     nbchg     number of charges before each index
c     winnbchg    window object corresponding to nbchg
c     chgloc    global-local charge correspondance
c     chglocnl  global-localnl charge correspondance
c     chgrecloc  global-local reciprocal charge correspondance
c     pchg      magnitude of the partial charges (e-)
c     winpchg    window object corresponding to pchg
c     pchg_orig   original magnitude of the partial charges (e-) (lambda dyn)
c     winpchg_orig    window object corresponding to pchg_orig
c     pchg0     original partial charge values for charge flux
c     winpchg0    window object corresponding to pchg0
c
c     nionlocnlb First multiple of BLOCK_SIZE after nionlocnl
c     nionlocnlb_pair  total number of charge pair blocks interaction
c     nionlocnlb2_pair  total number of charge pair blocks interaction from C2 nblist
c     nshortionlocnlb2_pair  total number of charge pair blocks interaction in short range interaction list
c
#include "tinker_precision.h"
      module charge
      implicit none
      integer nion,nionloc,nionbloc,nionlocnl,nionrecloc
      integer nionlocloop,nionlocnlloop
      integer nionlocnlb,nionlocnlb_pair,nionlocnlb2_pair
     &       ,nshortionlocnlb2_pair
      integer, allocatable :: chgloc(:),chglocnl(:)
      integer, allocatable :: chgrecloc(:)
      integer, pointer :: iion(:)
      integer, pointer :: jion(:),kion(:)
      integer, pointer :: chglist(:)
      integer, pointer :: nbchg(:)
      integer winiion,winjion,winkion,winchglist,winnbchg
     &       ,winpchg,winpchg0,winpchg_orig
      real(t_p), pointer :: pchg(:),pchg_orig(:),pchg0(:)
      end
