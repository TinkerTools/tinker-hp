!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module mplpot  --  specifics of atomic multipole functions  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     m2scale   scale factor for 1-2 multipole energy interactions
!     m3scale   scale factor for 1-3 multipole energy interactions
!     m4scale   scale factor for 1-4 multipole energy interactions
!     m5scale   scale factor for 1-5 multipole energy interactions
!     pentyp       type of penetration damping (NONE, GORDON1, GORDON2)
!     mcorrect_ik      pair mscale interactions container
!     mcorrect_scale   mscale value of mcorrect_ik interaction
!     n_mscale         number of mscale interactions
!
!
#include "tinker_macro.h"
module mplpot
   implicit none
   real(t_p) m2scale,m3scale
   real(t_p) m4scale,m5scale
   character*7 pentyp
   integer pentyp_i
   integer n_mscale
   integer  ,allocatable:: mcorrect_ik(:,:)
   real(t_p),allocatable:: mcorrect_scale(:)
   enum,bind(C)
      enumerator PT_GORDON1,PT_GORDON2,PT_NONE
   end enum
!$acc declare create(m2scale,m3scale,m4scale,m5scale)
end
