!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module repel  --  Pauli repulsion for current structure  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     nrep      total number of repulsion sites in the system
!     nreploc    local number of repulsion sites in the system
!     sizpr     Pauli repulsion size parameter value at each site
!     dmppr     Pauli repulsion alpha damping value at each site
!     elepr     Pauli repulsion valence electrons at each site
!     winsizepr window object corresponding to sizepr
!     windmppr window object corresponding to damppr
!     winelepr window object corresponding to elepr
!
!
#include "tinker_macro.h"
module repel
   implicit none
   integer nrep,nreploc
   real(t_p), pointer :: sizpr(:)
   real(t_p), pointer :: dmppr(:)
   real(t_p), pointer :: elepr(:)
   integer :: winsizpr,windmppr,winelepr
end
