!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module opdist  --  out-of-plane distances in current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     nopdist   total number of out-of-plane distances in the system
!     nopdistloc   local number of out-of-plane distances in the system
!     iopd      numbers of the atoms in each out-of-plane distance
!     winiopd    window object corresponding to iopd
!     nbopdist number of angle used in out-of-plane distance before each atom
!     winnbopdist    window object corresponding to nbopdist
!     opdk      force constant values for out-of-plane distance
!     winopdk    window object corresponding to opdk
!
!
#include "tinker_macro.h"
module opdist
   implicit none
   integer nopdist,nopdistloc
   integer, pointer :: iopd(:,:),nbopdist(:)
   real(t_p), pointer ::  opdk(:)
   integer :: winiopd,winnbopdist,winopdk
   save
end
