!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module chgpot  --  specifics of charge-charge functional form  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     electric   energy factor in kcal/mole for current force field
!     dielec     dielectric constant for electrostatic interactions
!     ebuffer    electrostatic buffering constant added to distance
!     c2scale    factor by which 1-2 charge interactions are scaled
!     c3scale    factor by which 1-3 charge interactions are scaled
!     c4scale    factor by which 1-4 charge interactions are scaled
!     c5scale    factor by which 1-5 charge interactions are scaled
!     neutnbr    logical flag governing use of neutral group neighbors
!     neutcut    logical flag governing use of neutral group cutoffs
!     ccorrect_ik      pair cscale interactions container
!     ccorrect_scale   vscale value of ccorrect_ik interaction
!     n_cscale         number of cscale interactions
!
!
#include "tinker_macro.h"
module chgpot
   implicit none
   real(t_p) electric
   real(t_p) dielec,ebuffer
   real(t_p) c2scale,c3scale
   real(t_p) c4scale,c5scale
   integer n_cscale
   integer  ,allocatable:: ccorrect_ik(:,:)
   real(t_p),allocatable:: ccorrect_scale(:)
   logical neutnbr,neutcut
end
