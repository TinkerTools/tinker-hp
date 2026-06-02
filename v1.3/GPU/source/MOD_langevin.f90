!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ########################################################################
!     ##                                                                    ##
!     ##  module langevin  --  parameters and arrays for langevin dynamics  ##
!     ##                                                                    ##
!     ########################################################################
!
!     gamma : friction parameter in ps-1
!     Rn : white noise
!
!
#include "tinker_macro.h"
module langevin
   implicit none
   logical :: use_noselangevin
   logical :: use_noselangevin_massive
   real(r_p) gamma
   real(r_p), allocatable :: gamma_friction(:)
   real(r_p) nose,nose_mass
   !DIR$ ATTRIBUTES ALIGN:64 :: Rn
   real(t_p), allocatable :: Rn(:,:),noses(:,:)
end
