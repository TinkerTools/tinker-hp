c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################################
c     ##                                                                    ##
c     ##  module langevin  --  parameters and arrays for langevin dynamics  ##
c     ##                                                                    ##
c     ########################################################################
c
c     gamma : friction parameter in ps-1
c     Rn : white noise 
c
c
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
