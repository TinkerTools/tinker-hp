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
      real(t_p) gamma
      !DIR$ ATTRIBUTES ALIGN:64 :: Rn
      real(t_p), allocatable :: Rn(:,:)
      end
