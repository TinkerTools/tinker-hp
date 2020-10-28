c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module scales  --  parameter scale factors for optimization  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     scale      multiplicative factor for each optimization parameter
c     set_scale  logical flag to show if scale factors have been set
c
c
#include "tinker_precision.h"
      module scales
      implicit none
      real(r_p), pointer :: scale(:)
      logical set_scale
      save
      end
