c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module minima  --  general parameters for minimizations  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     fctmin    value below which function is deemed optimized
c     maxiter   maximum number of iterations during optimization
c     nextiter  iteration number to use for the first iteration
c
c
#include "tinker_precision.h"
      module minima
      implicit none
      integer maxiter,nextiter
      real(r_p) fctmin
      save
      end
