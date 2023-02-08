c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module linmin  --  parameters for line search minimization  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     stpmin   minimum step length in current line search direction
c     stpmax   maximum step length in current line search direction
c     cappa    stringency of line search (0=tight < cappa < 1=loose)
c     slpmax   projected gradient above which stepsize is reduced
c     angmax   maximum angle between search direction and -gradient
c     intmax   maximum number of interpolations during line search
c
c
#include "tinker_macro.h"
      module  linmin
      implicit none
      integer intmax
      real(r_p) stpmin,stpmax
      real(r_p) cappa,slpmax
      real(r_p) angmax
      save 
      end
