c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kanang  --  forcefield parameters for angle-angle terms  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     anan   angle-angle cross term parameters for each atom class
c
c
#include "tinker_precision.h"
      module kanang
      use sizes ,only: maxclass
      implicit none
      real(t_p) anan(3,maxclass)
      save
      end
