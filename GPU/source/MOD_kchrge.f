c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module kchrge  --  forcefield parameters for partial charges  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     chg   partial charge parameters for each atom type
c
c
#include "tinker_precision.h"
      module kchrge
      use sizes ,only: maxtyp
      implicit none
      real(t_p) chg(maxtyp)
      save
      end
