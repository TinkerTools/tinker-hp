c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kcpen  --  charge penetration forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     cpele     valence electron magnitude for each atom class
c     cpalp     alpha charge penetration parameter for each atom class
c
c
#include "tinker_precision.h"
      module kcpen
      use sizes ,only: maxclass
      implicit none
      real(t_p) cpele(maxclass)
      real(t_p) cpalp(maxclass)
      end
