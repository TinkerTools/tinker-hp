c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kctrn  --  charge transfer forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ctchg     charge transfer magnitude for each atom class
c     ctdmp     alpha charge transfer parameter for each atom class
c
c
#include "tinker_precision.h"
      module kctrn
      use sizes ,only: maxtyp
      implicit none
      real(t_p) ctchg(maxtyp)
      real(t_p) ctdmp(maxtyp)
      end
