c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  module kopdst  --  forcefield parameters for out-plane distance  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     maxnopd   maximum number of out-of-plane distance entries
c
c     opds      force constant parameters for out-of-plane distance
c     kopd      string of atom classes for out-of-plane distance
c
c
#include "tinker_macro.h"
      module kopdst
      implicit none
      integer maxnopd
      parameter (maxnopd=500)
      real(t_p) opds(maxnopd)
      character*16 kopd(maxnopd)
      save
      end
