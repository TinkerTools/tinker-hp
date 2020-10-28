c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kopbnd  --  forcefield parameters for out-of-plane bend  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxnopb   maximum number of out-of-plane bending entries
c
c     opbn      force constant parameters for out-of-plane bending
c     kopb      signature of atom classes for out-of-plane bending
c     kopb_sys  system signature of atom classes for out-of-plane bending
c
c
#include "tinker_precision.h"
      module kopbnd
      implicit none
      integer maxnopb
      parameter (maxnopb=500)
      real(t_p) opbn(maxnopb)
      integer(8) kopb(maxnopb)
      integer(8) kopb_sys(0:maxnopb)
      !DIR$ ATTRIBUTES ALIGN:64:: jopb
      logical, allocatable :: jopb(:)
!$acc declare create(kopb,kopb_sys)
      end
