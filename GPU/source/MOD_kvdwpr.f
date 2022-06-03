c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kvdwpr  --  forcefield parameters for special vdw terms  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxnvp   maximum number of special van der Waals pair entries
c     vdwpr_l  switch to detect the use of van der Waals special pairs
c
c     radpr    radius parameter for special van der Waals pairs
c     epspr    well depth parameter for special van der Waals pairs
c     kvpr     string of atom classes for special van der Waals pairs
c
c
#include "tinker_precision.h"
      module kvdwpr
      implicit none
      integer maxnvp
      parameter (maxnvp=500)
      logical vdwpr_l
      real(t_p) radpr(maxnvp),epspr(maxnvp)
      character*8 kvpr(maxnvp)
      save
      end
