c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kitors  --  forcefield parameters for improper torsions  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxnti   maximum number of improper torsion parameter entries
c
c     ti1      torsional parameters for improper 1-fold rotation
c     ti2      torsional parameters for improper 2-fold rotation
c     ti3      torsional parameters for improper 3-fold rotation
c     kti      string of atom classes for improper torsional parameters
c
c
#include "tinker_precision.h"
      module kitors
      implicit none
      integer maxnti
      parameter (maxnti=500)
      real(t_p) ti1(2,maxnti),ti2(2,maxnti),ti3(2,maxnti)
      character*16 kti(maxnti)
      save
      end
