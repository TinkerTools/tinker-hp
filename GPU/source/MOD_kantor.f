c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  module kantor  --  angle-torsion forcefield parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxnat   maximum number of angle-torsion parameter entries
c
c     atcon    torsional amplitude parameters for angle-torsion
c     kat      id of atom classes for angle-torsion terms of the model
c     kat_sys  id of atom classes for angle-torsion terms of the system
c
c
#include "tinker_precision.h"
      module kantor
      implicit none
      integer maxnat
      parameter (maxnat=500)
      real(t_p) atcon(6,maxnat)
      integer(8) kat(maxnat)
      integer(8) kat_sys(0:maxnat)
      end
