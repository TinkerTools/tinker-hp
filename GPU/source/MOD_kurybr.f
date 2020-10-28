c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  module kurybr  --  forcefield parameters for Urey-Bradley terms  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     maxnu   maximum number of Urey-Bradley parameter entries
c
c     ucon    force constant parameters for Urey-Bradley terms
c     dst13   ideal 1-3 distance parameters for Urey-Bradley terms
c     ku      integer signature of atom classes for Urey-Bradley terms
c     ku_sys  integer signature of atom classes for Urey-Bradley terms of the simulated system
c
c
#include "tinker_precision.h"
      module kurybr
      implicit none
      integer maxnu
      parameter (maxnu=2000)
      real(t_p) ucon(maxnu),dst13(maxnu)
      integer(8) ku(maxnu)
      integer(8) ku_sys(0:maxnu)
      save
!$acc declare create(ku,ku_sys)
      end
