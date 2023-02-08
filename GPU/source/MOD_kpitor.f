c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kpitor  --  forcefield parameters for pi-orbit torsions  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxnpt   maximum number of pi-orbital torsion parameter entries
c
c     ptcon    force constant parameters for pi-orbital torsions
c     kpt      integer siganture of atom classes for pi-orbital torsion terms
c     kpt_sys  system siganture of atom classes for pi-orbital torsion terms
c
c
#include "tinker_macro.h"
      module kpitor
      implicit none
      integer maxnpt
      parameter (maxnpt=500)
      real(t_p) ptcon(maxnpt)
      integer(8) kpt(maxnpt)
      integer(8) kpt_sys(0:maxnpt)
      save
!$acc declare create(kpt,kpt_sys)
      end
