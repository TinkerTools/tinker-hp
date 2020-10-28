c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kiprop  --  forcefield parameters for improper dihedral  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxndi   maximum number of improper dihedral parameter entries
c
c     dcon     force constant parameters for improper dihedrals
c     tdi      ideal dihedral angle values for improper dihedrals
c     kdi      string of atom classes for improper dihedral angles
c
c
#include "tinker_precision.h"
      module kiprop
      implicit none
      integer maxndi
      parameter (maxndi=500)
      real(t_p) dcon(maxndi),tdi(maxndi)
      integer(8) kdi(maxndi)
      integer(8) kdi_sys(0:maxndi)
!$acc declare create(kdi_sys)
      end
