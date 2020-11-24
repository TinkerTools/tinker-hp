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
      module kiprop
      implicit none
      integer maxndi
      parameter (maxndi=500)
      real*8 dcon(maxndi),tdi(maxndi)
      character*16 kdi(maxndi)
      save
      end
