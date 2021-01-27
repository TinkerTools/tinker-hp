c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module ksttor  --  forcefield parameters for stretch-torsions  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxnbt   maximum number of stretch-torsion parameter entries
c
c     btcon    force constant parameters for stretch-torsion
c     kbt      string of atom classes for stretch-torsion terms
c
c
      module ksttor
      implicit none
      integer maxnbt
      parameter (maxnbt=500)
      real*8 btcon(9,maxnbt)
      character*16 kbt(maxnbt)
      save
      end
