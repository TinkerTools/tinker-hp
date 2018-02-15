c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  ksttor.i  --  forcefield parameters for stretch-torsions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnbt   maximum number of stretch-torsion parameter entries
c
c     btcon    force constant parameters for stretch-torsion
c     kbt      string of atom classes for stretch-torsion terms
c
c
      integer maxnbt
      parameter (maxnbt=500)
      real*8 btcon
      character*16 kbt
      common /ksttor/ btcon(3,maxnbt),kbt(maxnbt)
