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
c     kpt      string of atom classes for pi-orbital torsion terms
c
c
      module kpitor
      implicit none
      integer maxnpt
      parameter (maxnpt=500)
      real*8 ptcon(maxnpt)
      character*8 kpt(maxnpt)
      save
      end
