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
c     kat      string of atom classes for angle-torsion terms
c
c
      module kantor
      implicit none
      integer maxnat
      parameter (maxnat=500)
      real*8 atcon(6,maxnat)
      character*16 kat(maxnat)
      save
      end
