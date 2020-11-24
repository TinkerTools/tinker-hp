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
c     ku      string of atom classes for Urey-Bradley terms
c
c
      module kurybr
      implicit none
      integer maxnu
      parameter (maxnu=2000)
      real*8 ucon(maxnu),dst13(maxnu)
      character*12 ku(maxnu)
      save
      end
