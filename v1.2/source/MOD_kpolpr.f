c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  module kpolpr  --  special Thole forcefield parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c     thlpr    Thole damping values for special polarization pairs
c     thdpr    Thole direct damping for special polarization pairs
c     kppr     string of atom types for special polarization pairs
c
      module kpolpr
      use sizes
      implicit none
      integer maxnpp
      real*8 :: thlpr(maxtyp)
      real*8 :: thdpr(maxtyp)
      parameter (maxnpp=100)
      character*8 :: kppr(maxnpp)
      save
      end
