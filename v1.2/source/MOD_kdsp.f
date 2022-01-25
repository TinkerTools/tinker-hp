c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module kdsp  --  damped dispersion forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     dspsix   C6 dispersion coefficient for each atom class
c     dspdmp   alpha dispersion parameter for each atom class
c
c
      module kdsp
      use sizes
      implicit none
      real*8 dspsix(maxtyp)
      real*8 dspdmp(maxtyp)
      save
      end
