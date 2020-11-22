c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kanang  --  forcefield parameters for angle-angle terms  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     anan   angle-angle cross term parameters for each atom class
c
c
      module kanang
      use sizes
      implicit none
      real*8 anan(3,maxclass)
      save
      end
