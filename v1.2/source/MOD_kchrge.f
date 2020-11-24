c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module kchrge  --  forcefield parameters for partial charges  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     chg   partial charge parameters for each atom type
c
c
      module kchrge
      use sizes
      implicit none
      real*8 chg(maxtyp)
      save
      end
