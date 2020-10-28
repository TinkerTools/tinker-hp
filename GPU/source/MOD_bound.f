c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module bound  --  control of periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     use_bounds    flag to use periodic boundary conditions
c     use_polymer   flag to mark presence of infinite polymer
c
c
      module bound
      implicit none
      logical use_bounds
      logical use_polymer
      end
