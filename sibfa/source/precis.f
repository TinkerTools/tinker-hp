c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module precis  --  values of machine precision tolerances  ##
c     ##                                                             ##
c     #################################################################
c
c
c     tiny    the smallest positive floating point value
c     small   the smallest relative floating point spacing
c     huge    the largest relative floating point spacing
c
c
      module precis
      implicit none
      real*8 tiny,small,huge
      save
      end
