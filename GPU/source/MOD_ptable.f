c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module ptable  --  atomic symbols for the chemical elements  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     elemnt   atomic symbol for each chemical element
c
c
      module ptable
      use sizes
      implicit none
      character*3 elemnt(maxele)
      save
      end
