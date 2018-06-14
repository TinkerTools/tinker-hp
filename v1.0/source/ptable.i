c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  ptable.i  --  atomic symbols for the chemical elements  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     elemnt   atomic symbol for each chemical element
c
c
      character*3 elemnt
      common /ptable/ elemnt(maxele)
