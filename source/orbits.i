c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  orbits.i  --  orbital energies for conjugated pisystem  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     q       number of pi-electrons contributed by each atom
c     w       ionization potential of each pisystem atom
c     em      repulsion integral for each pisystem atom
c
c
      real*8 q,w,em
      common /orbits/ q(maxatm),w(maxatm),em(maxatm)
