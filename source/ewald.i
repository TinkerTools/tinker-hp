c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  ewald.i  --  parameters and options for Ewald summation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     aewald     Ewald convergence coefficient value (Ang-1)
c     boundary   Ewald boundary condition; none, tinfoil or vacuum
c
c
      real*8 aewald
      character*7 boundary
      common /ewald/ aewald,boundary
