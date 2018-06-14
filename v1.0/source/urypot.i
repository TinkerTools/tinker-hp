c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  urypot.i  --  specifics of Urey-Bradley functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     cury       cubic coefficient in Urey-Bradley potential
c     qury       quartic coefficient in Urey-Bradley potential
c     ureyunit   convert Urey-Bradley energy to kcal/mole
c
c
      real*8 cury,qury
      real*8 ureyunit
      common /urypot/ cury,qury,ureyunit
