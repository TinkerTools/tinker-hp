c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module params  --  contents of force field parameter file  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nprm      number of nonblank lines in the parameter file
c     prmline   contents of each individual parameter file line
c
c
      module params
      use sizes
      implicit none
      integer nprm
      character*240 prmline(maxprm)
      save
      end
