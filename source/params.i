c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  params.i  --  contents of force field parameter file  ##
c     ##                                                        ##
c     ############################################################
c
c
c     nprm      number of nonblank lines in the parameter file
c     prmline   contents of each individual parameter file line
c
c
      integer nprm
      character*120 prmline
      common /params/ nprm,prmline(maxprm)
