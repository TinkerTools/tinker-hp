c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  align.i  --  information for superposition of structures  ##
c     ##                                                            ##
c     ################################################################
c
c
c     wfit    weights assigned to atom pairs during superposition
c     nfit    number of atoms to use in superimposing two structures
c     ifit    atom numbers of pairs of atoms to be superimposed
c
c
c      integer nfit,ifit
c      real*8 wfit
c      common /align/ wfit(maxatm),nfit,ifit(2,maxatm)
