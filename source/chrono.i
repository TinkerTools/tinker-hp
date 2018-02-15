c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  chrono.i  --  clock time values for the current program  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     twall   current processor wall clock time in seconds
c     tcpu    elapsed cpu time from start of program in seconds
c
c
      real*8 twall,tcpu
      common /chrono/ twall,tcpu
