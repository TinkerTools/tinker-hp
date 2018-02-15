c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  keys.i  --  contents of current keyword parameter file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c
      integer nkey
      character*120 keyline
      common /keys/ nkey,keyline(maxkey)
