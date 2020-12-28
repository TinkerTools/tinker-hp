c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module keys  --  contents of current keyword parameter file  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c
      module keys
      use sizes
      implicit none
      integer nkey
      character*240 keyline(maxkey)
      save
      end
