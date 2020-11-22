c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kopbnd  --  forcefield parameters for out-of-plane bend  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxnopb   maximum number of out-of-plane bending entries
c
c     opbn      force constant parameters for out-of-plane bending
c     kopb      string of atom classes for out-of-plane bending
c
c
      module kopbnd
      implicit none
      integer maxnopb
      parameter (maxnopb=500)
      real*8 opbn(maxnopb)
      character*16 kopb(maxnopb)
      !DIR$ ATTRIBUTES ALIGN:64:: jopb
      logical, allocatable :: jopb(:)
      end
