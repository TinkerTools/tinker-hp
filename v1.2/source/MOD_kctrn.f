c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kctrn  --  charge transfer forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ctchg     charge transfer magnitude for each atom class
c     ctdmp     alpha charge transfer parameter for each atom class
c
c
      module kctrn
      use sizes
      implicit none
      real*8 :: ctchg(maxtyp)
      real*8 :: ctdmp(maxtyp)
      save
      end
