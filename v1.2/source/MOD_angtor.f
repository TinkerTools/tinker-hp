c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  module angtor  --  angle-torsions in current structure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nangtor   total number of angle-torsion interactions
c     nangtorloc   local number of angle-torsion interactions
c     nbangtor   number of angle-torsion interactions before each atom
c     winnbangtor    window object corresponding to nbangtor
c     iat       torsion and angle numbers used in angle-torsion
c     winiat    window object corresponding to iat
c     kant      1-, 2- and 3-fold angle-torsion force constants
c     winkant   window object corresponding to kant
c
c
      module angtor
      implicit none
      integer nangtor,nangtorloc
      integer, pointer :: iat(:,:),nbangtor(:)
      real*8, pointer :: kant(:,:)
      integer :: winiat,winnbangtor,winkant
      save
      end
