c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module angang  --  angle-angle terms in current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     kaa       force constant for angle-angle cross terms
c     nangang   total number of angle-angle interactions
c     nangangloc   local number of angle-angle interactions
c     nbangang   total number of angle-angle interactions before each atom in the global index
c     iaa       angle numbers used in each angle-angle term
c
c
      module angang
      implicit none
      integer nangang,nangangloc
      integer, pointer :: nbangang(:)
      integer, pointer :: iaa(:,:)
      real*8, pointer ::  kaa(:)
      save
      end
