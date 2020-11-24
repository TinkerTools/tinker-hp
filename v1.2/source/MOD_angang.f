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
c     winkaa    window object corresponding to kaa
c     nangang   total number of angle-angle interactions
c     nangangloc   local number of angle-angle interactions
c     nbangang   total number of angle-angle interactions before each atom in the global index
c     winnbangang    window object corresponding to nbangang
c     iaa       angle numbers used in each angle-angle term
c     winiaa    window object corresponding to iaa
c
c
      module angang
      implicit none
      integer nangang,nangangloc
      integer, pointer :: nbangang(:)
      integer, pointer :: iaa(:,:)
      real*8, pointer ::  kaa(:)
      integer winnbangang,winiaa,winkaa
      save
      end
