c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module bitor  --  bitorsions within the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nbitor  total number of bitorsions in the system
c     nbitorloc  local number of bitorsions in the system
c     ibitor  numbers of the atoms in each bitorsion
c     nbbitors  numbers of bitorsions before each angle in the global index
c
c
      module bitor
      implicit none
      integer nbitor,nbitorloc
      integer, pointer :: ibitor(:,:), nbbitors(:)
      save
      end
