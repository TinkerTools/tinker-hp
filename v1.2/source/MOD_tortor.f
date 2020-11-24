c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module tortor  --  torsion-torsions in the current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     ntortor   total number of torsion-torsion interactions
c     ntortorloc   local number of torsion-torsion interactions
c     itt       atoms and parameter indices for torsion-torsion
c     winitt    window object corresponding to itt
c     nbtortor   number of atoms before each torsion-torsion
c     winnbtortor    window object corresponding to nbtortor
c
c
      module tortor
      implicit none
      integer ntortor,ntortorloc
      integer, pointer :: itt(:,:),nbtortor(:)
      integer :: winitt,winnbtortor
      save
      end
