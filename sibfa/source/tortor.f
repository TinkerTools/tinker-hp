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
c     nbtortor   number of atoms before each torsion-torsion
c
c
      module tortor
      implicit none
      integer ntortor,ntortorloc
      integer, pointer :: itt(:,:),nbtortor(:)
      save
      end
