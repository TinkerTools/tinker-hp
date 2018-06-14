c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  tortor.i  --  torsion-torsions in the current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ntortor   total number of torsion-torsion interactions
c     itt       atoms and parameter indices for torsion-torsion
c
c
      integer ntortor,ntortorloc
      integer, pointer :: itt(:,:),nbtortor(:)
      common /tortor/ ntortor,itt,ntortorloc,nbtortor
