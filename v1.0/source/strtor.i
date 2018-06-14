c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  strtor.i  --  stretch-torsions in the current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     kst       1-, 2- and 3-fold stretch-torsion force constants
c     nstrtor   total number of stretch-torsion interactions
c     ist       torsion and bond numbers used in stretch-torsion
c
c
      integer nstrtor,nstrtorloc
      integer, pointer :: ist(:,:),nbstrtor(:)
      real*8, pointer :: kst(:,:)
      common /strtor/ kst,nstrtor,ist,nstrtorloc,nbstrtor
