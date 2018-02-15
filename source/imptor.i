c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  imptor.i  --  improper torsions in the current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     itors1   1-fold amplitude and phase for each improper torsion
c     itors2   2-fold amplitude and phase for each improper torsion
c     itors3   3-fold amplitude and phase for each improper torsion
c     nitors   total number of improper torsional angles in the system
c     iitors   numbers of the atoms in each improper torsional angle
c
c
      integer nitors,nitorsloc
      integer, pointer :: iitors(:,:),nbimptor(:)
      real*8, pointer :: itors1(:,:),itors2(:,:),itors3(:,:)
      common /imptor/ itors1,itors2,itors3,nitors,iitors,nbimptor
