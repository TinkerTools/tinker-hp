c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  pitors.i  --  pi-orbital torsions in the current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     kpit      2-fold pi-orbital torsional force constants
c     npitors   total number of pi-orbital torsional interactions
c     ipit      numbers of the atoms in each pi-orbital torsion
c
c
      integer npitors,npitorsloc
      integer, pointer :: ipit(:,:), nbpitors(:)
      real*8, pointer :: kpit(:)
      common /pitors/ kpit,npitors,ipit,nbpitors,npitorsloc
