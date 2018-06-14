c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  angle.i  --  bond angles within the current structure  ##
c     ##                                                         ##
c     #############################################################
c
c
c     ak       harmonic angle force constant (kcal/mole/rad**2)
c     anat     ideal bond angle or phase shift angle (degrees)
c     afld     periodicity for Fourier bond angle term
c     nangle   total number of bond angles in the system
c     iang     numbers of the atoms in each bond angle
c
c
      integer nangle,nangleloc
      integer, pointer :: iang(:,:)
      real*8, pointer ::  ak(:),anat(:),afld(:)
c      common /angle/ ak(maxang),anat(maxang),afld(maxang),nangle,
c     &               iang(4,maxang)
      common /angle/ ak,anat,afld,nangle,iang,nangleloc
