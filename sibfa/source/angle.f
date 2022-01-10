c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module angle  --  bond angles within the current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     ak       harmonic angle force constant (kcal/mole/rad**2)
c     anat     ideal bond angle or phase shift angle (degrees)
c     afld     periodicity for Fourier bond angle term
c     nangle   total number of bond angles in the system
c     iang     numbers of the atoms in each bond angle
c     nangleloc numbers of bond angles in the local domain
c     angleloc correspondance between global and local bond angles
c
c
      module angle
      implicit none
      integer nangle,nangleloc
      integer, allocatable :: angleloc(:)
      integer, pointer :: iang(:,:)
      real*8, pointer ::  ak(:),anat(:),afld(:)
      save
      end
