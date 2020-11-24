c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  module pitors  --  pi-orbital torsions in the current structure  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     npitors   total number of pi-orbital torsional interactions
c     npitorsloc   local number of pi-orbital torsional interactions
c     nbpitors  number of pi-orbital torsional interactions before each atom
c     winnbpitors    window object corresponding to nbpitors
c     kpit      2-fold pi-orbital torsional force constants
c     winkpit    window object corresponding to kpit
c     ipit      numbers of the atoms in each pi-orbital torsion
c     winipit    window object corresponding to ipit
c
c
      module pitors
      implicit none
      integer npitors,npitorsloc
      integer, pointer :: ipit(:,:), nbpitors(:)
      real*8, pointer :: kpit(:)
      integer :: winipit,winnbpitors,winkpit
      save 
      end
