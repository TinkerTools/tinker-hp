c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module tors  --  torsional angles within the current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     tors1   1-fold amplitude and phase for each torsional angle
c     tors2   2-fold amplitude and phase for each torsional angle
c     tors3   3-fold amplitude and phase for each torsional angle
c     tors4   4-fold amplitude and phase for each torsional angle
c     tors5   5-fold amplitude and phase for each torsional angle
c     tors6   6-fold amplitude and phase for each torsional angle
c     ntors   total number of torsional angles in the system
c     ntorsloc   local number of torsional angles in the system
c     nbtors   number of torsional angles before each atom
c     itors   numbers of the atoms in each torsional angle
c
c
      module tors
      implicit none
      integer ntors,ntorsloc
      integer, pointer :: nbtors(:)
      integer, pointer :: itors(:,:)
      real*8, pointer :: tors1(:,:),tors2(:,:),tors3(:,:)
      real*8, pointer :: tors4(:,:),tors5(:,:),tors6(:,:)
      save 
      end
