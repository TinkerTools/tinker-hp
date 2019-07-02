c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module moldyn  --  velocity and acceleration on MD trajectory  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     v       current velocity of each atom along the x,y,z-axes
c     a       current acceleration of each atom along x,y,z-axes
c     aalt    alternate acceleration of each atom along x,y,z-axes
c
c
      module moldyn
      implicit none
      real*8,allocatable :: v(:,:),a(:,:),aalt(:,:)
      save
      end
