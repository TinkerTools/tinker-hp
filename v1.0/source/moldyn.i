c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  moldyn.i  --  velocity and acceleration on MD trajectory  ##
c     ##                                                            ##
c     ################################################################
c
c
c     v       current velocity of each atom along the x,y,z-axes
c     a       current acceleration of each atom along x,y,z-axes
c     aalt    alternate acceleration of each atom along x,y,z-axes
c
c
      real*8,pointer :: v(:,:),a(:,:),aalt(:,:)
      common /moldyn/ v,a,aalt
