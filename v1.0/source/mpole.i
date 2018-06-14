c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  mpole.i  --  multipole components for current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     polsiz    number of multipole components at each atom
c     pollist   multipole site for each atom (0=no multipole)
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     polaxe    local axis type for each multipole site
c
c
      integer maxpole
      parameter (maxpole=13)
      integer npole,npoleloc,npolebloc,npolerecloc,npolelocnl
      integer, pointer :: ipole(:), polsiz(:), pollist(:),poleloc(:)
      integer, pointer :: polelocnl(:)
      integer, pointer :: nbpole(:),polerecloc(:)
      integer, pointer :: zaxis(:),xaxis(:),yaxis(:)
      real*8, pointer :: pole(:,:),rpole(:,:)
      real*8, pointer :: alphapen(:),betapen(:),gammapen(:)
      character*8, pointer ::  polaxe(:)
      common /mpole/ pole,rpole,npole,npolelocnl,
     &               ipole,polsiz,pollist,poleloc,
     &               zaxis,xaxis,yaxis,
     &               polaxe,nbpole,npoleloc,
     &               npolebloc,npolerecloc,polerecloc,polelocnl,
     &               alphapen,betapen,gammapen
