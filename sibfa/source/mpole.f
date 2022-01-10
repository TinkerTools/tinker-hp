c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module mpole  --  multipole components for current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c
c     npole     total number of multipole sites in the system
c     npoleloc  local number of multipole sites in the system
c     npolebloc  local+neighbor number of multipole sites in the system
c     npolelocnl  localnl number of multipole sites in the system
c     npolerecloc  local number of reciprocal multipole sites in the system
c
c     poleloc   global-local correspondance for multipoles
c     polelocnl global-localnl correspondace for multipoles
c     polerecloc   global-local correspondance for reciprocal multipoles
c     nbpole    number of multipoles before each atom
c
c     ipole     number of the atom for each multipole site
c     polsiz    number of multipole components at each atom
c     pollist   multipole site for each atom (0=no multipole)
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     polaxe    local axis type for each multipole site
c
c     alphapen  alpha parameters for sibfa like charge penetration
c     betapen   beta parameters for sibfa like charge penetration
c     gammapen  gamma parameters for sibfa like charge penetration
c     vdwpen    vdw radius parameters for sibfa like charge penetration
c
c     gamma1pen  gamma parameters for original sibfa like charge penetration
c     deltapen  delta parameters for  original sibfa like charge penetration
c     khipen   khi parameters for  original sibfa like charge penetration
c     omegapen  omega parameters for original sibfa like charge penetration
c
c
      module mpole
      implicit none
      integer valemtp(18)
      integer maxpole
      parameter (maxpole=13)
      integer npole,npoleloc,npolebloc,npolerecloc,npolelocnl
      integer, pointer :: ipole(:), polsiz(:), pollist(:)
      integer, allocatable :: poleloc(:),polelocnl(:),polerecloc(:)
      integer, pointer :: nbpole(:)
      integer, pointer :: zaxis(:),xaxis(:),yaxis(:)
      real*8, allocatable :: rpole(:,:)
      real*8, pointer :: pole(:,:)
      real*8, pointer :: alphapen(:),betapen(:),gammapen(:)
      real*8 gamma1pen,deltapen,khipen,omegapen
      real*8, pointer :: vdwpen(:)
      real*8 vmxx,vmyy,vmzz
      real*8 vmxy,vmxz,vmyz
      character*8, pointer ::  polaxe(:)
      data valemtp/ 1,2,1,2,3,4,5,6,7,8,
     $  1,2,3,4,5,6,7,8/
      save
      end
