!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module mpole  --  multipole components for current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
module mpole
   implicit none
   integer :: maxpole !<max components (monopole=1,dipole=4,quadrupole=13)
   parameter (maxpole=13)
   integer :: npole !<total number of multipole sites in the system
   integer :: npoleloc !<local number of multipole sites in the system
   integer :: npolebloc !<local+neighbor number of multipole sites in the system
   integer :: npolerecloc !<local number of reciprocal multipole sites in the system
   integer :: npolelocnl !<localnl number of multipole sites in the system
   integer :: winipole !<window object corresponding to ipole
   integer :: winpolsiz !<window object corresponding to polsiz
   integer :: winpollist !<window object corresponding to pollist
   integer :: winnbpole !<window object corresponding to nbpole
   integer :: winmono0 !<window object corresponding to mono0
   integer :: winpole !<window object corresponding to pole
   integer :: winpole_orig !<window object corresponding to pole_orig
   integer :: winalphapen !<window object corresponding to alphapen
   integer :: winbetapen !<window object corresponding to betapen
   integer :: wingammapen !<window object corresponding to gammapen
   integer :: winpolaxe !<window object corresponding to polaxe
   integer, pointer :: ipole(:) !<number of the atom for each multipole site
   integer, pointer :: polsiz(:) !<number of multipole components at each atom
   integer, pointer :: pollist(:) !<multipole site for each atom (0=no multipole)
   integer, pointer :: nbpole(:) !<number of multipoles before each atom
   integer, allocatable :: poleloc(:) !<global-local correspondance for multipoles
   integer, allocatable :: polelocnl(:) !<global-localnl correspondance for multipoles
   integer, allocatable :: polerecloc(:) !<global-local correspondance for reciprocal multipoles
   integer, allocatable :: zaxis(:) !<number of the z-axis defining atom for each site
   integer, allocatable :: xaxis(:) !<number of the x-axis defining atom for each site
   integer, allocatable :: yaxis(:) !<number of the y-axis defining atom for each site
   real*8 :: vmxx !<temporary virial component
   real*8 :: vmyy !<temporary virial component
   real*8 :: vmzz !<temporary virial component
   real*8 :: vmxy !<temporary virial component
   real*8 :: vmxz !<temporary virial component
   real*8 :: vmyz !<temporary virial component
   real*8, pointer :: mono0(:) !<original atomic monopole values for charge flux
   real*8, pointer :: pole(:,:) !<multipole values for each site in the local frame
   real*8, pointer :: pole_orig(:,:) !<original multipole values for each site in the local frame (lambdadyn)
   real*8, pointer :: alphapen(:) !<alpha parameters for sibfa like charge penetration
   real*8, pointer :: betapen(:) !<beta parameters for sibfa like charge penetration
   real*8, pointer :: gammapen(:) !<gamma parameters for sibfa like charge penetration
   real*8, allocatable :: rpole(:,:) !<multipoles rotated to the global coordinate system
   character*8, pointer ::  polaxe(:) !<local axis type for each multipole site
   save
end
