!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module disp  --  damped dispersion for current structure  ##
!     ##                                                            ##
!     ################################################################
!
!
module disp
   implicit none
   integer :: ndisp !<total number of dispersion sites in the system
   integer :: ndisploc !<local number of dispersion sites in the system
   integer :: ndispbloc !<local+neighbor number of dispersion sites in the system
   integer :: ndisplocnl !<localnl number of dispersion sites in the system
   integer :: ndisprecloc !<local reciprocal number of dispersion sites in the system
   integer :: winidisp !<window object corresponding to idisp
   integer :: wincsix !<window object corresponding to csix
   integer :: winadisp !<window object corresponding to adisp
   integer :: windisplist !<window object corresponding to displist
   integer :: winnbdisp !<window object corresponding to nbdisp
   integer, pointer :: idisp(:) !<number of the atom for each dispersion site
   integer, pointer :: displist(:) !<number of the dispersion site for each atom
   integer, pointer :: nbdisp(:) !<number of dispersion sites before each atom
   integer, allocatable :: displocnl(:) !<glob-locnl dispersion correspondance
   integer, allocatable :: disprecloc(:) !<global-local reciprocal dispersion correspondance
   real*8 ::  csixpr !<pairwise sum of C6 dispersion coefficients
   real*8, pointer :: csix(:) !<C6 dispersion coefficient value at each site
   real*8, pointer :: adisp(:) !<alpha dispersion damping value at each site
   save
end
