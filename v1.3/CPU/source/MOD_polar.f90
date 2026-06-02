!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module polar  --  polarizabilities and induced dipole moments  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!
module polar
   implicit none
   integer npolar !<total number of polarizable sites in the system
   integer :: winpolarity !<window object corresponding to polariy
   integer :: winpolarity_orig !<window object corresponding to polariy_orig
   integer :: winthole !<window object corresponding to thole
   integer :: winpdamp !<window object corresponding to pdamp
   integer :: wintholed !<window object corresponding to tholed
   integer :: winjpolar !<window object corresponding to jpolar
   integer, pointer :: jpolar(:) !<index into polarization parameter matrix for each atom
   real*8, parameter :: tinypol = 1e-5 !<negligeable polarisability to allow convergence
   real*8, pointer :: polarity(:) !<dipole polarizability for each multipole site (Ang**3)
   real*8, pointer :: polarity_orig(:) !<original dipole polarizability for each multipole site (Ang**3) (lambdadyn)
   real*8, pointer :: thole(:) !<Thole polarizability damping value for each site
   real*8, pointer :: pdamp(:) !<value of polarizability scale factor for each site
   real*8, pointer :: tholed(:) !<Thole direct polarization damping value for each atom
   real*8, allocatable :: uind(:,:) !<induced dipole components at each multipole site
   real*8, allocatable :: uinp(:,:) !<induced dipoles in field used for energy interactions
   real*8, allocatable :: thlval(:,:) !<Thole damping parameter value for each atom type pair
   real*8, allocatable :: thdval(:,:) !<alternate Thole direct damping value for atom type pair
   save
end
