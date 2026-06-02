!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module repel  --  Pauli repulsion for current structure  ##
!     ##                                                           ##
!     ###############################################################
!
!
!
module repel
   implicit none
   integer :: nrep !<total number of repulsion sites in the system
   integer :: nreploc !<local number of repulsion sites in the system
   real*8, pointer :: sizpr(:) !<Pauli repulsion size parameter value at each site
   real*8, pointer :: dmppr(:) !<Pauli repulsion alpha damping value at each site
   real*8, pointer :: elepr(:) !<Pauli repulsion valence electrons at each site
   integer :: winsizpr !<window object corresponding to sizepr
   integer :: windmppr !<window object corresponding to damppr
   integer :: winelepr !<window object corresponding to elepr
   save
end
