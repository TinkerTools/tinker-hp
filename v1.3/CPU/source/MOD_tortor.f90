!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module tortor  --  torsion-torsions in the current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
module tortor
   implicit none
   integer :: ntortor !<total number of torsion-torsion interactions
   integer :: ntortorloc !<local number of torsion-torsion interactions
   integer :: winitt !<window object corresponding to itt
   integer :: winnbtortor !<window object corresponding to nbtortor
   integer, pointer :: itt(:,:) !<atoms and parameter indices for torsion-torsion
   integer, pointer :: nbtortor(:) !<number of atoms before each torsion-torsion
   save
end
