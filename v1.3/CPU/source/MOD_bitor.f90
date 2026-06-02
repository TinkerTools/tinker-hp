!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module bitor  --  bitorsions within the current structure  ##
!     ##                                                             ##
!     #################################################################
!
!
module bitor
   implicit none
   integer :: nbitor !<total number of bitorsions in the system
   integer :: nbitorloc !<local number of bitorsions in the system
   integer :: winibitor !<window object corresponding to nbitor
   integer :: winnbbitors !<window object corresponding to ibitor
   integer, pointer :: ibitor(:,:) !<numbers of the atoms in each bitorsion
   integer, pointer :: nbbitors(:) !<numbers of bitorsions before each angle in the global index
   save
end
