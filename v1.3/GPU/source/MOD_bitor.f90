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
!     nbitor  total number of bitorsions in the system
!     winnbitor    window object corresponding to nbitor
!     nbitorloc  local number of bitorsions in the system
!     ibitor  numbers of the atoms in each bitorsion
!     winibitor    window object corresponding to ibitor
!     nbitor_pe  total number of bitorsions per process element in the system
!     nbbitors  numbers of bitorsions before each angle in the global index
!
!
module bitor
   implicit none
   integer nbitor,nbitorloc,nbitor_pe
   integer, pointer :: ibitor(:,:), nbbitors(:)
   integer :: winibitor,winnbbitors
end
