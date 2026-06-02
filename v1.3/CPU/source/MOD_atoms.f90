!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module atoms  --  number, position and type of current atoms  ##
!     ##                                                                ##
!     ####################################################################
!
!
module atoms
   implicit none
   integer :: n !<total number of atoms in the current system
   logical :: pbcunwrap
   integer, allocatable :: type(:) !<atom type number for each atom in the system
   integer, allocatable :: pbcwrapindex(:,:) !<array containing the indexes (n1,n2,n3) of the "real" cell of each atom
   real*8, allocatable :: x(:)  !<current x-coordinate for each atom in the system
   real*8, allocatable :: y(:)  !<current y-coordinate for each atom in the system
   real*8, allocatable :: z(:)  !<current z-coordinate for each atom in the system
   real*8, allocatable :: xwrite(:) !<x-coordinate to be printed (potentially unwrapped) for each atom in the system
   real*8, allocatable :: ywrite(:) !<y-coordinate to be printed (potentially unwrapped) for each atom in the system
   real*8, allocatable :: zwrite(:) !<z-coordinate to be printed (potentially unwrapped) for each atom in the system
   real*8, allocatable :: xold(:) !<last x-coordinate for each atom in the system
   real*8, allocatable :: yold(:) !<last y-coordinate for each atom in the system
   real*8, allocatable :: zold(:) !<last z-coordinate for each atom in the system
   save
end
