!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module atmtyp  --  atomic properties for each current atom  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     story     !<descriptive type for each atom in system
!
!
module atmtyp
   implicit none
   integer :: winclass !<window object corresponding to class
   integer :: winatomic !<window object corresponding to atomic
   integer :: winvalence !<window object corresponding to valence
   integer :: winmass !<window object corresponding to mass
   integer :: winstory !<window object corresponding to story
   integer, allocatable :: tag(:) !<integer atom labels from input coordinates file
   integer, pointer :: class(:) !<atom class number for each atom in the system
   integer, pointer :: atomic(:) !<atomic number for each atom in the system
   integer, pointer :: valence(:) !<valence number for each atom in the system
   real*8, pointer :: mass(:) !<atomic weight for each atom in the system
   character*3, allocatable :: name(:) !<atom name for each atom in the system
   character*24, pointer :: story(:) !<descriptive type for each atom in system
   save
end
