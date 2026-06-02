!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module usage  --  atoms active during energy computation  ##
!     ##                                                            ##
!     ################################################################
!
!
!
module usage
   implicit none
   integer :: nuse !<total number of active atoms in energy calculation
   integer :: winuse !<window object corresponding to use
   integer :: winiuse !<window object corresponding to iuse
   integer, pointer :: iuse(:) !<numbers of the atoms active in energy calculation
   logical, pointer :: use(:) !<true if an atom is active, false if inactive
   save
end
