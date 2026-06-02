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
!     nuse   total number of active atoms in energy calculation
!     iuse   numbers of the atoms active in energy calculation
!     winiuse window object corresponding to iuse
!     use    true if an atom is active, false if inactive
!     useAll true if all atoms are being used
!     winuse window object corresponding to use
!
!
module usage
   implicit none
   integer nuse
   integer, pointer :: iuse(:)
   logical, pointer :: use(:)
   integer winuse,winiuse
   logical useAll
end
