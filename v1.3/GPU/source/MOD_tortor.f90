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
!     ntortor   total number of torsion-torsion interactions
!     ntortorloc   local number of torsion-torsion interactions
!     itt       atoms and parameter indices for torsion-torsion
!     winitt    window object corresponding to itt
!     nbtortor   number of atoms before each torsion-torsion
!     winnbtortor    window object corresponding to nbtortor
!
!
module tortor
   implicit none
   integer ntortor,ntortorloc
   integer, pointer :: itt(:,:),nbtortor(:)
   integer :: winitt,winnbtortor
end
