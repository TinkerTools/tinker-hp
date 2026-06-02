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
!     mass      atomic weight for each atom in the system
!     tag       integer atom labels from input coordinates file
!     class     atom class number for each atom in the system
!     atomic    atomic number for each atom in the system
!     valence   valence number for each atom in the system
!     name      atom name for each atom in the system
!     story     descriptive type for each atom in system
!     win*      window object corresponding to *
!
!
#include "tinker_macro.h"
module atmtyp
   implicit none
   integer,   pointer :: tag(:),class(:)
   integer,   pointer :: atomic(:),valence(:)
   real(r_p), pointer :: mass(:)
   character*3,  pointer :: name(:)
   character*24, pointer :: story(:)

   integer wintag,winatomic,winvalence,winmass&
      &,winstory,winclass
end
