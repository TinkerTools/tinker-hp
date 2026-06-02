!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module molcul  --  individual molecules within current system  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     molmass   molecular weight for each molecule in the system
!     winmolmass    window object corresponding to molmass
!     totmass   total weight of all the molecules in the system
!     nmol      total number of separate molecules in the system
!     nmoleloc  local number of separate molecules in the system
!     kmol      contiguous list of the atoms in each molecule
!     winkmol    window object corresponding to kmol
!     imol      first and last atom of each molecule in the list
!     winimol    window object corresponding to imol
!     molcule   number of the molecule to which each atom belongs
!     winmolcule    window object corresponding to molcule
!
!
#include "tinker_macro.h"
module molcul
   implicit none
   integer nmol,nmoleloc
   integer :: winmolcule,winkmol,winimol
   integer :: winmolmass
   integer, pointer :: molcule(:),kmol(:),imol(:,:)
   real(r_p) totmass
   real(t_p), pointer :: molmass(:)
end
