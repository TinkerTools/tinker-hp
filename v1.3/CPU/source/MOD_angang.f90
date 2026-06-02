!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module angang  --  angle-angle terms in current structure  ##
!     ##                                                             ##
!     #################################################################
!
!
module angang
   implicit none
   integer :: nangang        !<total number of angle-angle interactions
   integer :: nangangloc     !<local number of angle-angle interactions
   integer :: winnbangang    !<window object corresponding to nbangang
   integer :: winiaa         !<window object corresponding to kaa
   integer :: winkaa         !<window object corresponding to kaa
   integer, pointer :: nbangang(:)!<total number of angle-angle interactions before each atom in the global index
   integer, pointer :: iaa(:,:)!<angle numbers used in each angle-angle term
   real*8, pointer ::  kaa(:)!<force constant for angle-angle cross terms
   save
end
