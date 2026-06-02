!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module bond  --  covalent bonds in the current structure  ##
!     ##                                                            ##
!     ################################################################
!
!
module bond
   implicit none
   integer :: nbond !<total number of bond stretches in the system
   integer :: nbondloc !<local number of bond stretches in the system
   integer :: winibnd !<window object corresponding to ibnd
   integer :: winbk !<window object corresponding to bk
   integer :: winbl !<window object corresponding to bl
   integer :: winba !<window object corresponding to ba
   integer, pointer :: ibnd(:,:) !<numbers of the atoms in each bond stretch
   real*8, pointer ::  bk(:) !<bond stretch force constants (kcal/mole/Ang**2)
   real*8, pointer :: bl(:) !<ideal bond length values in Angstroms
   real*8, pointer :: ba(:) !<a parameter for morse potential
   save
end
